import pronto
import pandas as pd #we would use it to count the frequencies
import subprocess #to run some linux instructions

#-----------------output directories------
output_directory='/root/HPO_project/HPO_new/New_analysis/orion_queries'
input_directory='/root/HPO_project/resources'
query_file='/root/HPO_project/HPO_new/log_pcf_ds_20210608_20240826.txt'
# Load the HPO ontology
hpo = pronto.Ontology(f'{input_directory}/hp.obo')

#---------------------- FUNCTION TO GENERATE THE FREQUENCY FILES WITH AN HPO LIST--------------------------------------
def frequency_files(input_list, output_file_name):
    # Convert the list to a Series
    serie = pd.Series(input_list)
    # Count the frequency of each value
    HPO_frequency = serie.value_counts()
    # Convert the Series to a DataFrame with two columns (HPO and Count)
    HPO_frequency_df = HPO_frequency.reset_index()  # Helps to eliminate misunderstandings with headers
    # Rename the columns
    HPO_frequency_df.columns = ['HPO', 'Count']
    
    print("HPO Frequency DataFrame:")
    print(HPO_frequency_df)
    
    # Loop through the DataFrame rows
    for index, row in HPO_frequency_df.iterrows():
        term = row['HPO']  # Obtain the value of the 'HPO' column for each row

        try:
            hpo_code = hpo[term]  # Look up the term in the HPO ontology
            hpo_name = hpo_code.name  # Get the name of the HPO term
            
            # Calculate the frequency of each HPO
            count_sum = HPO_frequency_df['Count'].sum()  # Total number of HPO
            frequency = row['Count'] / count_sum  # Calculate the frequency
            
            # Save the HPO name and frequency to new columns
            HPO_frequency_df.at[index, 'HPO_Name'] = hpo_name
            HPO_frequency_df.at[index, 'HPO_frequency'] = frequency

        except KeyError:
            print(f"Warning: HPO term '{term}' not found in the HPO ontology.")
            HPO_frequency_df.at[index, 'HPO_Name'] = 'Unknown'
            HPO_frequency_df.at[index, 'HPO_frequency'] = 0  # or some other placeholder value
    
    # Save the resulting DataFrame to a TSV file
    HPO_frequency_df.to_csv(f'{output_directory}/{output_file_name}', sep='\t', index=False)



# Read the CSV file into a DataFrame
df = pd.read_csv(f'{output_directory}/hpo_dedup_merged.csv', sep=',', header=None)
# Flatten the DataFrame and split each cell by comma
queries_list = [item for sublist in df.values.flatten() for item in str(sublist).split(',')]


#----GENERATION OF INDIVIDUAL HPO FREQUENCY FILE----------------------------
frequency_files(queries_list, 'IndividualHPO_frequency.tsv')


#----OBTAINING FREQUENCY FOR GROUPS COMBINATION------------------------------------------
df_queries2 = pd.read_csv(f'{output_directory}/hpo_dedup_merged.csv', header=None)
# Primero ordenamos cada consulta para que la combinación de grupos se pueda contar incluso si estaban en un orden diferente
df_queries2[0] = df_queries2[0].apply(lambda x: sorted(x.split(',')))
# Guardar el DataFrame en un nuevo archivo CSV
df_queries2.to_csv(f'{output_directory}/queries_sorted.csv', index=False, header=False)

# Comando de Linux para eliminar comillas no deseadas del archivo
df_queries2 = pd.read_csv(f'{output_directory}/queries_sorted.csv', sep=',', header=None)
sed_command = f"sed -i 's/\\[//g; s/\\]//g' {output_directory}/queries_sorted.csv"
sed_command2 = f"sed -i \"s/'//g; s/ //g\" {output_directory}/queries_sorted.csv"

# Ejecutar el comando usando subprocess
subprocess.run(sed_command, shell=True)
subprocess.run(sed_command2, shell=True)

df_queries2 = pd.read_csv(f'{output_directory}/queries_sorted.csv', header=None)
# Luego hacemos una lista con cada combinación
list_of_values = df_queries2.values.flatten().tolist()
# Convertir la lista a una Serie
serie = pd.Series(list_of_values)
# Contar la frecuencia de cada valor
HPO_frequency = serie.value_counts()
print(HPO_frequency)
HPO_frequency_df = HPO_frequency.reset_index()  # Ayuda a eliminar malentendidos con los encabezados
# Renombrar las columnas
HPO_frequency_df.columns = ['HPO', 'Count']
HPO_frequency_df.to_csv(f'{output_directory}/combination_HPO_freq', sep='\t', index=False)

# Function to get names of HPO terms
def get_names(hpo_terms):
    terms = hpo_terms.split(',')  # Split the HPO terms
    names = [hpo[term].name if term in hpo else "Unknown" for term in terms]  # Get names
    return ', '.join(names)  # Join names with commas

# Add "names" and "frequency" columns
HPO_frequency_df['names'] = HPO_frequency_df['HPO'].apply(get_names)
HPO_frequency_df['frequency'] = HPO_frequency_df['Count'] / HPO_frequency_df['Count'].sum()

# Save the DataFrame to a new file
HPO_frequency_df.to_csv(f'{output_directory}/combination_HPO_freq.tsv', sep='\t', index=False)


########## PHENOTYPIC ABNORMALITY-------------------------------------
#----------------------------------------------#
#----------------------------------------------#
def get_descendants_iterative(term):
    """Gets all descendants of a term iteratively."""
    to_visit = list(term.subclasses())  # Start with direct children
    descendants = set(to_visit)  # A set to store descendants without duplicates

    while to_visit:
        current = to_visit.pop(0)  # Take the first term from the list
        children = current.subclasses()  # Get its direct children
        for child in children:
            if child not in descendants:
                descendants.add(child)
                to_visit.append(child)  # Add new terms to explore

    return descendants

phenotype_term = hpo['HP:0000118']  # HPO term for 'phenotypic abnormality'
phenotypes = phenotype_term.subclasses(distance=1)  # Get direct subclasses at level 1

branch_count = {}  # Initialize the dictionary to store the results

for phenotype in phenotypes:
    # Get all descendants of the current term using the iterative function
    all_descendants = get_descendants_iterative(phenotype)
    # Store the number of descendants for each term in the dictionary
    branch_count[phenotype.id] = len(all_descendants)


print(branch_count)


#------------------------to calculate the frequency of terms in a specific HPO

def get_descendants_iterative(term):
    """Gets all descendants of a term iteratively."""
    to_visit = list(term.subclasses())  # Start with direct children
    descendants = set(to_visit)  # A set to store descendants without duplicates

    while to_visit:
        current = to_visit.pop(0)  # Take the first term from the list
        children = current.subclasses()  # Get its direct children
        for child in children:
            if child not in descendants:
                descendants.add(child)
                to_visit.append(child)  # Add new terms to explore

    return descendants

phenotype_term = hpo['HP:0000118']  # HPO term for 'phenotypic abnormality'
phenotypes = phenotype_term.subclasses(distance=1)  # Get direct subclasses at level 1

branch_count = {}  # Initialize the dictionary to store the results

for phenotype in phenotypes:
    # Get all descendants of the current term using the iterative function
    all_descendants = get_descendants_iterative(phenotype)
    # Store the number of descendants for each term in the dictionary
    branch_count[phenotype.id] = len(all_descendants)

#    print(f"ID: {phenotype_id}, Total descendants: {branch_count[phenotype_id]-1}")   

#------function to count the number of terms a branch has (to count the frequency of branch in the correct ratio)-----------------
def frequency_root(input_list, output_file_name):
    # Convert the list to a Series
    serie = pd.Series(input_list)
    # Count the frequency of each value
    HPO_frequency = serie.value_counts()
    # Convert the Series to a DataFrame with two columns (HPO and Count)
    HPO_frequency_df = HPO_frequency.reset_index()  # Helps to eliminate misunderstandings with headers
    # Rename the columns
    HPO_frequency_df.columns = ['HPO', 'Count']
    
    # Loop through the DataFrame rows
    for index, row in HPO_frequency_df.iterrows():
        term = row['HPO']  # Obtain the value of the 'HPO' column for each row

        try:
            hpo_code = hpo[term]  # Look up the term in the HPO ontology
            hpo_name = hpo_code.name  # Get the name of the HPO term
            
            # Calculate the frequency of each HPO 
            count_sum =  branch_count[hpo_code.id] # Total number of HPO
            frequency = row['Count'] / count_sum  # Calculate the frequency
            
            # Save the HPO name and frequency to new columns
            HPO_frequency_df.at[index, 'HPO_Name'] = hpo_name
            HPO_frequency_df.at[index, 'HPO_frequency'] = frequency

        except KeyError:
            print(f"Warning: HPO term '{term}' not found in the HPO ontology.")
            HPO_frequency_df.at[index, 'HPO_Name'] = 'Unknown'
            HPO_frequency_df.at[index, 'HPO_frequency'] = 0  # or some other placeholder value
    
    # Save the resulting DataFrame to a TSV file
    HPO_frequency_df.to_csv(f'{output_directory}/{output_file_name}', sep='\t', index=False)


#----------------Phenotypic abnormality dictionary-------------------------------------
phenotype_term= hpo['HP:0000118'] #this is the HPO for the category 'phenotypic_abnormality'
phenotypes=phenotype_term.subclasses(distance=1) #we obtain the next level of subclasses, which are the phenotypes
phenotypic_abnormality={} #we initialize our dictionary to save each phenotypic abnormality with its HPO
for phenotype in phenotypes: #filling the dictionary
    phenotypic_abnormality[phenotype.id]=phenotype.name
del phenotypic_abnormality['HP:0000118'] #the subclass comand also saves the current level, so it was eliminated


#-----------------------Found the phenotypic abnormality of each HPO------------------------------------
#serie = pd.Series(queries_list)
#we  already have saved all the HPO's in a list so we just need to translate the list into its corresponding phenotypic abnormality
#we define our list of phenotypic abnormality
phe_abnor=[]
for value in queries_list:
    if value in hpo:
        term = hpo[value]
        superclasses = term.superclasses(distance=30) #obtains all previous levels of branches (set at 20)
        # Iterate over the superclasses to find the corresponding phenotypic abnormality to the desired HPO
        for superclass in superclasses:
            if superclass.id in phenotypic_abnormality:
                #print(f'The phenotypic abnormality for {term.name} is {phenotypic_abnormality[superclass.id]}')
                phe_abnor.append(superclass.id)
            else:
                continue
    else:print(f'HPO {value} does not exist in the dataset, skipping.')


#----obtaining the frequency file for phenotypic_abnormality-----------------
# Ejecutar la función predefinida "frequency files"
frequency_root(phe_abnor, 'phenotypic_abnormality_frequency.tsv')
