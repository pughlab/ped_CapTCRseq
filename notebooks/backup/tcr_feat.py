import pandas as pd
import numpy as np

# Define amino acid letters
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Create a DataFrame with 5 columns and random values
df = pd.DataFrame(np.random.rand(20, 5), 
                  columns=['I', 'II', 'III', 'IV', 'V'],
                  index=amino_acids)

# Display the table
print(df)
