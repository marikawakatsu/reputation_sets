
number = 7

# function is_even(number::Int64)
#     while number > 0
#         number -= 2
#     end
#     number == 0 ? (return true) : (return false)
# end
#
# is_even(number)

# function is_even(number::Int64)
#     if number > 0
#         if number == 1
#             return false
#         else
#             return is_even(number - 2)
#         end
#     else
#         if number == 0
#             return true
#         else
#             return is_even(-number)
#         end
#     end
# end

function is_even(number::Int64)
    number > 0 ? (number == 1 ? (return false) : (return is_even(number - 2))) : (number == 0 ? (return true) : (return is_even(-number)))
end

number = 519051

is_even(number)
