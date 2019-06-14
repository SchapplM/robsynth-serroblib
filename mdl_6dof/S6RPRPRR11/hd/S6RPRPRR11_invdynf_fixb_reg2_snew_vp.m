% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 20:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRPRR11_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR11_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_invdynf_fixb_reg2_snew_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:26:27
% EndTime: 2019-05-05 20:26:51
% DurationCPUTime: 25.50s
% Computational Cost: add. (296604->423), mult. (944235->689), div. (0->0), fcn. (817045->16), ass. (0->368)
t3024 = sin(pkin(12));
t3025 = sin(pkin(7));
t3029 = cos(pkin(7));
t3026 = sin(pkin(6));
t3028 = cos(pkin(12));
t3030 = cos(pkin(6));
t3034 = sin(qJ(1));
t3038 = cos(qJ(1));
t3016 = t3034 * g(1) - t3038 * g(2);
t3039 = qJD(1) ^ 2;
t3117 = t3026 * t3039;
t3051 = -qJ(2) * t3117 - t3016;
t3132 = pkin(9) * t3024;
t3055 = -pkin(2) * t3028 - t3025 * t3132;
t3101 = -t3030 * g(3) + qJDD(2);
t3119 = t3025 * t3030;
t3019 = t3024 ^ 2;
t3021 = t3028 ^ 2;
t3144 = -t3019 - t3021;
t3042 = ((-pkin(1) + t3055) * qJDD(1) + t3051) * t3026 + (t3024 * t3030 * pkin(2) + (t3026 * t3029 * t3144 - t3028 * t3119) * pkin(9)) * t3117 + t3101;
t3120 = t3024 * t3026;
t3109 = t3025 * t3120;
t3118 = t3026 * t3028;
t3116 = t3028 * t3029;
t3143 = t3026 * t3116 + t3119;
t3043 = (pkin(2) * t3118 * t3120 + (t3030 * t3143 + t3109 * t3120) * pkin(9)) * qJD(1);
t3048 = qJDD(1) * pkin(1) - t3051;
t3046 = t3030 * t3048;
t3045 = (-t3026 * g(3) + t3046) * t3028;
t3017 = -t3038 * g(1) - t3034 * g(2);
t3053 = t3039 * pkin(1) - t3017;
t3111 = qJDD(1) * t3030;
t3110 = pkin(2) * t3111;
t3112 = qJDD(1) * t3026;
t3147 = t3029 * (t3110 + t3045 + ((-pkin(9) * t3029 - qJ(2)) * t3112 + t3053) * t3024 + (-0.2e1 * qJD(2) * t3120 + t3043) * qJD(1)) + t3025 * t3042;
t3102 = t3029 * t3112;
t3146 = t3025 * t3111 + t3028 * t3102;
t3033 = sin(qJ(3));
t3037 = cos(qJ(3));
t2985 = (t3033 * t3119 + (t3024 * t3037 + t3033 * t3116) * t3026) * qJD(1);
t3115 = t3029 * t3030;
t3054 = t3025 * t3118 - t3115;
t2992 = qJD(1) * t3054 - qJD(3);
t3023 = sin(pkin(13));
t3027 = cos(pkin(13));
t2972 = -t3023 * t2985 - t2992 * t3027;
t2973 = t2985 * t3027 - t2992 * t3023;
t3032 = sin(qJ(5));
t3036 = cos(qJ(5));
t2935 = -t3036 * t2972 + t2973 * t3032;
t2932 = qJD(6) + t2935;
t3145 = qJD(6) + t2932;
t3050 = -qJDD(1) * t3054 + qJDD(3);
t2937 = t2972 * t3032 + t2973 * t3036;
t2983 = (t3033 * t3120 - t3037 * t3143) * qJD(1);
t2982 = qJD(5) + t2983;
t3031 = sin(qJ(6));
t3035 = cos(qJ(6));
t2920 = t2937 * t3031 - t3035 * t2982;
t3142 = t2920 ^ 2;
t2922 = t2937 * t3035 + t2982 * t3031;
t3141 = t2922 ^ 2;
t3140 = t2932 ^ 2;
t3139 = t2935 ^ 2;
t3138 = t2937 ^ 2;
t3137 = t2972 ^ 2;
t3136 = t2973 ^ 2;
t3135 = t2982 ^ 2;
t2962 = t2983 ^ 2;
t3134 = t2985 ^ 2;
t3133 = t2992 ^ 2;
t3022 = t3030 ^ 2;
t3131 = qJD(1) * qJD(2);
t3129 = t2920 * t2922;
t3128 = t2935 * t2937;
t3127 = t2972 * t2973;
t3126 = t2973 * t2983;
t3125 = t2983 * t2972;
t3124 = t2983 * t2992;
t3123 = t2985 * t2983;
t3122 = t2985 * t2992;
t3020 = t3026 ^ 2;
t3121 = t3020 * t3039;
t3114 = qJD(5) - t2982;
t3113 = qJD(6) - t2932;
t3047 = qJ(2) * t3112 - t3053;
t2971 = -g(3) * t3120 + t3024 * t3046 + t3028 * t3047 + 0.2e1 * t3118 * t3131;
t2931 = t3143 * qJDD(1) * pkin(9) + (-t3022 * pkin(2) + (t3055 * t3118 + t3115 * t3132) * t3026) * t3039 + t2971;
t2892 = t3037 * t2931 + t3033 * t3147;
t2958 = pkin(3) * t2983 - qJ(4) * t2985;
t2876 = -pkin(3) * t3133 + qJ(4) * t3050 - t2983 * t2958 + t2892;
t2910 = -t3025 * (-g(3) * t3118 - t3024 * t3047 + t3028 * t3046 - t3102 * t3132 + t3110) + t3029 * t3042 + (0.2e1 * qJD(2) * t3109 - t3025 * t3043) * qJD(1);
t3105 = t3024 * t3112;
t3096 = t2985 * qJD(3) + t3033 * t3105 - t3037 * t3146;
t2938 = t3096 - t3122;
t2961 = -t2983 * qJD(3) + t3033 * t3146 + t3037 * t3105;
t3097 = -t2961 - t3124;
t2884 = pkin(3) * t2938 + qJ(4) * t3097 + t2910;
t2839 = 0.2e1 * qJD(4) * t2972 + t3027 * t2876 + t3023 * t2884;
t2948 = pkin(4) * t2983 - pkin(10) * t2973;
t2949 = -t2961 * t3023 + t3027 * t3050;
t2831 = -pkin(4) * t3137 + pkin(10) * t2949 - t2948 * t2983 + t2839;
t2838 = -0.2e1 * qJD(4) * t2973 - t3023 * t2876 + t3027 * t2884;
t2950 = t3027 * t2961 + t3023 * t3050;
t2915 = -t2950 + t3125;
t2916 = t3096 + t3127;
t3044 = pkin(4) * t2916 + pkin(10) * t2915 + t2838;
t2796 = t3036 * t2831 + t3032 * t3044;
t3106 = t3030 * t3117;
t3103 = t3028 * t3112;
t2795 = -t2831 * t3032 + t3036 * t3044;
t3063 = -t3032 * t2949 - t3036 * t2950;
t2897 = -qJD(5) * t2935 - t3063;
t3100 = t2982 * t2935 - t2897;
t3056 = qJDD(5) + t3096;
t3099 = -t3031 * t2897 + t3035 * t3056;
t3098 = -t3036 * t2949 + t3032 * t2950;
t3094 = t3033 * t2931 - t3037 * t3147;
t2908 = pkin(5) * t2935 - pkin(11) * t2937;
t2792 = -pkin(5) * t3135 + pkin(11) * t3056 - t2935 * t2908 + t2796;
t2875 = -t3050 * pkin(3) - t3133 * qJ(4) + t2985 * t2958 + qJDD(4) + t3094;
t2858 = -t2949 * pkin(4) - t3137 * pkin(10) + t2973 * t2948 + t2875;
t2879 = (qJD(5) + t2982) * t2937 + t3098;
t2819 = pkin(5) * t2879 + pkin(11) * t3100 + t2858;
t2776 = -t2792 * t3031 + t2819 * t3035;
t2777 = t2792 * t3035 + t2819 * t3031;
t2761 = -t2776 * t3031 + t2777 * t3035;
t2791 = -pkin(5) * t3056 - pkin(11) * t3135 + t2908 * t2937 - t2795;
t2749 = t2761 * t3032 - t2791 * t3036;
t2750 = t2761 * t3036 + t2791 * t3032;
t2738 = t2749 * t3027 + t2750 * t3023;
t2739 = -t2749 * t3023 + t2750 * t3027;
t2760 = t2776 * t3035 + t2777 * t3031;
t3092 = t2739 * t3033 - t2760 * t3037;
t2731 = -t3025 * t2738 + t3029 * t3092;
t2734 = t2739 * t3037 + t2760 * t3033;
t3093 = t2731 * t3028 + t2734 * t3024;
t2773 = t2795 * t3036 + t2796 * t3032;
t2774 = -t2795 * t3032 + t2796 * t3036;
t2754 = t2773 * t3027 + t2774 * t3023;
t2755 = -t2773 * t3023 + t2774 * t3027;
t3090 = t2755 * t3033 - t2858 * t3037;
t2743 = -t3025 * t2754 + t3029 * t3090;
t2753 = t2755 * t3037 + t2858 * t3033;
t3091 = t2743 * t3028 + t2753 * t3024;
t2853 = -t2922 * t3113 + t3099;
t3049 = -t3035 * t2897 - t3031 * t3056;
t2855 = t2920 * t3113 + t3049;
t2822 = t2853 * t3035 - t2855 * t3031;
t2877 = -t3141 - t3142;
t2807 = t2822 * t3032 - t2877 * t3036;
t2808 = t2822 * t3036 + t2877 * t3032;
t2782 = t2807 * t3027 + t2808 * t3023;
t2783 = -t2807 * t3023 + t2808 * t3027;
t2821 = t2853 * t3031 + t2855 * t3035;
t3085 = t2783 * t3033 - t2821 * t3037;
t2759 = -t3025 * t2782 + t3029 * t3085;
t2772 = t2783 * t3037 + t2821 * t3033;
t3089 = t2759 * t3028 + t2772 * t3024;
t3052 = -qJD(5) * t2937 - qJDD(6) - t3098;
t2864 = -t3052 - t3129;
t2886 = -t3140 - t3142;
t2835 = -t2864 * t3031 + t2886 * t3035;
t2852 = t2922 * t3145 - t3099;
t2813 = t2835 * t3032 - t2852 * t3036;
t2814 = t2835 * t3036 + t2852 * t3032;
t2787 = t2813 * t3027 + t2814 * t3023;
t2788 = -t2813 * t3023 + t2814 * t3027;
t2834 = t2864 * t3035 + t2886 * t3031;
t3083 = t2788 * t3033 - t2834 * t3037;
t2763 = -t3025 * t2787 + t3029 * t3083;
t2775 = t2788 * t3037 + t2834 * t3033;
t3088 = t2763 * t3028 + t2775 * t3024;
t2865 = t3052 - t3129;
t2887 = -t3140 - t3141;
t2837 = t2865 * t3035 - t2887 * t3031;
t2854 = -t2920 * t3145 - t3049;
t2817 = t2837 * t3032 - t2854 * t3036;
t2818 = t2837 * t3036 + t2854 * t3032;
t2789 = t2817 * t3027 + t2818 * t3023;
t2790 = -t2817 * t3023 + t2818 * t3027;
t2836 = t2865 * t3031 + t2887 * t3035;
t3082 = t2790 * t3033 - t2836 * t3037;
t2765 = -t3025 * t2789 + t3029 * t3082;
t2778 = t2790 * t3037 + t2836 * t3033;
t3087 = t2765 * t3028 + t2778 * t3024;
t2811 = t2838 * t3027 + t2839 * t3023;
t2812 = -t2838 * t3023 + t2839 * t3027;
t3079 = t2812 * t3033 - t2875 * t3037;
t2780 = -t3025 * t2811 + t3029 * t3079;
t2799 = t2812 * t3037 + t2875 * t3033;
t3086 = t2780 * t3028 + t2799 * t3024;
t2880 = -t2937 * t3114 - t3098;
t2882 = t2935 * t3114 + t3063;
t2846 = t2880 * t3032 + t2882 * t3036;
t2847 = t2880 * t3036 - t2882 * t3032;
t2815 = t2846 * t3027 + t2847 * t3023;
t2816 = -t2846 * t3023 + t2847 * t3027;
t2890 = -t3138 - t3139;
t3078 = t2816 * t3033 - t2890 * t3037;
t2785 = -t3025 * t2815 + t3029 * t3078;
t2806 = t2816 * t3037 + t2890 * t3033;
t3084 = t2785 * t3028 + t2806 * t3024;
t2900 = t3056 - t3128;
t2907 = -t3135 - t3139;
t2866 = t2900 * t3036 + t2907 * t3032;
t2867 = -t2900 * t3032 + t2907 * t3036;
t2832 = t2866 * t3027 + t2867 * t3023;
t2833 = -t2866 * t3023 + t2867 * t3027;
t3077 = t2833 * t3033 - t2879 * t3037;
t2794 = -t3025 * t2832 + t3029 * t3077;
t2823 = t2833 * t3037 + t2879 * t3033;
t3081 = t2794 * t3028 + t2823 * t3024;
t2901 = -t3056 - t3128;
t2911 = -t3135 - t3138;
t2871 = t2901 * t3032 + t2911 * t3036;
t2872 = t2901 * t3036 - t2911 * t3032;
t2840 = t2871 * t3027 + t2872 * t3023;
t2841 = -t2871 * t3023 + t2872 * t3027;
t3076 = t2841 * t3033 + t3037 * t3100;
t2801 = -t3025 * t2840 + t3029 * t3076;
t2824 = t2841 * t3037 - t3033 * t3100;
t3080 = t2801 * t3028 + t2824 * t3024;
t3070 = t2892 * t3033 - t3037 * t3094;
t2843 = -t3025 * t2910 + t3029 * t3070;
t2863 = t2892 * t3037 + t3033 * t3094;
t3075 = t2843 * t3028 + t2863 * t3024;
t2913 = t2949 + t3126;
t2888 = t2913 * t3023 + t2915 * t3027;
t2889 = t2913 * t3027 - t2915 * t3023;
t2918 = -t3136 - t3137;
t3071 = t2889 * t3033 - t2918 * t3037;
t2849 = -t3025 * t2888 + t3029 * t3071;
t2869 = t2889 * t3037 + t2918 * t3033;
t3074 = t2849 * t3028 + t2869 * t3024;
t2927 = -t2962 - t3137;
t2895 = t2916 * t3027 + t2927 * t3023;
t2896 = -t2916 * t3023 + t2927 * t3027;
t2912 = -t2949 + t3126;
t3068 = t2896 * t3033 - t2912 * t3037;
t2851 = -t3025 * t2895 + t3029 * t3068;
t2870 = t2896 * t3037 + t2912 * t3033;
t3073 = t2851 * t3028 + t2870 * t3024;
t2917 = -t3096 + t3127;
t2933 = -t2962 - t3136;
t2898 = t2917 * t3023 + t2933 * t3027;
t2899 = t2917 * t3027 - t2933 * t3023;
t2914 = t2950 + t3125;
t3067 = t2899 * t3033 - t2914 * t3037;
t2857 = -t3025 * t2898 + t3029 * t3067;
t2874 = t2899 * t3037 + t2914 * t3033;
t3072 = t2857 * t3028 + t2874 * t3024;
t2945 = -t3134 - t2962;
t2939 = -t3096 - t3122;
t2941 = -t2961 + t3124;
t3064 = t2939 * t3033 + t2941 * t3037;
t2894 = -t3025 * t2945 + t3029 * t3064;
t2909 = t2939 * t3037 - t2941 * t3033;
t3069 = t2894 * t3028 + t2909 * t3024;
t2952 = -t2962 - t3133;
t2954 = t3050 - t3123;
t3062 = t2952 * t3033 + t2954 * t3037;
t2903 = -t3025 * t2938 + t3029 * t3062;
t2919 = t2952 * t3037 - t2954 * t3033;
t3066 = t2903 * t3028 + t2919 * t3024;
t2953 = -t3050 - t3123;
t2955 = -t3133 - t3134;
t3061 = t2953 * t3033 + t2955 * t3037;
t2905 = t3025 * t3097 + t3029 * t3061;
t2923 = t2953 * t3037 - t2955 * t3033;
t3065 = t2905 * t3028 + t2923 * t3024;
t2970 = t3045 + ((-qJ(2) * qJDD(1) - 0.2e1 * t3131) * t3026 + t3053) * t3024;
t3060 = t2970 * t3028 + t2971 * t3024;
t3007 = t3028 * t3106;
t2994 = t3007 - t3105;
t3006 = t3024 * t3106;
t2995 = t3006 + t3103;
t3059 = t2994 * t3028 + t2995 * t3024;
t3005 = t3028 * t3024 * t3121;
t2997 = t3005 + t3111;
t3001 = (-t3020 * t3021 - t3022) * t3039;
t3058 = t2997 * t3028 + t3001 * t3024;
t2998 = t3005 - t3111;
t3000 = (-t3019 * t3020 - t3022) * t3039;
t3057 = t2998 * t3024 + t3000 * t3028;
t3014 = -qJDD(1) * t3034 - t3038 * t3039;
t3013 = qJDD(1) * t3038 - t3034 * t3039;
t2999 = t3144 * t3121;
t2996 = t3006 - t3103;
t2993 = t3007 + t3105;
t2986 = -t3026 * t3048 + t3101;
t2978 = t2998 * t3028 - t3000 * t3024;
t2977 = -t2997 * t3024 + t3001 * t3028;
t2974 = -t2994 * t3024 + t2995 * t3028;
t2966 = -t3026 * t2993 + t3030 * t3057;
t2965 = -t3026 * t2996 + t3030 * t3058;
t2964 = t3030 * t2993 + t3026 * t3057;
t2963 = t3030 * t2996 + t3026 * t3058;
t2960 = -t3026 * t2999 + t3030 * t3059;
t2959 = t3030 * t2999 + t3026 * t3059;
t2934 = -t2970 * t3024 + t2971 * t3028;
t2925 = -t3026 * t2986 + t3030 * t3060;
t2924 = t3030 * t2986 + t3026 * t3060;
t2904 = t3025 * t3061 - t3029 * t3097;
t2902 = t3029 * t2938 + t3025 * t3062;
t2893 = t3029 * t2945 + t3025 * t3064;
t2885 = -t2905 * t3024 + t2923 * t3028;
t2878 = -t2903 * t3024 + t2919 * t3028;
t2868 = -t2894 * t3024 + t2909 * t3028;
t2862 = -t3026 * t2904 + t3030 * t3065;
t2861 = t3030 * t2904 + t3026 * t3065;
t2860 = -t3026 * t2902 + t3030 * t3066;
t2859 = t3030 * t2902 + t3026 * t3066;
t2856 = t3029 * t2898 + t3025 * t3067;
t2850 = t3029 * t2895 + t3025 * t3068;
t2848 = t3029 * t2888 + t3025 * t3071;
t2845 = -t3026 * t2893 + t3030 * t3069;
t2844 = t3030 * t2893 + t3026 * t3069;
t2842 = t3029 * t2910 + t3025 * t3070;
t2829 = -t2857 * t3024 + t2874 * t3028;
t2826 = -t2851 * t3024 + t2870 * t3028;
t2825 = -t2849 * t3024 + t2869 * t3028;
t2820 = -t2843 * t3024 + t2863 * t3028;
t2810 = -t3026 * t2856 + t3030 * t3072;
t2809 = t3030 * t2856 + t3026 * t3072;
t2805 = -t3026 * t2850 + t3030 * t3073;
t2804 = t3030 * t2850 + t3026 * t3073;
t2803 = -t3026 * t2848 + t3030 * t3074;
t2802 = t3030 * t2848 + t3026 * t3074;
t2800 = t3029 * t2840 + t3025 * t3076;
t2798 = -t3026 * t2842 + t3030 * t3075;
t2797 = t3030 * t2842 + t3026 * t3075;
t2793 = t3029 * t2832 + t3025 * t3077;
t2786 = -t2801 * t3024 + t2824 * t3028;
t2784 = t3029 * t2815 + t3025 * t3078;
t2781 = -t2794 * t3024 + t2823 * t3028;
t2779 = t3029 * t2811 + t3025 * t3079;
t2771 = -t2785 * t3024 + t2806 * t3028;
t2770 = -t3026 * t2800 + t3030 * t3080;
t2769 = t3030 * t2800 + t3026 * t3080;
t2768 = -t2780 * t3024 + t2799 * t3028;
t2767 = -t3026 * t2793 + t3030 * t3081;
t2766 = t3030 * t2793 + t3026 * t3081;
t2764 = t3029 * t2789 + t3025 * t3082;
t2762 = t3029 * t2787 + t3025 * t3083;
t2758 = t3029 * t2782 + t3025 * t3085;
t2757 = -t3026 * t2784 + t3030 * t3084;
t2756 = t3030 * t2784 + t3026 * t3084;
t2752 = -t3026 * t2779 + t3030 * t3086;
t2751 = t3030 * t2779 + t3026 * t3086;
t2748 = -t2765 * t3024 + t2778 * t3028;
t2747 = -t2763 * t3024 + t2775 * t3028;
t2746 = -t2759 * t3024 + t2772 * t3028;
t2745 = -t3026 * t2764 + t3030 * t3087;
t2744 = t3030 * t2764 + t3026 * t3087;
t2742 = t3029 * t2754 + t3025 * t3090;
t2741 = -t3026 * t2762 + t3030 * t3088;
t2740 = t3030 * t2762 + t3026 * t3088;
t2737 = -t3026 * t2758 + t3030 * t3089;
t2736 = t3030 * t2758 + t3026 * t3089;
t2735 = -t2743 * t3024 + t2753 * t3028;
t2733 = -t3026 * t2742 + t3030 * t3091;
t2732 = t3030 * t2742 + t3026 * t3091;
t2730 = t3029 * t2738 + t3025 * t3092;
t2729 = -t2731 * t3024 + t2734 * t3028;
t2728 = -t3026 * t2730 + t3030 * t3093;
t2727 = t3030 * t2730 + t3026 * t3093;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t3014, -t3013, 0, -t3016 * t3034 + t3017 * t3038, 0, 0, 0, 0, 0, 0, -t2965 * t3034 + t2977 * t3038, -t2966 * t3034 + t2978 * t3038, -t2960 * t3034 + t2974 * t3038, -t2925 * t3034 + t2934 * t3038, 0, 0, 0, 0, 0, 0, -t2860 * t3034 + t2878 * t3038, -t2862 * t3034 + t2885 * t3038, -t2845 * t3034 + t2868 * t3038, -t2798 * t3034 + t2820 * t3038, 0, 0, 0, 0, 0, 0, -t2805 * t3034 + t2826 * t3038, -t2810 * t3034 + t2829 * t3038, -t2803 * t3034 + t2825 * t3038, -t2752 * t3034 + t2768 * t3038, 0, 0, 0, 0, 0, 0, -t2767 * t3034 + t2781 * t3038, -t2770 * t3034 + t2786 * t3038, -t2757 * t3034 + t2771 * t3038, -t2733 * t3034 + t2735 * t3038, 0, 0, 0, 0, 0, 0, -t2741 * t3034 + t2747 * t3038, -t2745 * t3034 + t2748 * t3038, -t2737 * t3034 + t2746 * t3038, -t2728 * t3034 + t2729 * t3038; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t3013, t3014, 0, t3016 * t3038 + t3017 * t3034, 0, 0, 0, 0, 0, 0, t2965 * t3038 + t2977 * t3034, t2966 * t3038 + t2978 * t3034, t2960 * t3038 + t2974 * t3034, t2925 * t3038 + t2934 * t3034, 0, 0, 0, 0, 0, 0, t2860 * t3038 + t2878 * t3034, t2862 * t3038 + t2885 * t3034, t2845 * t3038 + t2868 * t3034, t2798 * t3038 + t2820 * t3034, 0, 0, 0, 0, 0, 0, t2805 * t3038 + t2826 * t3034, t2810 * t3038 + t2829 * t3034, t2803 * t3038 + t2825 * t3034, t2752 * t3038 + t2768 * t3034, 0, 0, 0, 0, 0, 0, t2767 * t3038 + t2781 * t3034, t2770 * t3038 + t2786 * t3034, t2757 * t3038 + t2771 * t3034, t2733 * t3038 + t2735 * t3034, 0, 0, 0, 0, 0, 0, t2741 * t3038 + t2747 * t3034, t2745 * t3038 + t2748 * t3034, t2737 * t3038 + t2746 * t3034, t2728 * t3038 + t2729 * t3034; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2963, t2964, t2959, t2924, 0, 0, 0, 0, 0, 0, t2859, t2861, t2844, t2797, 0, 0, 0, 0, 0, 0, t2804, t2809, t2802, t2751, 0, 0, 0, 0, 0, 0, t2766, t2769, t2756, t2732, 0, 0, 0, 0, 0, 0, t2740, t2744, t2736, t2727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3039, -qJDD(1), 0, t3017, 0, 0, 0, 0, 0, 0, t2977, t2978, t2974, t2934, 0, 0, 0, 0, 0, 0, t2878, t2885, t2868, t2820, 0, 0, 0, 0, 0, 0, t2826, t2829, t2825, t2768, 0, 0, 0, 0, 0, 0, t2781, t2786, t2771, t2735, 0, 0, 0, 0, 0, 0, t2747, t2748, t2746, t2729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t3039, 0, t3016, 0, 0, 0, 0, 0, 0, t2965, t2966, t2960, t2925, 0, 0, 0, 0, 0, 0, t2860, t2862, t2845, t2798, 0, 0, 0, 0, 0, 0, t2805, t2810, t2803, t2752, 0, 0, 0, 0, 0, 0, t2767, t2770, t2757, t2733, 0, 0, 0, 0, 0, 0, t2741, t2745, t2737, t2728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2963, t2964, t2959, t2924, 0, 0, 0, 0, 0, 0, t2859, t2861, t2844, t2797, 0, 0, 0, 0, 0, 0, t2804, t2809, t2802, t2751, 0, 0, 0, 0, 0, 0, t2766, t2769, t2756, t2732, 0, 0, 0, 0, 0, 0, t2740, t2744, t2736, t2727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3001, t2998, t2995, t2971, 0, 0, 0, 0, 0, 0, t2919, t2923, t2909, t2863, 0, 0, 0, 0, 0, 0, t2870, t2874, t2869, t2799, 0, 0, 0, 0, 0, 0, t2823, t2824, t2806, t2753, 0, 0, 0, 0, 0, 0, t2775, t2778, t2772, t2734; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2997, t3000, t2994, t2970, 0, 0, 0, 0, 0, 0, t2903, t2905, t2894, t2843, 0, 0, 0, 0, 0, 0, t2851, t2857, t2849, t2780, 0, 0, 0, 0, 0, 0, t2794, t2801, t2785, t2743, 0, 0, 0, 0, 0, 0, t2763, t2765, t2759, t2731; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2996, t2993, t2999, t2986, 0, 0, 0, 0, 0, 0, t2902, t2904, t2893, t2842, 0, 0, 0, 0, 0, 0, t2850, t2856, t2848, t2779, 0, 0, 0, 0, 0, 0, t2793, t2800, t2784, t2742, 0, 0, 0, 0, 0, 0, t2762, t2764, t2758, t2730; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2952, t2953, t2939, t2892, 0, 0, 0, 0, 0, 0, t2896, t2899, t2889, t2812, 0, 0, 0, 0, 0, 0, t2833, t2841, t2816, t2755, 0, 0, 0, 0, 0, 0, t2788, t2790, t2783, t2739; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2954, t2955, t2941, -t3094, 0, 0, 0, 0, 0, 0, -t2912, -t2914, -t2918, -t2875, 0, 0, 0, 0, 0, 0, -t2879, t3100, -t2890, -t2858, 0, 0, 0, 0, 0, 0, -t2834, -t2836, -t2821, -t2760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2938, -t3097, t2945, t2910, 0, 0, 0, 0, 0, 0, t2895, t2898, t2888, t2811, 0, 0, 0, 0, 0, 0, t2832, t2840, t2815, t2754, 0, 0, 0, 0, 0, 0, t2787, t2789, t2782, t2738; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2927, t2917, t2913, t2839, 0, 0, 0, 0, 0, 0, t2867, t2872, t2847, t2774, 0, 0, 0, 0, 0, 0, t2814, t2818, t2808, t2750; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2916, t2933, t2915, t2838, 0, 0, 0, 0, 0, 0, t2866, t2871, t2846, t2773, 0, 0, 0, 0, 0, 0, t2813, t2817, t2807, t2749; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2912, t2914, t2918, t2875, 0, 0, 0, 0, 0, 0, t2879, -t3100, t2890, t2858, 0, 0, 0, 0, 0, 0, t2834, t2836, t2821, t2760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2907, t2901, t2880, t2796, 0, 0, 0, 0, 0, 0, t2835, t2837, t2822, t2761; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2900, t2911, t2882, t2795, 0, 0, 0, 0, 0, 0, -t2852, -t2854, -t2877, -t2791; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2879, -t3100, t2890, t2858, 0, 0, 0, 0, 0, 0, t2834, t2836, t2821, t2760; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2886, t2865, t2853, t2777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2864, t2887, t2855, t2776; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2852, t2854, t2877, t2791;];
f_new_reg  = t1;