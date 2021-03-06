% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 15:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRPRR13_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR13_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_invdynf_fixb_reg2_snew_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:58:38
% EndTime: 2019-05-07 15:59:07
% DurationCPUTime: 29.88s
% Computational Cost: add. (405571->418), mult. (1016524->672), div. (0->0), fcn. (867810->16), ass. (0->361)
t3094 = sin(pkin(7));
t3097 = cos(pkin(7));
t3095 = sin(pkin(6));
t3102 = sin(qJ(2));
t3107 = cos(qJ(2));
t3098 = cos(pkin(6));
t3103 = sin(qJ(1));
t3108 = cos(qJ(1));
t3082 = t3103 * g(1) - t3108 * g(2);
t3109 = qJD(1) ^ 2;
t3117 = pkin(9) * t3095 * t3109 + qJDD(1) * pkin(1) + t3082;
t3115 = t3098 * t3117;
t3083 = -g(1) * t3108 - g(2) * t3103;
t3171 = qJDD(1) * t3095;
t3116 = -pkin(1) * t3109 + pkin(9) * t3171 + t3083;
t3038 = -t3102 * t3116 + (-t3095 * g(3) + t3115) * t3107;
t3190 = qJD(1) * t3095;
t3168 = t3107 * t3190;
t3069 = qJD(2) * t3168 + t3102 * t3171;
t3087 = t3098 * qJDD(1) + qJDD(2);
t3191 = pkin(10) * t3094;
t3119 = (-pkin(2) * t3107 - t3102 * t3191) * t3190;
t3088 = qJD(1) * t3098 + qJD(2);
t3176 = t3097 * t3088;
t3170 = pkin(10) * t3176;
t3192 = t3088 ^ 2;
t3111 = t3087 * pkin(2) + (-t3069 * t3097 + t3094 * t3192) * pkin(10) + (-t3102 * t3119 + t3107 * t3170) * t3190 + t3038;
t3057 = -t3098 * g(3) - t3095 * t3117;
t3177 = t3095 * t3102;
t3169 = qJD(1) * t3177;
t3070 = -qJD(2) * t3169 + t3107 * t3171;
t3162 = t3097 * t3168;
t3113 = -t3069 * t3191 - t3070 * pkin(2) + (-pkin(10) * t3097 * t3169 + pkin(2) * t3088) * t3169 - (t3088 * t3094 + t3162) * pkin(10) * t3168 + t3057;
t3204 = t3094 * t3113 + t3097 * t3111;
t3106 = cos(qJ(3));
t3101 = sin(qJ(3));
t3175 = t3097 * t3101;
t3179 = t3094 * t3101;
t3051 = t3088 * t3179 + (t3102 * t3106 + t3107 * t3175) * t3190;
t3062 = t3094 * t3168 - qJD(3) - t3176;
t3093 = sin(pkin(13));
t3096 = cos(pkin(13));
t3036 = -t3093 * t3051 - t3062 * t3096;
t3037 = t3051 * t3096 - t3062 * t3093;
t3100 = sin(qJ(5));
t3105 = cos(qJ(5));
t3002 = -t3105 * t3036 + t3037 * t3100;
t3000 = qJD(6) + t3002;
t3203 = qJD(6) + t3000;
t3004 = t3036 * t3100 + t3037 * t3105;
t3178 = t3094 * t3106;
t3049 = -t3088 * t3178 + t3101 * t3169 - t3106 * t3162;
t3047 = qJD(5) + t3049;
t3099 = sin(qJ(6));
t3104 = cos(qJ(6));
t2986 = t3004 * t3099 - t3104 * t3047;
t3202 = t2986 ^ 2;
t2988 = t3004 * t3104 + t3047 * t3099;
t3201 = t2988 ^ 2;
t3200 = t3000 ^ 2;
t3199 = t3002 ^ 2;
t3198 = t3004 ^ 2;
t3197 = t3036 ^ 2;
t3196 = t3037 ^ 2;
t3195 = t3047 ^ 2;
t3031 = t3049 ^ 2;
t3194 = t3051 ^ 2;
t3193 = t3062 ^ 2;
t3189 = t2986 * t2988;
t3188 = t3002 * t3004;
t3187 = t3036 * t3037;
t3186 = t3037 * t3049;
t3185 = t3049 * t3036;
t3184 = t3049 * t3051;
t3183 = t3049 * t3062;
t3182 = t3051 * t3062;
t3181 = t3070 * t3097;
t3180 = t3095 ^ 2 * t3109;
t3174 = qJD(5) - t3047;
t3173 = qJD(6) - t3000;
t3172 = t3102 * t3115 + t3107 * t3116;
t3005 = -t3192 * pkin(2) + (t3087 * t3094 + t3181) * pkin(10) + (-t3102 * g(3) + (t3102 * t3170 + t3107 * t3119) * qJD(1)) * t3095 + t3172;
t2960 = t3106 * t3005 + t3101 * t3204;
t3030 = pkin(3) * t3049 - qJ(4) * t3051;
t3121 = -t3070 * t3094 + t3087 * t3097 + qJDD(3);
t2946 = -pkin(3) * t3193 + qJ(4) * t3121 - t3049 * t3030 + t2960;
t2975 = -t3094 * t3111 + t3097 * t3113;
t3161 = t3051 * qJD(3) + t3101 * t3069 - t3087 * t3178 - t3106 * t3181;
t2991 = t3161 - t3182;
t3022 = -t3049 * qJD(3) + t3106 * t3069 + t3070 * t3175 + t3087 * t3179;
t3164 = -t3022 - t3183;
t2948 = pkin(3) * t2991 + qJ(4) * t3164 + t2975;
t2904 = 0.2e1 * qJD(4) * t3036 + t3096 * t2946 + t3093 * t2948;
t3008 = -t3022 * t3093 + t3096 * t3121;
t3017 = pkin(4) * t3049 - pkin(11) * t3037;
t2896 = -pkin(4) * t3197 + pkin(11) * t3008 - t3017 * t3049 + t2904;
t2903 = -0.2e1 * qJD(4) * t3037 - t3093 * t2946 + t3096 * t2948;
t3009 = t3096 * t3022 + t3093 * t3121;
t2979 = -t3009 + t3185;
t2980 = t3161 + t3187;
t3114 = pkin(4) * t2980 + pkin(11) * t2979 + t2903;
t2861 = t3105 * t2896 + t3100 * t3114;
t2860 = -t2896 * t3100 + t3105 * t3114;
t3129 = -t3100 * t3008 - t3105 * t3009;
t2958 = -qJD(5) * t3002 - t3129;
t3167 = t3047 * t3002 - t2958;
t3122 = qJDD(5) + t3161;
t3166 = -t3099 * t2958 + t3104 * t3122;
t3165 = -t3105 * t3008 + t3100 * t3009;
t3163 = t3088 * t3168;
t3160 = t3101 * t3005 - t3106 * t3204;
t2974 = pkin(5) * t3002 - pkin(12) * t3004;
t2857 = -pkin(5) * t3195 + pkin(12) * t3122 - t3002 * t2974 + t2861;
t2944 = -t3121 * pkin(3) - t3193 * qJ(4) + t3051 * t3030 + qJDD(4) + t3160;
t2925 = -t3008 * pkin(4) - t3197 * pkin(11) + t3037 * t3017 + t2944;
t2939 = (qJD(5) + t3047) * t3004 + t3165;
t2884 = pkin(5) * t2939 + pkin(12) * t3167 + t2925;
t2840 = -t2857 * t3099 + t2884 * t3104;
t2841 = t2857 * t3104 + t2884 * t3099;
t2826 = -t2840 * t3099 + t2841 * t3104;
t2856 = -pkin(5) * t3122 - pkin(12) * t3195 + t2974 * t3004 - t2860;
t2814 = t2826 * t3100 - t2856 * t3105;
t2815 = t2826 * t3105 + t2856 * t3100;
t2803 = t2814 * t3096 + t2815 * t3093;
t2804 = -t2814 * t3093 + t2815 * t3096;
t2825 = t2840 * t3104 + t2841 * t3099;
t3158 = t2804 * t3101 - t2825 * t3106;
t2796 = -t3094 * t2803 + t3097 * t3158;
t2799 = t2804 * t3106 + t2825 * t3101;
t3159 = t2796 * t3107 + t2799 * t3102;
t2838 = t2860 * t3105 + t2861 * t3100;
t2839 = -t2860 * t3100 + t2861 * t3105;
t2821 = t2838 * t3096 + t2839 * t3093;
t2822 = -t2838 * t3093 + t2839 * t3096;
t3156 = t2822 * t3101 - t2925 * t3106;
t2806 = -t3094 * t2821 + t3097 * t3156;
t2818 = t2822 * t3106 + t2925 * t3101;
t3157 = t2806 * t3107 + t2818 * t3102;
t2918 = -t2988 * t3173 + t3166;
t3118 = -t3104 * t2958 - t3099 * t3122;
t2920 = t2986 * t3173 + t3118;
t2886 = t2918 * t3104 - t2920 * t3099;
t2950 = -t3201 - t3202;
t2874 = t2886 * t3100 - t2950 * t3105;
t2875 = t2886 * t3105 + t2950 * t3100;
t2847 = t2874 * t3096 + t2875 * t3093;
t2848 = -t2874 * t3093 + t2875 * t3096;
t2885 = t2918 * t3099 + t2920 * t3104;
t3151 = t2848 * t3101 - t2885 * t3106;
t2824 = -t3094 * t2847 + t3097 * t3151;
t2837 = t2848 * t3106 + t2885 * t3101;
t3155 = t2824 * t3107 + t2837 * t3102;
t3120 = -qJD(5) * t3004 - qJDD(6) - t3165;
t2929 = -t3120 - t3189;
t2951 = -t3200 - t3202;
t2900 = -t2929 * t3099 + t2951 * t3104;
t2917 = t2988 * t3203 - t3166;
t2880 = t2900 * t3100 - t2917 * t3105;
t2881 = t2900 * t3105 + t2917 * t3100;
t2852 = t2880 * t3096 + t2881 * t3093;
t2853 = -t2880 * t3093 + t2881 * t3096;
t2899 = t2929 * t3104 + t2951 * t3099;
t3149 = t2853 * t3101 - t2899 * t3106;
t2828 = -t3094 * t2852 + t3097 * t3149;
t2842 = t2853 * t3106 + t2899 * t3101;
t3154 = t2828 * t3107 + t2842 * t3102;
t2930 = t3120 - t3189;
t2954 = -t3200 - t3201;
t2902 = t2930 * t3104 - t2954 * t3099;
t2919 = -t2986 * t3203 - t3118;
t2882 = t2902 * t3100 - t2919 * t3105;
t2883 = t2902 * t3105 + t2919 * t3100;
t2854 = t2882 * t3096 + t2883 * t3093;
t2855 = -t2882 * t3093 + t2883 * t3096;
t2901 = t2930 * t3099 + t2954 * t3104;
t3148 = t2855 * t3101 - t2901 * t3106;
t2830 = -t3094 * t2854 + t3097 * t3148;
t2843 = t2855 * t3106 + t2901 * t3101;
t3153 = t2830 * t3107 + t2843 * t3102;
t2876 = t2903 * t3096 + t2904 * t3093;
t2877 = -t2903 * t3093 + t2904 * t3096;
t3145 = t2877 * t3101 - t2944 * t3106;
t2845 = -t3094 * t2876 + t3097 * t3145;
t2862 = t2877 * t3106 + t2944 * t3101;
t3152 = t2845 * t3107 + t2862 * t3102;
t2940 = -t3004 * t3174 - t3165;
t2942 = t3002 * t3174 + t3129;
t2907 = t2940 * t3100 + t2942 * t3105;
t2908 = t2940 * t3105 - t2942 * t3100;
t2878 = t2907 * t3096 + t2908 * t3093;
t2879 = -t2907 * t3093 + t2908 * t3096;
t2957 = -t3198 - t3199;
t3144 = t2879 * t3101 - t2957 * t3106;
t2850 = -t3094 * t2878 + t3097 * t3144;
t2869 = t2879 * t3106 + t2957 * t3101;
t3150 = t2850 * t3107 + t2869 * t3102;
t2963 = t3122 - t3188;
t2972 = -t3195 - t3199;
t2932 = t2963 * t3105 + t2972 * t3100;
t2933 = -t2963 * t3100 + t2972 * t3105;
t2897 = t2932 * t3096 + t2933 * t3093;
t2898 = -t2932 * t3093 + t2933 * t3096;
t3143 = t2898 * t3101 - t2939 * t3106;
t2859 = -t3094 * t2897 + t3097 * t3143;
t2887 = t2898 * t3106 + t2939 * t3101;
t3147 = t2859 * t3107 + t2887 * t3102;
t2964 = -t3122 - t3188;
t2982 = -t3195 - t3198;
t2937 = t2964 * t3100 + t2982 * t3105;
t2938 = t2964 * t3105 - t2982 * t3100;
t2905 = t2937 * t3096 + t2938 * t3093;
t2906 = -t2937 * t3093 + t2938 * t3096;
t3142 = t2906 * t3101 + t3106 * t3167;
t2868 = -t3094 * t2905 + t3097 * t3142;
t2890 = t2906 * t3106 - t3101 * t3167;
t3146 = t2868 * t3107 + t2890 * t3102;
t2977 = t3008 + t3186;
t2952 = t2977 * t3093 + t2979 * t3096;
t2953 = t2977 * t3096 - t2979 * t3093;
t2984 = -t3196 - t3197;
t3137 = t2953 * t3101 - t2984 * t3106;
t2912 = -t3094 * t2952 + t3097 * t3137;
t2934 = t2953 * t3106 + t2984 * t3101;
t3141 = t2912 * t3107 + t2934 * t3102;
t3135 = t2960 * t3101 - t3106 * t3160;
t2914 = -t3094 * t2975 + t3097 * t3135;
t2928 = t2960 * t3106 + t3101 * t3160;
t3140 = t2914 * t3107 + t2928 * t3102;
t2990 = -t3031 - t3197;
t2961 = t2980 * t3096 + t2990 * t3093;
t2962 = -t2980 * t3093 + t2990 * t3096;
t2976 = -t3008 + t3186;
t3134 = t2962 * t3101 - t2976 * t3106;
t2916 = -t3094 * t2961 + t3097 * t3134;
t2935 = t2962 * t3106 + t2976 * t3101;
t3139 = t2916 * t3107 + t2935 * t3102;
t2981 = -t3161 + t3187;
t3001 = -t3031 - t3196;
t2965 = t2981 * t3093 + t3001 * t3096;
t2966 = t2981 * t3096 - t3001 * t3093;
t2978 = t3009 + t3185;
t3133 = t2966 * t3101 - t2978 * t3106;
t2922 = -t3094 * t2965 + t3097 * t3133;
t2936 = t2966 * t3106 + t2978 * t3101;
t3138 = t2922 * t3107 + t2936 * t3102;
t3013 = -t3194 - t3031;
t2992 = -t3161 - t3182;
t2994 = -t3022 + t3183;
t3130 = t2992 * t3101 + t2994 * t3106;
t2956 = -t3094 * t3013 + t3097 * t3130;
t2973 = t2992 * t3106 - t2994 * t3101;
t3136 = t2956 * t3107 + t2973 * t3102;
t3016 = t3121 - t3184;
t3021 = -t3031 - t3193;
t3127 = t3016 * t3106 + t3021 * t3101;
t2968 = -t3094 * t2991 + t3097 * t3127;
t2983 = -t3016 * t3101 + t3021 * t3106;
t3132 = t2968 * t3107 + t2983 * t3102;
t3015 = -t3121 - t3184;
t3027 = -t3193 - t3194;
t3128 = t3015 * t3101 + t3027 * t3106;
t2970 = t3094 * t3164 + t3097 * t3128;
t2985 = t3015 * t3106 - t3027 * t3101;
t3131 = t2970 * t3107 + t2985 * t3102;
t3039 = -g(3) * t3177 + t3172;
t3126 = t3038 * t3107 + t3039 * t3102;
t3054 = t3163 - t3069;
t3073 = t3088 * t3169;
t3055 = t3070 + t3073;
t3125 = t3054 * t3107 + t3055 * t3102;
t3091 = t3102 ^ 2;
t3063 = -t3091 * t3180 - t3192;
t3081 = t3107 * t3102 * t3180;
t3068 = t3081 - t3087;
t3124 = t3063 * t3107 + t3068 * t3102;
t3067 = t3081 + t3087;
t3092 = t3107 ^ 2;
t3071 = -t3092 * t3180 - t3192;
t3123 = t3067 * t3107 + t3071 * t3102;
t3080 = -qJDD(1) * t3103 - t3108 * t3109;
t3079 = qJDD(1) * t3108 - t3103 * t3109;
t3072 = (-t3091 - t3092) * t3180;
t3056 = -t3070 + t3073;
t3053 = t3163 + t3069;
t3042 = -t3067 * t3102 + t3071 * t3107;
t3040 = -t3063 * t3102 + t3068 * t3107;
t3033 = -t3054 * t3102 + t3055 * t3107;
t3029 = -t3095 * t3056 + t3098 * t3123;
t3028 = t3098 * t3056 + t3095 * t3123;
t3026 = -t3095 * t3053 + t3098 * t3124;
t3025 = t3098 * t3053 + t3095 * t3124;
t3024 = -t3095 * t3072 + t3098 * t3125;
t3023 = t3098 * t3072 + t3095 * t3125;
t3014 = -t3038 * t3102 + t3039 * t3107;
t2996 = -t3095 * t3057 + t3098 * t3126;
t2995 = t3098 * t3057 + t3095 * t3126;
t2969 = t3094 * t3128 - t3097 * t3164;
t2967 = t3097 * t2991 + t3094 * t3127;
t2955 = t3097 * t3013 + t3094 * t3130;
t2949 = -t2970 * t3102 + t2985 * t3107;
t2945 = -t2968 * t3102 + t2983 * t3107;
t2931 = -t2956 * t3102 + t2973 * t3107;
t2927 = -t3095 * t2969 + t3098 * t3131;
t2926 = t3098 * t2969 + t3095 * t3131;
t2924 = -t3095 * t2967 + t3098 * t3132;
t2923 = t3098 * t2967 + t3095 * t3132;
t2921 = t3097 * t2965 + t3094 * t3133;
t2915 = t3097 * t2961 + t3094 * t3134;
t2913 = t3097 * t2975 + t3094 * t3135;
t2911 = t3097 * t2952 + t3094 * t3137;
t2910 = -t3095 * t2955 + t3098 * t3136;
t2909 = t3098 * t2955 + t3095 * t3136;
t2894 = -t2922 * t3102 + t2936 * t3107;
t2891 = -t2916 * t3102 + t2935 * t3107;
t2889 = -t2912 * t3102 + t2934 * t3107;
t2888 = -t2914 * t3102 + t2928 * t3107;
t2873 = -t3095 * t2921 + t3098 * t3138;
t2872 = t3098 * t2921 + t3095 * t3138;
t2871 = -t3095 * t2915 + t3098 * t3139;
t2870 = t3098 * t2915 + t3095 * t3139;
t2867 = t3097 * t2905 + t3094 * t3142;
t2866 = -t3095 * t2911 + t3098 * t3141;
t2865 = t3098 * t2911 + t3095 * t3141;
t2864 = -t3095 * t2913 + t3098 * t3140;
t2863 = t3098 * t2913 + t3095 * t3140;
t2858 = t3097 * t2897 + t3094 * t3143;
t2851 = -t2868 * t3102 + t2890 * t3107;
t2849 = t3097 * t2878 + t3094 * t3144;
t2846 = -t2859 * t3102 + t2887 * t3107;
t2844 = t3097 * t2876 + t3094 * t3145;
t2836 = -t3095 * t2867 + t3098 * t3146;
t2835 = t3098 * t2867 + t3095 * t3146;
t2834 = -t2850 * t3102 + t2869 * t3107;
t2833 = -t2845 * t3102 + t2862 * t3107;
t2832 = -t3095 * t2858 + t3098 * t3147;
t2831 = t3098 * t2858 + t3095 * t3147;
t2829 = t3097 * t2854 + t3094 * t3148;
t2827 = t3097 * t2852 + t3094 * t3149;
t2823 = t3097 * t2847 + t3094 * t3151;
t2820 = -t3095 * t2849 + t3098 * t3150;
t2819 = t3098 * t2849 + t3095 * t3150;
t2817 = -t3095 * t2844 + t3098 * t3152;
t2816 = t3098 * t2844 + t3095 * t3152;
t2813 = -t2830 * t3102 + t2843 * t3107;
t2812 = -t2828 * t3102 + t2842 * t3107;
t2811 = -t2824 * t3102 + t2837 * t3107;
t2810 = -t3095 * t2829 + t3098 * t3153;
t2809 = t3098 * t2829 + t3095 * t3153;
t2808 = -t3095 * t2827 + t3098 * t3154;
t2807 = t3098 * t2827 + t3095 * t3154;
t2805 = t3097 * t2821 + t3094 * t3156;
t2802 = -t3095 * t2823 + t3098 * t3155;
t2801 = t3098 * t2823 + t3095 * t3155;
t2800 = -t2806 * t3102 + t2818 * t3107;
t2798 = -t3095 * t2805 + t3098 * t3157;
t2797 = t3098 * t2805 + t3095 * t3157;
t2795 = t3097 * t2803 + t3094 * t3158;
t2794 = -t2796 * t3102 + t2799 * t3107;
t2793 = -t3095 * t2795 + t3098 * t3159;
t2792 = t3098 * t2795 + t3095 * t3159;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t3080, -t3079, 0, -t3082 * t3103 + t3083 * t3108, 0, 0, 0, 0, 0, 0, -t3029 * t3103 + t3042 * t3108, -t3026 * t3103 + t3040 * t3108, -t3024 * t3103 + t3033 * t3108, -t2996 * t3103 + t3014 * t3108, 0, 0, 0, 0, 0, 0, -t2924 * t3103 + t2945 * t3108, -t2927 * t3103 + t2949 * t3108, -t2910 * t3103 + t2931 * t3108, -t2864 * t3103 + t2888 * t3108, 0, 0, 0, 0, 0, 0, -t2871 * t3103 + t2891 * t3108, -t2873 * t3103 + t2894 * t3108, -t2866 * t3103 + t2889 * t3108, -t2817 * t3103 + t2833 * t3108, 0, 0, 0, 0, 0, 0, -t2832 * t3103 + t2846 * t3108, -t2836 * t3103 + t2851 * t3108, -t2820 * t3103 + t2834 * t3108, -t2798 * t3103 + t2800 * t3108, 0, 0, 0, 0, 0, 0, -t2808 * t3103 + t2812 * t3108, -t2810 * t3103 + t2813 * t3108, -t2802 * t3103 + t2811 * t3108, -t2793 * t3103 + t2794 * t3108; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t3079, t3080, 0, t3082 * t3108 + t3083 * t3103, 0, 0, 0, 0, 0, 0, t3029 * t3108 + t3042 * t3103, t3026 * t3108 + t3040 * t3103, t3024 * t3108 + t3033 * t3103, t2996 * t3108 + t3014 * t3103, 0, 0, 0, 0, 0, 0, t2924 * t3108 + t2945 * t3103, t2927 * t3108 + t2949 * t3103, t2910 * t3108 + t2931 * t3103, t2864 * t3108 + t2888 * t3103, 0, 0, 0, 0, 0, 0, t2871 * t3108 + t2891 * t3103, t2873 * t3108 + t2894 * t3103, t2866 * t3108 + t2889 * t3103, t2817 * t3108 + t2833 * t3103, 0, 0, 0, 0, 0, 0, t2832 * t3108 + t2846 * t3103, t2836 * t3108 + t2851 * t3103, t2820 * t3108 + t2834 * t3103, t2798 * t3108 + t2800 * t3103, 0, 0, 0, 0, 0, 0, t2808 * t3108 + t2812 * t3103, t2810 * t3108 + t2813 * t3103, t2802 * t3108 + t2811 * t3103, t2793 * t3108 + t2794 * t3103; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t3028, t3025, t3023, t2995, 0, 0, 0, 0, 0, 0, t2923, t2926, t2909, t2863, 0, 0, 0, 0, 0, 0, t2870, t2872, t2865, t2816, 0, 0, 0, 0, 0, 0, t2831, t2835, t2819, t2797, 0, 0, 0, 0, 0, 0, t2807, t2809, t2801, t2792; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3109, -qJDD(1), 0, t3083, 0, 0, 0, 0, 0, 0, t3042, t3040, t3033, t3014, 0, 0, 0, 0, 0, 0, t2945, t2949, t2931, t2888, 0, 0, 0, 0, 0, 0, t2891, t2894, t2889, t2833, 0, 0, 0, 0, 0, 0, t2846, t2851, t2834, t2800, 0, 0, 0, 0, 0, 0, t2812, t2813, t2811, t2794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t3109, 0, t3082, 0, 0, 0, 0, 0, 0, t3029, t3026, t3024, t2996, 0, 0, 0, 0, 0, 0, t2924, t2927, t2910, t2864, 0, 0, 0, 0, 0, 0, t2871, t2873, t2866, t2817, 0, 0, 0, 0, 0, 0, t2832, t2836, t2820, t2798, 0, 0, 0, 0, 0, 0, t2808, t2810, t2802, t2793; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t3028, t3025, t3023, t2995, 0, 0, 0, 0, 0, 0, t2923, t2926, t2909, t2863, 0, 0, 0, 0, 0, 0, t2870, t2872, t2865, t2816, 0, 0, 0, 0, 0, 0, t2831, t2835, t2819, t2797, 0, 0, 0, 0, 0, 0, t2807, t2809, t2801, t2792; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3071, t3068, t3055, t3039, 0, 0, 0, 0, 0, 0, t2983, t2985, t2973, t2928, 0, 0, 0, 0, 0, 0, t2935, t2936, t2934, t2862, 0, 0, 0, 0, 0, 0, t2887, t2890, t2869, t2818, 0, 0, 0, 0, 0, 0, t2842, t2843, t2837, t2799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3067, t3063, t3054, t3038, 0, 0, 0, 0, 0, 0, t2968, t2970, t2956, t2914, 0, 0, 0, 0, 0, 0, t2916, t2922, t2912, t2845, 0, 0, 0, 0, 0, 0, t2859, t2868, t2850, t2806, 0, 0, 0, 0, 0, 0, t2828, t2830, t2824, t2796; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3056, t3053, t3072, t3057, 0, 0, 0, 0, 0, 0, t2967, t2969, t2955, t2913, 0, 0, 0, 0, 0, 0, t2915, t2921, t2911, t2844, 0, 0, 0, 0, 0, 0, t2858, t2867, t2849, t2805, 0, 0, 0, 0, 0, 0, t2827, t2829, t2823, t2795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3021, t3015, t2992, t2960, 0, 0, 0, 0, 0, 0, t2962, t2966, t2953, t2877, 0, 0, 0, 0, 0, 0, t2898, t2906, t2879, t2822, 0, 0, 0, 0, 0, 0, t2853, t2855, t2848, t2804; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3016, t3027, t2994, -t3160, 0, 0, 0, 0, 0, 0, -t2976, -t2978, -t2984, -t2944, 0, 0, 0, 0, 0, 0, -t2939, t3167, -t2957, -t2925, 0, 0, 0, 0, 0, 0, -t2899, -t2901, -t2885, -t2825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2991, -t3164, t3013, t2975, 0, 0, 0, 0, 0, 0, t2961, t2965, t2952, t2876, 0, 0, 0, 0, 0, 0, t2897, t2905, t2878, t2821, 0, 0, 0, 0, 0, 0, t2852, t2854, t2847, t2803; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2990, t2981, t2977, t2904, 0, 0, 0, 0, 0, 0, t2933, t2938, t2908, t2839, 0, 0, 0, 0, 0, 0, t2881, t2883, t2875, t2815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2980, t3001, t2979, t2903, 0, 0, 0, 0, 0, 0, t2932, t2937, t2907, t2838, 0, 0, 0, 0, 0, 0, t2880, t2882, t2874, t2814; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2976, t2978, t2984, t2944, 0, 0, 0, 0, 0, 0, t2939, -t3167, t2957, t2925, 0, 0, 0, 0, 0, 0, t2899, t2901, t2885, t2825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2972, t2964, t2940, t2861, 0, 0, 0, 0, 0, 0, t2900, t2902, t2886, t2826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2963, t2982, t2942, t2860, 0, 0, 0, 0, 0, 0, -t2917, -t2919, -t2950, -t2856; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2939, -t3167, t2957, t2925, 0, 0, 0, 0, 0, 0, t2899, t2901, t2885, t2825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2951, t2930, t2918, t2841; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2929, t2954, t2920, t2840; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2917, t2919, t2950, t2856;];
f_new_reg  = t1;
