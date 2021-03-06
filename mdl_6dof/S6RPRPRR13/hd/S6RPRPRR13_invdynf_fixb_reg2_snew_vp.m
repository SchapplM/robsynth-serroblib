% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRPRR13
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRPRR13_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynf_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:02:38
% EndTime: 2019-05-05 21:02:58
% DurationCPUTime: 22.15s
% Computational Cost: add. (104819->396), mult. (334086->585), div. (0->0), fcn. (284014->14), ass. (0->332)
t3024 = sin(qJ(3));
t3016 = sin(pkin(12));
t3018 = sin(pkin(6));
t3114 = t3016 * t3018;
t3100 = qJD(1) * t3114;
t3021 = cos(pkin(6));
t3028 = cos(qJ(3));
t3019 = cos(pkin(12));
t3020 = cos(pkin(7));
t3110 = t3018 * t3020;
t3096 = t3019 * t3110;
t3017 = sin(pkin(7));
t3112 = t3017 * t3028;
t3140 = t3021 * t3112 + t3028 * t3096;
t2967 = -qJD(1) * t3140 + t3024 * t3100;
t3111 = t3018 * t3019;
t3099 = qJD(1) * t3111;
t2982 = -t3020 * t3021 * qJD(1) + t3017 * t3099 - qJD(3);
t3107 = qJD(3) - t2982;
t3108 = t3020 * t3024;
t3113 = t3017 * t3024;
t3036 = t3018 * (t3016 * t3028 + t3019 * t3108) + t3021 * t3113;
t3143 = t3036 * qJDD(1);
t2917 = -t2967 * t3107 + t3143;
t3102 = qJDD(1) * t3018;
t3091 = t3019 * t3102;
t3101 = qJDD(1) * t3021;
t3042 = t3017 * t3091 - t3020 * t3101 - qJDD(3);
t2969 = qJD(1) * t3036;
t3119 = t2967 * t2969;
t3039 = t3042 - t3119;
t2944 = t2969 ^ 2;
t3128 = t2982 ^ 2;
t3093 = -t3128 - t2944;
t3053 = t3024 * t3039 + t3028 * t3093;
t2879 = t3020 * t2917 + t3017 * t3053;
t2882 = t3017 * t2917 - t3020 * t3053;
t2898 = t3024 * t3093 - t3028 * t3039;
t3059 = t2882 * t3019 + t2898 * t3016;
t2833 = t3018 * t2879 + t3021 * t3059;
t2860 = t2882 * t3016 - t2898 * t3019;
t3025 = sin(qJ(1));
t3029 = cos(qJ(1));
t3160 = t2833 * t3025 + t2860 * t3029;
t3159 = t2833 * t3029 - t2860 * t3025;
t2831 = -t3021 * t2879 + t3018 * t3059;
t3106 = qJD(3) + t2982;
t2916 = -t2967 * t3106 + t3143;
t3154 = t2916 * t3024;
t3153 = t2916 * t3028;
t2943 = -t2967 * qJD(3) + t3143;
t3152 = qJDD(5) + t2943;
t2966 = t2967 ^ 2;
t2929 = -t3128 - t2966;
t2932 = t3042 + t3119;
t2894 = t2929 * t3028 + t2932 * t3024;
t3151 = t2894 * t3016;
t3150 = t2894 * t3019;
t3056 = -t2929 * t3024 + t2932 * t3028;
t3149 = t3017 * t3056;
t3148 = t3020 * t3056;
t3127 = -2 * qJD(4);
t3136 = -t2944 - t2966;
t3142 = t3017 * t3136;
t3030 = qJD(1) ^ 2;
t3109 = t3018 * t3030;
t3141 = t3020 * t3136;
t3023 = sin(qJ(5));
t3027 = cos(qJ(5));
t2952 = -t3027 * t2967 - t2982 * t3023;
t2949 = qJD(6) + t2952;
t3137 = qJD(6) + t2949;
t2954 = t2967 * t3023 - t2982 * t3027;
t2965 = qJD(5) + t2969;
t3022 = sin(qJ(6));
t3026 = cos(qJ(6));
t2925 = t2954 * t3022 - t3026 * t2965;
t3134 = t2925 ^ 2;
t2927 = t2954 * t3026 + t2965 * t3022;
t3133 = t2927 ^ 2;
t3132 = t2949 ^ 2;
t3131 = t2952 ^ 2;
t3130 = t2954 ^ 2;
t3129 = t2965 ^ 2;
t3126 = pkin(9) * t3016;
t3125 = pkin(9) * qJDD(1);
t3124 = t2925 * t2927;
t3121 = t2954 * t2952;
t3118 = t2967 * t2982;
t3117 = t2982 * t2969;
t3007 = t3025 * g(1) - g(2) * t3029;
t2988 = qJDD(1) * pkin(1) + qJ(2) * t3109 + t3007;
t3116 = t2988 * t3021;
t3013 = t3018 ^ 2;
t3115 = t3013 * t3030;
t3105 = qJD(5) - t2965;
t3104 = qJD(5) + t2965;
t3103 = qJD(6) - t2949;
t3092 = t3016 * t3102;
t3043 = qJDD(1) * t3140 - t3024 * t3092;
t2942 = qJD(3) * t2969 - t3043;
t2955 = pkin(4) * t2969 + pkin(10) * t2982;
t2981 = (t3017 * t3021 + t3096) * qJD(1) * pkin(9);
t3008 = -g(1) * t3029 - g(2) * t3025;
t2989 = -pkin(1) * t3030 + qJ(2) * t3102 + t3008;
t3046 = -pkin(2) * t3019 - t3017 * t3126;
t3045 = t3020 * t3125 + t3046 * t3109;
t3047 = -g(3) * t3111 - 0.2e1 * qJD(2) * t3100 + t3019 * t3116;
t2908 = (pkin(2) * qJDD(1) + qJD(1) * t2981) * t3021 + (-t3018 * t3045 - t2989) * t3016 + t3047;
t2990 = (pkin(2) * t3021 - t3110 * t3126) * qJD(1);
t3090 = -t3021 * g(3) + qJDD(2);
t2928 = (-t2988 + t3046 * qJDD(1) + (-t2981 * t3019 + t2990 * t3016) * qJD(1)) * t3018 + t3090;
t2890 = -t3017 * t2908 + t3020 * t2928;
t3038 = t2969 * t3127 + (-t2943 - t3118) * qJ(4) - pkin(3) * t3117 + t2890;
t2839 = -t2966 * pkin(4) - t2969 * t2955 + (pkin(3) + pkin(10)) * t2942 + t3038;
t2939 = pkin(3) * t2967 - qJ(4) * t2969;
t3094 = 0.2e1 * qJD(2) * t3099 + t3019 * t2989 + t3016 * t3116;
t2909 = (-qJD(1) * t2990 + t3017 * t3125) * t3021 + (-g(3) * t3016 + t3019 * t3045) * t3018 + t3094;
t3085 = -t3020 * t3028 * t2908 + t3024 * t2909 - t2928 * t3112;
t2850 = pkin(3) * t3042 - t3128 * qJ(4) + t2969 * t2939 + qJDD(4) + t3085;
t3032 = t2932 * pkin(10) + (t2943 - t3118) * pkin(4) + t2850;
t2803 = t3027 * t2839 + t3023 * t3032;
t3095 = t3021 * t3109;
t2863 = t2908 * t3108 + t3028 * t2909 + t2928 * t3113;
t2802 = -t2839 * t3023 + t3027 * t3032;
t3052 = t3023 * t2942 - t3027 * t3042;
t2902 = -qJD(5) * t2952 + t3052;
t3089 = -t3022 * t2902 + t3026 * t3152;
t3088 = -t3027 * t2942 - t3023 * t3042;
t2911 = pkin(5) * t2952 - pkin(11) * t2954;
t2796 = -t3129 * pkin(5) + pkin(11) * t3152 - t2952 * t2911 + t2803;
t3040 = -t3128 * pkin(3) - qJ(4) * t3042 - t2967 * t2939 + t2863;
t2837 = -t2942 * pkin(4) - t2966 * pkin(10) + (t3127 - t2955) * t2982 + t3040;
t2884 = t2954 * t3104 + t3088;
t2806 = (t2952 * t2965 - t2902) * pkin(11) + t2884 * pkin(5) + t2837;
t2776 = -t2796 * t3022 + t2806 * t3026;
t2777 = t2796 * t3026 + t2806 * t3022;
t2766 = -t2776 * t3022 + t2777 * t3026;
t2795 = -pkin(5) * t3152 - pkin(11) * t3129 + t2911 * t2954 - t2802;
t2757 = t2766 * t3027 + t2795 * t3023;
t2756 = t2766 * t3023 - t2795 * t3027;
t2765 = t2776 * t3026 + t2777 * t3022;
t3083 = -t2756 * t3028 + t2765 * t3024;
t2749 = -t3017 * t2757 + t3020 * t3083;
t2752 = t2756 * t3024 + t2765 * t3028;
t3084 = t2749 * t3019 + t2752 * t3016;
t2781 = -t2802 * t3023 + t2803 * t3027;
t2780 = t2802 * t3027 + t2803 * t3023;
t3079 = -t2780 * t3028 + t2837 * t3024;
t2763 = -t3017 * t2781 + t3020 * t3079;
t2773 = t2780 * t3024 + t2837 * t3028;
t3082 = t2763 * t3019 + t2773 * t3016;
t2857 = -t2927 * t3103 + t3089;
t3031 = -t3026 * t2902 - t3022 * t3152;
t2859 = t2925 * t3103 + t3031;
t2815 = t2857 * t3026 - t2859 * t3022;
t2878 = -t3133 - t3134;
t2805 = t2815 * t3027 + t2878 * t3023;
t2804 = t2815 * t3023 - t2878 * t3027;
t2814 = t2857 * t3022 + t2859 * t3026;
t3076 = -t2804 * t3028 + t2814 * t3024;
t2772 = -t3017 * t2805 + t3020 * t3076;
t2791 = t2804 * t3024 + t2814 * t3028;
t3081 = t2772 * t3019 + t2791 * t3016;
t3041 = -qJD(5) * t2954 - qJDD(6) - t3088;
t2868 = -t3041 - t3124;
t2883 = -t3132 - t3134;
t2842 = -t2868 * t3022 + t2883 * t3026;
t2856 = t2927 * t3137 - t3089;
t2810 = t2842 * t3027 + t2856 * t3023;
t2809 = t2842 * t3023 - t2856 * t3027;
t2841 = t2868 * t3026 + t2883 * t3022;
t3074 = -t2809 * t3028 + t2841 * t3024;
t2779 = -t3017 * t2810 + t3020 * t3074;
t2793 = t2809 * t3024 + t2841 * t3028;
t3080 = t2779 * t3019 + t2793 * t3016;
t2869 = t3041 - t3124;
t2891 = -t3132 - t3133;
t2847 = t2869 * t3026 - t2891 * t3022;
t2858 = -t2925 * t3137 - t3031;
t2812 = t2847 * t3027 + t2858 * t3023;
t2811 = t2847 * t3023 - t2858 * t3027;
t2846 = t2869 * t3022 + t2891 * t3026;
t3073 = -t2811 * t3028 + t2846 * t3024;
t2783 = -t3017 * t2812 + t3020 * t3073;
t2797 = t2811 * t3024 + t2846 * t3028;
t3078 = t2783 * t3019 + t2797 * t3016;
t2851 = t2942 * pkin(3) + t3038;
t2849 = t2982 * t3127 + t3040;
t3069 = t2849 * t3024 - t2850 * t3028;
t2799 = -t3017 * t2851 + t3020 * t3069;
t2813 = t2849 * t3028 + t2850 * t3024;
t3077 = t2799 * t3019 + t2813 * t3016;
t2885 = -t2954 * t3105 - t3088;
t2887 = t2952 * t3105 - t3052;
t2853 = t2885 * t3027 - t2887 * t3023;
t2852 = t2885 * t3023 + t2887 * t3027;
t2896 = -t3130 - t3131;
t3068 = -t2852 * t3028 + t2896 * t3024;
t2808 = -t3017 * t2853 + t3020 * t3068;
t2840 = t2852 * t3024 + t2896 * t3028;
t3075 = t2808 * t3019 + t2840 * t3016;
t3067 = t2863 * t3024 - t3028 * t3085;
t2817 = -t3017 * t2890 + t3020 * t3067;
t2834 = t2863 * t3028 + t3024 * t3085;
t3072 = t2817 * t3019 + t2834 * t3016;
t2892 = t3152 - t3121;
t2903 = -t3129 - t3131;
t2871 = -t2892 * t3023 + t2903 * t3027;
t2870 = t2892 * t3027 + t2903 * t3023;
t3064 = -t2870 * t3028 + t2884 * t3024;
t2823 = -t3017 * t2871 + t3020 * t3064;
t2845 = t2870 * t3024 + t2884 * t3028;
t3071 = t2823 * t3019 + t2845 * t3016;
t2893 = -t3121 - t3152;
t2907 = -t3129 - t3130;
t2873 = t2893 * t3027 - t2907 * t3023;
t2872 = t2893 * t3023 + t2907 * t3027;
t2886 = -t2952 * t3104 + t3052;
t3063 = -t2872 * t3028 + t2886 * t3024;
t2825 = -t3017 * t2873 + t3020 * t3063;
t2848 = t2872 * t3024 + t2886 * t3028;
t3070 = t2825 * t3019 + t2848 * t3016;
t2913 = t2942 + t3117;
t3058 = -t2913 * t3024 - t3153;
t2866 = t3020 * t3058 - t3142;
t2888 = -t2913 * t3028 + t3154;
t3066 = t2866 * t3019 + t2888 * t3016;
t2915 = -t2969 * t3106 + t3043;
t3057 = t2915 * t3024 - t3153;
t2867 = t3020 * t3057 - t3142;
t2889 = t2915 * t3028 + t3154;
t3065 = t2867 * t3019 + t2889 * t3016;
t2912 = t2942 - t3117;
t2876 = -t3017 * t2912 - t3148;
t3062 = t2876 * t3019 + t3151;
t2914 = -t2969 * t3107 + t3043;
t2877 = -t3017 * t2914 + t3148;
t3061 = t2877 * t3019 - t3151;
t2950 = -t3016 * t2989 + t3047;
t2951 = -g(3) * t3114 + t3094;
t3051 = t2950 * t3019 + t2951 * t3016;
t3000 = t3019 * t3095;
t2985 = t3000 - t3092;
t2999 = t3016 * t3095;
t2986 = t2999 + t3091;
t3050 = t2985 * t3019 + t2986 * t3016;
t2998 = t3019 * t3016 * t3115;
t2991 = t2998 + t3101;
t3014 = t3019 ^ 2;
t3015 = t3021 ^ 2;
t2995 = (-t3013 * t3014 - t3015) * t3030;
t3049 = t2991 * t3019 + t2995 * t3016;
t2992 = t2998 - t3101;
t3012 = t3016 ^ 2;
t2994 = (-t3012 * t3013 - t3015) * t3030;
t3048 = t2992 * t3016 + t2994 * t3019;
t3004 = -qJDD(1) * t3025 - t3029 * t3030;
t3003 = qJDD(1) * t3029 - t3025 * t3030;
t2993 = (-t3012 - t3014) * t3115;
t2987 = t2999 - t3091;
t2984 = t3000 + t3092;
t2971 = -t3018 * t2988 + t3090;
t2962 = t2992 * t3019 - t2994 * t3016;
t2961 = -t2991 * t3016 + t2995 * t3019;
t2958 = -t2985 * t3016 + t2986 * t3019;
t2948 = -t3018 * t2984 + t3021 * t3048;
t2947 = -t3018 * t2987 + t3021 * t3049;
t2946 = t3021 * t2984 + t3018 * t3048;
t2945 = t3021 * t2987 + t3018 * t3049;
t2941 = -t3018 * t2993 + t3021 * t3050;
t2940 = t3021 * t2993 + t3018 * t3050;
t2910 = -t2950 * t3016 + t2951 * t3019;
t2900 = -t3018 * t2971 + t3021 * t3051;
t2899 = t3021 * t2971 + t3018 * t3051;
t2875 = t3020 * t2914 + t3149;
t2874 = t3020 * t2912 - t3149;
t2865 = t3017 * t3057 + t3141;
t2864 = t3017 * t3058 + t3141;
t2855 = -t2877 * t3016 - t3150;
t2854 = -t2876 * t3016 + t3150;
t2844 = -t2867 * t3016 + t2889 * t3019;
t2843 = -t2866 * t3016 + t2888 * t3019;
t2829 = -t3018 * t2875 + t3021 * t3061;
t2828 = -t3018 * t2874 + t3021 * t3062;
t2827 = t3021 * t2875 + t3018 * t3061;
t2826 = t3021 * t2874 + t3018 * t3062;
t2824 = t3020 * t2873 + t3017 * t3063;
t2822 = t3020 * t2871 + t3017 * t3064;
t2821 = -t3018 * t2865 + t3021 * t3065;
t2820 = -t3018 * t2864 + t3021 * t3066;
t2819 = t3021 * t2865 + t3018 * t3065;
t2818 = t3021 * t2864 + t3018 * t3066;
t2816 = t3020 * t2890 + t3017 * t3067;
t2807 = t3020 * t2853 + t3017 * t3068;
t2801 = -t2825 * t3016 + t2848 * t3019;
t2800 = -t2823 * t3016 + t2845 * t3019;
t2798 = t3020 * t2851 + t3017 * t3069;
t2794 = -t2817 * t3016 + t2834 * t3019;
t2792 = -t2808 * t3016 + t2840 * t3019;
t2790 = -t3018 * t2824 + t3021 * t3070;
t2789 = t3021 * t2824 + t3018 * t3070;
t2788 = -t3018 * t2822 + t3021 * t3071;
t2787 = t3021 * t2822 + t3018 * t3071;
t2786 = -t2799 * t3016 + t2813 * t3019;
t2785 = -t3018 * t2816 + t3021 * t3072;
t2784 = t3021 * t2816 + t3018 * t3072;
t2782 = t3020 * t2812 + t3017 * t3073;
t2778 = t3020 * t2810 + t3017 * t3074;
t2775 = -t3018 * t2807 + t3021 * t3075;
t2774 = t3021 * t2807 + t3018 * t3075;
t2771 = t3020 * t2805 + t3017 * t3076;
t2770 = -t3018 * t2798 + t3021 * t3077;
t2769 = t3021 * t2798 + t3018 * t3077;
t2768 = -t2783 * t3016 + t2797 * t3019;
t2767 = -t2779 * t3016 + t2793 * t3019;
t2764 = -t2772 * t3016 + t2791 * t3019;
t2762 = t3020 * t2781 + t3017 * t3079;
t2761 = -t3018 * t2782 + t3021 * t3078;
t2760 = t3021 * t2782 + t3018 * t3078;
t2759 = -t3018 * t2778 + t3021 * t3080;
t2758 = t3021 * t2778 + t3018 * t3080;
t2755 = -t3018 * t2771 + t3021 * t3081;
t2754 = t3021 * t2771 + t3018 * t3081;
t2753 = -t2763 * t3016 + t2773 * t3019;
t2751 = -t3018 * t2762 + t3021 * t3082;
t2750 = t3021 * t2762 + t3018 * t3082;
t2748 = t3020 * t2757 + t3017 * t3083;
t2747 = -t2749 * t3016 + t2752 * t3019;
t2746 = -t3018 * t2748 + t3021 * t3084;
t2745 = t3021 * t2748 + t3018 * t3084;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t3004, -t3003, 0, -t3007 * t3025 + t3008 * t3029, 0, 0, 0, 0, 0, 0, -t2947 * t3025 + t2961 * t3029, -t2948 * t3025 + t2962 * t3029, -t2941 * t3025 + t2958 * t3029, -t2900 * t3025 + t2910 * t3029, 0, 0, 0, 0, 0, 0, -t2828 * t3025 + t2854 * t3029, t3160, -t2821 * t3025 + t2844 * t3029, -t2785 * t3025 + t2794 * t3029, 0, 0, 0, 0, 0, 0, -t2820 * t3025 + t2843 * t3029, -t2829 * t3025 + t2855 * t3029, -t3160, -t2770 * t3025 + t2786 * t3029, 0, 0, 0, 0, 0, 0, -t2788 * t3025 + t2800 * t3029, -t2790 * t3025 + t2801 * t3029, -t2775 * t3025 + t2792 * t3029, -t2751 * t3025 + t2753 * t3029, 0, 0, 0, 0, 0, 0, -t2759 * t3025 + t2767 * t3029, -t2761 * t3025 + t2768 * t3029, -t2755 * t3025 + t2764 * t3029, -t2746 * t3025 + t2747 * t3029; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t3003, t3004, 0, t3007 * t3029 + t3008 * t3025, 0, 0, 0, 0, 0, 0, t2947 * t3029 + t2961 * t3025, t2948 * t3029 + t2962 * t3025, t2941 * t3029 + t2958 * t3025, t2900 * t3029 + t2910 * t3025, 0, 0, 0, 0, 0, 0, t2828 * t3029 + t2854 * t3025, -t3159, t2821 * t3029 + t2844 * t3025, t2785 * t3029 + t2794 * t3025, 0, 0, 0, 0, 0, 0, t2820 * t3029 + t2843 * t3025, t2829 * t3029 + t2855 * t3025, t3159, t2770 * t3029 + t2786 * t3025, 0, 0, 0, 0, 0, 0, t2788 * t3029 + t2800 * t3025, t2790 * t3029 + t2801 * t3025, t2775 * t3029 + t2792 * t3025, t2751 * t3029 + t2753 * t3025, 0, 0, 0, 0, 0, 0, t2759 * t3029 + t2767 * t3025, t2761 * t3029 + t2768 * t3025, t2755 * t3029 + t2764 * t3025, t2746 * t3029 + t2747 * t3025; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2945, t2946, t2940, t2899, 0, 0, 0, 0, 0, 0, t2826, -t2831, t2819, t2784, 0, 0, 0, 0, 0, 0, t2818, t2827, t2831, t2769, 0, 0, 0, 0, 0, 0, t2787, t2789, t2774, t2750, 0, 0, 0, 0, 0, 0, t2758, t2760, t2754, t2745; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3030, -qJDD(1), 0, t3008, 0, 0, 0, 0, 0, 0, t2961, t2962, t2958, t2910, 0, 0, 0, 0, 0, 0, t2854, t2860, t2844, t2794, 0, 0, 0, 0, 0, 0, t2843, t2855, -t2860, t2786, 0, 0, 0, 0, 0, 0, t2800, t2801, t2792, t2753, 0, 0, 0, 0, 0, 0, t2767, t2768, t2764, t2747; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t3030, 0, t3007, 0, 0, 0, 0, 0, 0, t2947, t2948, t2941, t2900, 0, 0, 0, 0, 0, 0, t2828, -t2833, t2821, t2785, 0, 0, 0, 0, 0, 0, t2820, t2829, t2833, t2770, 0, 0, 0, 0, 0, 0, t2788, t2790, t2775, t2751, 0, 0, 0, 0, 0, 0, t2759, t2761, t2755, t2746; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2945, t2946, t2940, t2899, 0, 0, 0, 0, 0, 0, t2826, -t2831, t2819, t2784, 0, 0, 0, 0, 0, 0, t2818, t2827, t2831, t2769, 0, 0, 0, 0, 0, 0, t2787, t2789, t2774, t2750, 0, 0, 0, 0, 0, 0, t2758, t2760, t2754, t2745; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2995, t2992, t2986, t2951, 0, 0, 0, 0, 0, 0, t2894, -t2898, t2889, t2834, 0, 0, 0, 0, 0, 0, t2888, -t2894, t2898, t2813, 0, 0, 0, 0, 0, 0, t2845, t2848, t2840, t2773, 0, 0, 0, 0, 0, 0, t2793, t2797, t2791, t2752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2991, t2994, t2985, t2950, 0, 0, 0, 0, 0, 0, t2876, -t2882, t2867, t2817, 0, 0, 0, 0, 0, 0, t2866, t2877, t2882, t2799, 0, 0, 0, 0, 0, 0, t2823, t2825, t2808, t2763, 0, 0, 0, 0, 0, 0, t2779, t2783, t2772, t2749; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2987, t2984, t2993, t2971, 0, 0, 0, 0, 0, 0, t2874, t2879, t2865, t2816, 0, 0, 0, 0, 0, 0, t2864, t2875, -t2879, t2798, 0, 0, 0, 0, 0, 0, t2822, t2824, t2807, t2762, 0, 0, 0, 0, 0, 0, t2778, t2782, t2771, t2748; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2929, t3039, t2915, t2863, 0, 0, 0, 0, 0, 0, -t2913, -t2929, -t3039, t2849, 0, 0, 0, 0, 0, 0, t2884, t2886, t2896, t2837, 0, 0, 0, 0, 0, 0, t2841, t2846, t2814, t2765; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2932, t3093, -t2916, -t3085, 0, 0, 0, 0, 0, 0, -t2916, t2932, -t3093, -t2850, 0, 0, 0, 0, 0, 0, -t2870, -t2872, -t2852, -t2780, 0, 0, 0, 0, 0, 0, -t2809, -t2811, -t2804, -t2756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2912, t2917, t3136, t2890, 0, 0, 0, 0, 0, 0, t3136, t2914, -t2917, t2851, 0, 0, 0, 0, 0, 0, t2871, t2873, t2853, t2781, 0, 0, 0, 0, 0, 0, t2810, t2812, t2805, t2757; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3136, t2914, -t2917, t2851, 0, 0, 0, 0, 0, 0, t2871, t2873, t2853, t2781, 0, 0, 0, 0, 0, 0, t2810, t2812, t2805, t2757; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2913, t2929, t3039, -t2849, 0, 0, 0, 0, 0, 0, -t2884, -t2886, -t2896, -t2837, 0, 0, 0, 0, 0, 0, -t2841, -t2846, -t2814, -t2765; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2916, -t2932, t3093, t2850, 0, 0, 0, 0, 0, 0, t2870, t2872, t2852, t2780, 0, 0, 0, 0, 0, 0, t2809, t2811, t2804, t2756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2903, t2893, t2885, t2803, 0, 0, 0, 0, 0, 0, t2842, t2847, t2815, t2766; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2892, t2907, t2887, t2802, 0, 0, 0, 0, 0, 0, -t2856, -t2858, -t2878, -t2795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2884, t2886, t2896, t2837, 0, 0, 0, 0, 0, 0, t2841, t2846, t2814, t2765; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2883, t2869, t2857, t2777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2868, t2891, t2859, t2776; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2856, t2858, t2878, t2795;];
f_new_reg  = t1;
