% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRPPR9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 06:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRPPR9_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR9_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:43:09
% EndTime: 2019-05-07 06:43:24
% DurationCPUTime: 15.72s
% Computational Cost: add. (67562->352), mult. (148889->481), div. (0->0), fcn. (118832->12), ass. (0->269)
t2920 = cos(pkin(6));
t2912 = qJD(1) * t2920 + qJD(2);
t2922 = sin(qJ(3));
t2926 = cos(qJ(3));
t2918 = sin(pkin(6));
t2923 = sin(qJ(2));
t2981 = t2918 * t2923;
t2973 = qJD(1) * t2981;
t2886 = t2912 * t2922 + t2926 * t2973;
t2927 = cos(qJ(2));
t2980 = t2918 * t2927;
t2972 = qJD(1) * t2980;
t2975 = qJDD(1) * t2918;
t2893 = qJD(2) * t2972 + t2923 * t2975;
t2911 = t2920 * qJDD(1) + qJDD(2);
t2967 = t2922 * t2893 - t2926 * t2911;
t2849 = qJD(3) * t2886 + t2967;
t2903 = -qJD(3) + t2972;
t2917 = sin(pkin(11));
t2919 = cos(pkin(11));
t2865 = t2917 * t2886 + t2903 * t2919;
t2867 = t2886 * t2919 - t2903 * t2917;
t2986 = t2865 * t2867;
t2797 = t2849 + t2986;
t2864 = t2867 ^ 2;
t2884 = -t2926 * t2912 + t2922 * t2973;
t2999 = t2884 ^ 2;
t3007 = -t2864 - t2999;
t2773 = t2797 * t2919 + t2917 * t3007;
t2943 = -t2926 * t2893 - t2922 * t2911;
t2850 = -qJD(3) * t2884 - t2943;
t2906 = qJD(2) * t2973;
t2974 = qJDD(1) * t2927;
t2964 = t2918 * t2974 - t2906;
t2942 = qJDD(3) - t2964;
t2936 = t2919 * t2850 + t2917 * t2942;
t2985 = t2865 * t2884;
t2934 = t2936 - t2985;
t2742 = t2773 * t2922 + t2926 * t2934;
t2744 = t2773 * t2926 - t2922 * t2934;
t2771 = t2797 * t2917 - t2919 * t3007;
t2954 = t2744 * t2923 - t2771 * t2927;
t2705 = -t2918 * t2742 + t2920 * t2954;
t2724 = t2744 * t2927 + t2771 * t2923;
t2924 = sin(qJ(1));
t2928 = cos(qJ(1));
t3052 = t2705 * t2924 - t2724 * t2928;
t3051 = t2705 * t2928 + t2724 * t2924;
t2703 = t2920 * t2742 + t2918 * t2954;
t2793 = t2936 + t2985;
t2970 = t2850 * t2917 - t2919 * t2942;
t2984 = t2884 * t2867;
t2941 = -t2970 + t2984;
t3005 = -t2793 * t2919 + t2917 * t2941;
t3001 = t2865 ^ 2;
t2805 = t2864 + t3001;
t3004 = t2793 * t2917 + t2919 * t2941;
t3022 = -t2805 * t2922 + t2926 * t3004;
t3026 = t2923 * t3005 + t2927 * t3022;
t3021 = t2805 * t2926 + t2922 * t3004;
t3028 = t2923 * t3022 - t2927 * t3005;
t3035 = -t2918 * t3021 + t2920 * t3028;
t3046 = -t2924 * t3035 + t2928 * t3026;
t2798 = t2849 - t2986;
t3006 = -t2999 - t3001;
t3014 = t2798 * t2919 + t2917 * t3006;
t2940 = t2970 + t2984;
t3013 = -t2798 * t2917 + t2919 * t3006;
t3019 = t2922 * t2940 + t2926 * t3013;
t3025 = t2923 * t3014 + t2927 * t3019;
t3020 = t2922 * t3013 - t2926 * t2940;
t3027 = t2923 * t3019 - t2927 * t3014;
t3036 = -t2918 * t3020 + t2920 * t3027;
t3045 = -t2924 * t3036 + t2928 * t3025;
t3044 = t2924 * t3026 + t2928 * t3035;
t3043 = t2924 * t3025 + t2928 * t3036;
t3038 = t2918 * t3027 + t2920 * t3020;
t3037 = t2918 * t3028 + t2920 * t3021;
t2880 = -qJD(6) + t2884;
t3008 = qJD(6) - t2880;
t2921 = sin(qJ(6));
t2925 = cos(qJ(6));
t2827 = -t2925 * t2865 + t2867 * t2921;
t3003 = t2827 ^ 2;
t2829 = t2865 * t2921 + t2867 * t2925;
t3002 = t2829 ^ 2;
t3000 = t2880 ^ 2;
t2998 = t2886 ^ 2;
t2997 = t2903 ^ 2;
t2996 = t2912 ^ 2;
t2995 = -2 * qJD(4);
t2994 = t2920 * g(3);
t2993 = qJD(1) * t2923;
t2992 = qJD(1) * t2927;
t2987 = t2827 * t2829;
t2983 = t2884 * t2886;
t2929 = qJD(1) ^ 2;
t2982 = t2918 ^ 2 * t2929;
t2979 = qJD(3) + t2903;
t2836 = pkin(4) * t2865 - qJ(5) * t2867;
t2978 = (2 * qJD(4)) + t2836;
t2977 = qJD(6) + t2880;
t2892 = (-pkin(2) * t2927 - pkin(9) * t2923) * t2918 * qJD(1);
t2905 = -g(1) * t2928 - g(2) * t2924;
t2889 = -pkin(1) * t2929 + pkin(8) * t2975 + t2905;
t2904 = t2924 * g(1) - t2928 * g(2);
t2937 = t2929 * t2918 * pkin(8) + qJDD(1) * pkin(1) + t2904;
t2935 = t2920 * t2937;
t2976 = t2927 * t2889 + t2923 * t2935;
t2819 = -t2996 * pkin(2) + t2911 * pkin(9) + (-g(3) * t2923 + t2892 * t2992) * t2918 + t2976;
t2930 = t2906 * pkin(2) - t2893 * pkin(9) - t2994 + (-t2912 * pkin(9) * t2992 + (t2912 * t2993 - t2974) * pkin(2) - t2937) * t2918;
t2784 = t2926 * t2819 + t2922 * t2930;
t2855 = pkin(3) * t2884 - qJ(4) * t2886;
t2762 = -pkin(3) * t2997 + qJ(4) * t2942 - t2884 * t2855 + t2784;
t2968 = t2923 * t2889 - t2927 * t2935;
t2818 = -t2911 * pkin(2) - t2996 * pkin(9) + (g(3) * t2927 + t2892 * t2993) * t2918 + t2968;
t2969 = -t2884 * t2903 - t2850;
t2931 = t2969 * qJ(4) + (-t2886 * t2903 + t2849) * pkin(3) + t2818;
t2728 = t2919 * t2762 + t2865 * t2995 + t2917 * t2931;
t2971 = t2917 * t2762 - t2919 * t2931;
t2783 = -t2922 * t2819 + t2926 * t2930;
t2966 = t2912 * t2972;
t2965 = -pkin(5) * t2884 - pkin(10) * t2867;
t2939 = -t2849 * pkin(4) - qJ(5) * t2999 + qJDD(5) + t2971;
t2701 = -t2849 * pkin(5) - t2793 * pkin(10) + (pkin(5) * t2865 + t2978) * t2867 + t2939;
t2716 = -pkin(4) * t2999 + t2849 * qJ(5) + 0.2e1 * qJD(5) * t2884 - t2865 * t2836 + t2728;
t2702 = -pkin(5) * t3001 + pkin(10) * t2970 + t2884 * t2965 + t2716;
t2679 = t2701 * t2925 - t2702 * t2921;
t2680 = t2701 * t2921 + t2702 * t2925;
t2668 = t2679 * t2925 + t2680 * t2921;
t2669 = -t2679 * t2921 + t2680 * t2925;
t2661 = t2668 * t2917 + t2669 * t2919;
t2761 = -t2942 * pkin(3) - t2997 * qJ(4) + t2886 * t2855 + qJDD(4) - t2783;
t2726 = t2970 * pkin(4) + (pkin(4) * t2884 - 0.2e1 * qJD(5)) * t2867 + t2761 - t2934 * qJ(5);
t2715 = pkin(5) * t2970 + pkin(10) * t3001 - t2867 * t2965 + t2726;
t2659 = t2661 * t2926 + t2715 * t2922;
t2660 = -t2668 * t2919 + t2669 * t2917;
t2963 = t2659 * t2923 - t2660 * t2927;
t2717 = t2978 * t2867 + t2939;
t2686 = t2716 * t2919 + t2717 * t2917;
t2678 = t2686 * t2926 + t2726 * t2922;
t2685 = t2716 * t2917 - t2717 * t2919;
t2962 = t2678 * t2923 - t2685 * t2927;
t2933 = -t2921 * t2936 + t2925 * t2970;
t2751 = -t2829 * t2977 + t2933;
t2932 = -t2921 * t2970 - t2925 * t2936;
t2753 = t2827 * t2977 + t2932;
t2720 = t2751 * t2921 + t2753 * t2925;
t2721 = t2751 * t2925 - t2753 * t2921;
t2690 = t2720 * t2917 + t2721 * t2919;
t2769 = -t3002 - t3003;
t2684 = t2690 * t2926 - t2769 * t2922;
t2689 = -t2720 * t2919 + t2721 * t2917;
t2961 = t2684 * t2923 - t2689 * t2927;
t2727 = t2867 * t2995 - t2971;
t2692 = -t2727 * t2917 + t2728 * t2919;
t2688 = t2692 * t2926 + t2761 * t2922;
t2691 = t2727 * t2919 + t2728 * t2917;
t2960 = t2688 * t2923 - t2691 * t2927;
t2938 = qJDD(6) - t2849;
t2776 = t2938 - t2987;
t2785 = -t3000 - t3003;
t2734 = t2776 * t2925 + t2785 * t2921;
t2735 = -t2776 * t2921 + t2785 * t2925;
t2714 = t2734 * t2917 + t2735 * t2919;
t2750 = t2829 * t3008 - t2933;
t2694 = t2714 * t2926 - t2750 * t2922;
t2713 = -t2734 * t2919 + t2735 * t2917;
t2959 = t2694 * t2923 - t2713 * t2927;
t2775 = -t2938 - t2987;
t2800 = -t3000 - t3002;
t2754 = t2775 * t2921 + t2800 * t2925;
t2755 = t2775 * t2925 - t2800 * t2921;
t2719 = t2754 * t2917 + t2755 * t2919;
t2752 = -t2827 * t3008 - t2932;
t2700 = t2719 * t2926 - t2752 * t2922;
t2718 = -t2754 * t2919 + t2755 * t2917;
t2958 = t2700 * t2923 - t2718 * t2927;
t2741 = -t2783 * t2922 + t2784 * t2926;
t2955 = t2741 * t2923 - t2818 * t2927;
t2823 = -t2886 * t2979 - t2967;
t2825 = t2884 * t2979 + t2943;
t2787 = t2823 * t2926 - t2825 * t2922;
t2839 = -t2998 - t2999;
t2950 = t2787 * t2923 - t2839 * t2927;
t2841 = t2942 - t2983;
t2851 = -t2999 - t2997;
t2802 = -t2841 * t2922 + t2851 * t2926;
t2822 = (qJD(3) - t2903) * t2886 + t2967;
t2949 = t2802 * t2923 - t2822 * t2927;
t2840 = -t2942 - t2983;
t2857 = -t2997 - t2998;
t2808 = t2840 * t2926 - t2857 * t2922;
t2948 = t2808 * t2923 + t2927 * t2969;
t2852 = -g(3) * t2980 - t2968;
t2853 = -g(3) * t2981 + t2976;
t2947 = t2852 * t2927 + t2853 * t2923;
t2869 = t2966 - t2893;
t2896 = t2912 * t2973;
t2870 = t2896 + t2964;
t2946 = t2869 * t2927 + t2870 * t2923;
t2915 = t2923 ^ 2;
t2878 = -t2915 * t2982 - t2996;
t2902 = t2927 * t2923 * t2982;
t2891 = t2902 - t2911;
t2945 = t2878 * t2927 + t2891 * t2923;
t2890 = t2902 + t2911;
t2916 = t2927 ^ 2;
t2894 = -t2916 * t2982 - t2996;
t2944 = t2890 * t2927 + t2894 * t2923;
t2901 = -qJDD(1) * t2924 - t2928 * t2929;
t2900 = qJDD(1) * t2928 - t2924 * t2929;
t2895 = (-t2915 - t2916) * t2982;
t2873 = -t2918 * t2937 - t2994;
t2871 = t2896 - t2964;
t2868 = t2966 + t2893;
t2862 = -t2890 * t2923 + t2894 * t2927;
t2854 = -t2878 * t2923 + t2891 * t2927;
t2838 = -t2869 * t2923 + t2870 * t2927;
t2835 = -t2918 * t2871 + t2920 * t2944;
t2834 = t2920 * t2871 + t2918 * t2944;
t2821 = -t2918 * t2868 + t2920 * t2945;
t2820 = t2920 * t2868 + t2918 * t2945;
t2817 = -t2918 * t2895 + t2920 * t2946;
t2816 = t2920 * t2895 + t2918 * t2946;
t2807 = t2840 * t2922 + t2857 * t2926;
t2806 = -t2852 * t2923 + t2853 * t2927;
t2801 = t2841 * t2926 + t2851 * t2922;
t2789 = -t2918 * t2873 + t2920 * t2947;
t2788 = t2920 * t2873 + t2918 * t2947;
t2786 = t2823 * t2922 + t2825 * t2926;
t2782 = t2808 * t2927 - t2923 * t2969;
t2781 = t2802 * t2927 + t2822 * t2923;
t2770 = t2787 * t2927 + t2839 * t2923;
t2759 = -t2918 * t2807 + t2920 * t2948;
t2758 = t2920 * t2807 + t2918 * t2948;
t2757 = -t2918 * t2801 + t2920 * t2949;
t2756 = t2920 * t2801 + t2918 * t2949;
t2740 = t2783 * t2926 + t2784 * t2922;
t2733 = -t2918 * t2786 + t2920 * t2950;
t2732 = t2920 * t2786 + t2918 * t2950;
t2731 = t2741 * t2927 + t2818 * t2923;
t2712 = -t2918 * t2740 + t2920 * t2955;
t2711 = t2920 * t2740 + t2918 * t2955;
t2699 = t2719 * t2922 + t2752 * t2926;
t2693 = t2714 * t2922 + t2750 * t2926;
t2687 = t2692 * t2922 - t2761 * t2926;
t2683 = t2690 * t2922 + t2769 * t2926;
t2682 = t2700 * t2927 + t2718 * t2923;
t2681 = t2694 * t2927 + t2713 * t2923;
t2677 = t2686 * t2922 - t2726 * t2926;
t2676 = t2688 * t2927 + t2691 * t2923;
t2675 = -t2918 * t2699 + t2920 * t2958;
t2674 = t2920 * t2699 + t2918 * t2958;
t2673 = t2684 * t2927 + t2689 * t2923;
t2672 = -t2918 * t2693 + t2920 * t2959;
t2671 = t2920 * t2693 + t2918 * t2959;
t2670 = t2678 * t2927 + t2685 * t2923;
t2667 = -t2918 * t2687 + t2920 * t2960;
t2666 = t2920 * t2687 + t2918 * t2960;
t2665 = -t2918 * t2683 + t2920 * t2961;
t2664 = t2920 * t2683 + t2918 * t2961;
t2663 = -t2918 * t2677 + t2920 * t2962;
t2662 = t2920 * t2677 + t2918 * t2962;
t2658 = t2661 * t2922 - t2715 * t2926;
t2657 = t2659 * t2927 + t2660 * t2923;
t2656 = -t2918 * t2658 + t2920 * t2963;
t2655 = t2920 * t2658 + t2918 * t2963;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2901, -t2900, 0, -t2904 * t2924 + t2905 * t2928, 0, 0, 0, 0, 0, 0, -t2835 * t2924 + t2862 * t2928, -t2821 * t2924 + t2854 * t2928, -t2817 * t2924 + t2838 * t2928, -t2789 * t2924 + t2806 * t2928, 0, 0, 0, 0, 0, 0, -t2757 * t2924 + t2781 * t2928, -t2759 * t2924 + t2782 * t2928, -t2733 * t2924 + t2770 * t2928, -t2712 * t2924 + t2731 * t2928, 0, 0, 0, 0, 0, 0, t3045, t3052, t3046, -t2667 * t2924 + t2676 * t2928, 0, 0, 0, 0, 0, 0, t3045, t3046, -t3052, -t2663 * t2924 + t2670 * t2928, 0, 0, 0, 0, 0, 0, -t2672 * t2924 + t2681 * t2928, -t2675 * t2924 + t2682 * t2928, -t2665 * t2924 + t2673 * t2928, -t2656 * t2924 + t2657 * t2928; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2900, t2901, 0, t2904 * t2928 + t2905 * t2924, 0, 0, 0, 0, 0, 0, t2835 * t2928 + t2862 * t2924, t2821 * t2928 + t2854 * t2924, t2817 * t2928 + t2838 * t2924, t2789 * t2928 + t2806 * t2924, 0, 0, 0, 0, 0, 0, t2757 * t2928 + t2781 * t2924, t2759 * t2928 + t2782 * t2924, t2733 * t2928 + t2770 * t2924, t2712 * t2928 + t2731 * t2924, 0, 0, 0, 0, 0, 0, t3043, -t3051, t3044, t2667 * t2928 + t2676 * t2924, 0, 0, 0, 0, 0, 0, t3043, t3044, t3051, t2663 * t2928 + t2670 * t2924, 0, 0, 0, 0, 0, 0, t2672 * t2928 + t2681 * t2924, t2675 * t2928 + t2682 * t2924, t2665 * t2928 + t2673 * t2924, t2656 * t2928 + t2657 * t2924; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2834, t2820, t2816, t2788, 0, 0, 0, 0, 0, 0, t2756, t2758, t2732, t2711, 0, 0, 0, 0, 0, 0, t3038, -t2703, t3037, t2666, 0, 0, 0, 0, 0, 0, t3038, t3037, t2703, t2662, 0, 0, 0, 0, 0, 0, t2671, t2674, t2664, t2655; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2929, -qJDD(1), 0, t2905, 0, 0, 0, 0, 0, 0, t2862, t2854, t2838, t2806, 0, 0, 0, 0, 0, 0, t2781, t2782, t2770, t2731, 0, 0, 0, 0, 0, 0, t3025, -t2724, t3026, t2676, 0, 0, 0, 0, 0, 0, t3025, t3026, t2724, t2670, 0, 0, 0, 0, 0, 0, t2681, t2682, t2673, t2657; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2929, 0, t2904, 0, 0, 0, 0, 0, 0, t2835, t2821, t2817, t2789, 0, 0, 0, 0, 0, 0, t2757, t2759, t2733, t2712, 0, 0, 0, 0, 0, 0, t3036, -t2705, t3035, t2667, 0, 0, 0, 0, 0, 0, t3036, t3035, t2705, t2663, 0, 0, 0, 0, 0, 0, t2672, t2675, t2665, t2656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2834, t2820, t2816, t2788, 0, 0, 0, 0, 0, 0, t2756, t2758, t2732, t2711, 0, 0, 0, 0, 0, 0, t3038, -t2703, t3037, t2666, 0, 0, 0, 0, 0, 0, t3038, t3037, t2703, t2662, 0, 0, 0, 0, 0, 0, t2671, t2674, t2664, t2655; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2894, t2891, t2870, t2853, 0, 0, 0, 0, 0, 0, t2802, t2808, t2787, t2741, 0, 0, 0, 0, 0, 0, t3019, -t2744, t3022, t2688, 0, 0, 0, 0, 0, 0, t3019, t3022, t2744, t2678, 0, 0, 0, 0, 0, 0, t2694, t2700, t2684, t2659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2890, t2878, t2869, t2852, 0, 0, 0, 0, 0, 0, -t2822, t2969, -t2839, -t2818, 0, 0, 0, 0, 0, 0, -t3014, t2771, -t3005, -t2691, 0, 0, 0, 0, 0, 0, -t3014, -t3005, -t2771, -t2685, 0, 0, 0, 0, 0, 0, -t2713, -t2718, -t2689, -t2660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2871, t2868, t2895, t2873, 0, 0, 0, 0, 0, 0, t2801, t2807, t2786, t2740, 0, 0, 0, 0, 0, 0, t3020, -t2742, t3021, t2687, 0, 0, 0, 0, 0, 0, t3020, t3021, t2742, t2677, 0, 0, 0, 0, 0, 0, t2693, t2699, t2683, t2658; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2851, t2840, t2823, t2784, 0, 0, 0, 0, 0, 0, t3013, -t2773, t3004, t2692, 0, 0, 0, 0, 0, 0, t3013, t3004, t2773, t2686, 0, 0, 0, 0, 0, 0, t2714, t2719, t2690, t2661; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2841, t2857, t2825, t2783, 0, 0, 0, 0, 0, 0, -t2940, -t2934, t2805, -t2761, 0, 0, 0, 0, 0, 0, -t2940, t2805, t2934, -t2726, 0, 0, 0, 0, 0, 0, t2750, t2752, t2769, -t2715; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2822, -t2969, t2839, t2818, 0, 0, 0, 0, 0, 0, t3014, -t2771, t3005, t2691, 0, 0, 0, 0, 0, 0, t3014, t3005, t2771, t2685, 0, 0, 0, 0, 0, 0, t2713, t2718, t2689, t2660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3006, -t2797, t2941, t2728, 0, 0, 0, 0, 0, 0, t3006, t2941, t2797, t2716, 0, 0, 0, 0, 0, 0, t2735, t2755, t2721, t2669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2798, t3007, -t2793, t2727, 0, 0, 0, 0, 0, 0, t2798, -t2793, -t3007, -t2717, 0, 0, 0, 0, 0, 0, -t2734, -t2754, -t2720, -t2668; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2940, t2934, -t2805, t2761, 0, 0, 0, 0, 0, 0, t2940, -t2805, -t2934, t2726, 0, 0, 0, 0, 0, 0, -t2750, -t2752, -t2769, t2715; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3006, t2941, t2797, t2716, 0, 0, 0, 0, 0, 0, t2735, t2755, t2721, t2669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2940, -t2805, -t2934, t2726, 0, 0, 0, 0, 0, 0, -t2750, -t2752, -t2769, t2715; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2798, t2793, t3007, t2717, 0, 0, 0, 0, 0, 0, t2734, t2754, t2720, t2668; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2785, t2775, t2751, t2680; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2776, t2800, t2753, t2679; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2750, t2752, t2769, -t2715;];
f_new_reg  = t1;
