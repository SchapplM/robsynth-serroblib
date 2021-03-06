% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 09:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRPRP12_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:37:49
% EndTime: 2019-05-07 09:38:02
% DurationCPUTime: 13.04s
% Computational Cost: add. (32606->331), mult. (70726->409), div. (0->0), fcn. (55260->10), ass. (0->245)
t2883 = sin(qJ(3));
t2887 = cos(qJ(3));
t2881 = cos(pkin(6));
t2938 = qJD(1) * t2881 + qJD(2);
t2880 = sin(pkin(6));
t2884 = sin(qJ(2));
t2953 = t2880 * t2884;
t2943 = qJD(1) * t2953;
t2842 = t2883 * t2943 - t2887 * t2938;
t2888 = cos(qJ(2));
t2952 = t2880 * t2888;
t2942 = qJD(1) * t2952;
t2946 = qJDD(1) * t2880;
t2853 = qJD(2) * t2942 + t2884 * t2946;
t2875 = t2881 * qJDD(1) + qJDD(2);
t2909 = t2887 * t2853 + t2883 * t2875;
t2903 = -qJD(3) * t2842 + t2909;
t2901 = -qJDD(5) - t2903;
t2866 = -qJD(3) + t2942;
t2882 = sin(qJ(5));
t2886 = cos(qJ(5));
t2821 = -t2886 * t2842 - t2866 * t2882;
t2823 = t2842 * t2882 - t2866 * t2886;
t2957 = t2821 * t2823;
t2763 = t2901 - t2957;
t2844 = t2883 * t2938 + t2887 * t2943;
t2840 = qJD(5) + t2844;
t2837 = t2840 ^ 2;
t2964 = t2823 ^ 2;
t2944 = t2837 + t2964;
t2734 = -t2882 * t2763 + t2886 * t2944;
t2935 = t2883 * t2853 - t2887 * t2875;
t2812 = qJD(3) * t2844 + t2935;
t2870 = qJD(2) * t2943;
t2945 = qJDD(1) * t2888;
t2933 = t2880 * t2945 - t2870;
t2908 = -qJDD(3) + t2933;
t2914 = t2882 * t2812 - t2886 * t2908;
t2948 = qJD(5) + t2840;
t2752 = -t2821 * t2948 + t2914;
t2712 = t2734 * t2887 + t2752 * t2883;
t2710 = t2734 * t2883 - t2752 * t2887;
t2740 = t2763 * t2886 + t2882 * t2944;
t2924 = t2710 * t2884 + t2740 * t2888;
t2673 = t2880 * t2712 + t2881 * t2924;
t2688 = t2710 * t2888 - t2740 * t2884;
t2885 = sin(qJ(1));
t2889 = cos(qJ(1));
t3029 = t2673 * t2885 - t2688 * t2889;
t3028 = t2673 * t2889 + t2688 * t2885;
t2671 = -t2881 * t2712 + t2880 * t2924;
t2949 = qJD(5) - t2840;
t2900 = t2821 * t2949 - t2914;
t2937 = t2886 * t2812 + t2882 * t2908;
t2905 = -t2823 * t2949 + t2937;
t2967 = -t2882 * t2900 + t2886 * t2905;
t2798 = t2821 ^ 2;
t2768 = -t2964 - t2798;
t2965 = t2882 * t2905 + t2886 * t2900;
t2988 = t2768 * t2887 + t2883 * t2965;
t2997 = t2884 * t2967 + t2888 * t2988;
t2989 = t2768 * t2883 - t2887 * t2965;
t2999 = t2884 * t2988 - t2888 * t2967;
t3010 = -t2880 * t2989 + t2881 * t2999;
t3023 = -t2885 * t3010 + t2889 * t2997;
t2894 = t2901 + t2957;
t2969 = -t2798 - t2837;
t2976 = t2882 * t2894 + t2886 * t2969;
t2904 = t2823 * t2948 - t2937;
t2975 = t2882 * t2969 - t2886 * t2894;
t2987 = t2883 * t2975 + t2887 * t2904;
t2996 = t2884 * t2976 + t2888 * t2987;
t2986 = t2883 * t2904 - t2887 * t2975;
t2998 = t2884 * t2987 - t2888 * t2976;
t3011 = -t2880 * t2986 + t2881 * t2998;
t3022 = -t2885 * t3011 + t2889 * t2996;
t3021 = t2885 * t2997 + t2889 * t3010;
t3020 = t2885 * t2996 + t2889 * t3011;
t3013 = t2880 * t2998 + t2881 * t2986;
t3012 = t2880 * t2999 + t2881 * t2989;
t2956 = t2842 * t2844;
t2804 = t2908 + t2956;
t2841 = t2842 ^ 2;
t2862 = t2866 ^ 2;
t2813 = -t2862 - t2841;
t2766 = t2804 * t2883 + t2813 * t2887;
t3007 = t2766 * t2884;
t3006 = t2766 * t2888;
t2765 = t2804 * t2887 - t2813 * t2883;
t3005 = t2880 * t2765;
t3002 = t2881 * t2765;
t2902 = t2908 - t2956;
t2963 = t2844 ^ 2;
t2940 = -t2862 - t2963;
t2774 = t2883 * t2940 - t2887 * t2902;
t2995 = t2774 * t2884;
t2994 = t2774 * t2888;
t2771 = t2883 * t2902 + t2887 * t2940;
t2993 = t2880 * t2771;
t2992 = t2881 * t2771;
t2890 = qJD(1) ^ 2;
t2984 = t2880 * t2890;
t2968 = -t2841 - t2963;
t2980 = t2884 * t2968;
t2977 = t2888 * t2968;
t2791 = t2842 * t2866 + t2903;
t2934 = t2938 ^ 2;
t2950 = qJD(3) + t2866;
t2899 = t2842 * t2950 - t2909;
t2973 = t2883 * t2899;
t2971 = t2887 * t2899;
t2962 = -2 * qJD(4);
t2961 = t2881 * g(3);
t2960 = (-pkin(2) * t2888 - pkin(9) * t2884) * t2984;
t2954 = t2880 ^ 2 * t2890;
t2951 = qJD(3) - t2866;
t2869 = -g(1) * t2889 - g(2) * t2885;
t2849 = -pkin(1) * t2890 + pkin(8) * t2946 + t2869;
t2868 = t2885 * g(1) - t2889 * g(2);
t2898 = qJDD(1) * pkin(1) + pkin(8) * t2984 + t2868;
t2895 = t2881 * t2898;
t2947 = t2888 * t2849 + t2884 * t2895;
t2781 = t2875 * pkin(9) - t2934 * pkin(2) + (-g(3) * t2884 + t2888 * t2960) * t2880 + t2947;
t2932 = qJD(1) * t2938;
t2906 = t2884 * t2932;
t2907 = t2888 * t2932;
t2782 = t2870 * pkin(2) - t2853 * pkin(9) - t2961 + (-pkin(9) * t2907 + (t2906 - t2945) * pkin(2) - t2898) * t2880;
t2746 = -t2883 * t2781 + t2887 * t2782;
t2817 = pkin(3) * t2842 - qJ(4) * t2844;
t2730 = pkin(3) * t2908 - t2862 * qJ(4) + t2844 * t2817 + qJDD(4) - t2746;
t2699 = -t2899 * pkin(4) + pkin(10) * t2804 + t2730;
t2824 = pkin(4) * t2844 + pkin(10) * t2866;
t2936 = t2884 * t2849 - t2888 * t2895;
t2780 = -t2875 * pkin(2) - t2934 * pkin(9) + (g(3) * t2888 + t2884 * t2960) * t2880 + t2936;
t2891 = t2812 * pkin(3) - qJ(4) * t2791 + t2780;
t2939 = -pkin(3) * t2866 + t2962;
t2706 = -t2841 * pkin(4) + t2812 * pkin(10) + (-t2824 + t2939) * t2844 + t2891;
t2683 = t2882 * t2699 + t2886 * t2706;
t2797 = pkin(5) * t2821 - qJ(6) * t2823;
t2669 = -pkin(5) * t2837 - qJ(6) * t2901 + 0.2e1 * qJD(6) * t2840 - t2797 * t2821 + t2683;
t2926 = t2699 * t2886 - t2706 * t2882;
t2892 = -pkin(5) * t2901 + qJ(6) * t2837 - t2797 * t2823 - qJDD(6) + t2926;
t2659 = t2882 * t2669 + t2886 * t2892;
t2747 = t2887 * t2781 + t2883 * t2782;
t2896 = -t2862 * pkin(3) - qJ(4) * t2908 - t2842 * t2817 + t2747;
t2705 = -t2812 * pkin(4) - t2841 * pkin(10) + (t2962 - t2824) * t2866 + t2896;
t2686 = t2904 * pkin(5) - qJ(6) * t2752 - 0.2e1 * qJD(6) * t2823 + t2705;
t2656 = t2659 * t2883 + t2686 * t2887;
t2660 = t2669 * t2886 - t2882 * t2892;
t2931 = t2656 * t2884 - t2660 * t2888;
t2661 = t2882 * t2683 + t2886 * t2926;
t2658 = t2661 * t2883 + t2705 * t2887;
t2662 = t2683 * t2886 - t2882 * t2926;
t2930 = t2658 * t2884 - t2662 * t2888;
t2729 = t2866 * t2962 + t2896;
t2692 = t2729 * t2887 + t2730 * t2883;
t2731 = t2844 * t2939 + t2891;
t2929 = t2692 * t2884 - t2731 * t2888;
t2716 = -t2746 * t2883 + t2747 * t2887;
t2921 = t2716 * t2884 - t2780 * t2888;
t2788 = -t2844 * t2950 - t2935;
t2757 = t2788 * t2887 - t2973;
t2920 = t2757 * t2884 - t2977;
t2831 = t2866 * t2844;
t2789 = -t2812 - t2831;
t2758 = t2789 * t2887 - t2973;
t2919 = t2758 * t2884 - t2977;
t2786 = t2844 * t2951 + t2935;
t2918 = -t2786 * t2888 + t3007;
t2787 = t2812 - t2831;
t2917 = t2787 * t2888 - t3007;
t2916 = -t2791 * t2888 - t2995;
t2790 = -t2842 * t2951 + t2909;
t2915 = t2790 * t2888 + t2995;
t2814 = -g(3) * t2952 - t2936;
t2815 = -g(3) * t2953 + t2947;
t2913 = t2814 * t2888 + t2815 * t2884;
t2827 = t2880 * t2907 - t2853;
t2856 = t2880 * t2906;
t2828 = t2856 + t2933;
t2912 = t2827 * t2888 + t2828 * t2884;
t2878 = t2884 ^ 2;
t2838 = -t2878 * t2954 - t2934;
t2865 = t2888 * t2884 * t2954;
t2851 = t2865 - t2875;
t2911 = t2838 * t2888 + t2851 * t2884;
t2850 = t2865 + t2875;
t2879 = t2888 ^ 2;
t2854 = -t2879 * t2954 - t2934;
t2910 = t2850 * t2888 + t2854 * t2884;
t2864 = -qJDD(1) * t2885 - t2889 * t2890;
t2863 = qJDD(1) * t2889 - t2885 * t2890;
t2855 = (-t2878 - t2879) * t2954;
t2832 = -t2880 * t2898 - t2961;
t2829 = t2856 - t2933;
t2826 = t2938 * t2942 + t2853;
t2820 = -t2850 * t2884 + t2854 * t2888;
t2816 = -t2838 * t2884 + t2851 * t2888;
t2799 = -t2827 * t2884 + t2828 * t2888;
t2796 = -t2880 * t2829 + t2881 * t2910;
t2795 = t2881 * t2829 + t2880 * t2910;
t2784 = -t2880 * t2826 + t2881 * t2911;
t2783 = t2881 * t2826 + t2880 * t2911;
t2779 = -t2880 * t2855 + t2881 * t2912;
t2778 = t2881 * t2855 + t2880 * t2912;
t2770 = -t2814 * t2884 + t2815 * t2888;
t2760 = -t2880 * t2832 + t2881 * t2913;
t2759 = t2881 * t2832 + t2880 * t2913;
t2756 = t2789 * t2883 + t2971;
t2755 = t2788 * t2883 + t2971;
t2745 = -t2790 * t2884 + t2994;
t2744 = t2791 * t2884 - t2994;
t2743 = -t2787 * t2884 - t3006;
t2742 = t2786 * t2884 + t3006;
t2733 = t2758 * t2888 + t2980;
t2732 = t2757 * t2888 + t2980;
t2728 = t2881 * t2915 + t2993;
t2727 = t2881 * t2916 - t2993;
t2726 = t2880 * t2915 - t2992;
t2725 = t2880 * t2916 + t2992;
t2724 = t2881 * t2917 - t3005;
t2723 = t2881 * t2918 + t3005;
t2722 = t2880 * t2917 + t3002;
t2721 = t2880 * t2918 - t3002;
t2715 = t2746 * t2887 + t2747 * t2883;
t2703 = -t2880 * t2756 + t2881 * t2919;
t2702 = -t2880 * t2755 + t2881 * t2920;
t2701 = t2881 * t2756 + t2880 * t2919;
t2700 = t2881 * t2755 + t2880 * t2920;
t2693 = t2716 * t2888 + t2780 * t2884;
t2691 = t2729 * t2883 - t2730 * t2887;
t2681 = -t2880 * t2715 + t2881 * t2921;
t2680 = t2881 * t2715 + t2880 * t2921;
t2679 = t2692 * t2888 + t2731 * t2884;
t2664 = -t2880 * t2691 + t2881 * t2929;
t2663 = t2881 * t2691 + t2880 * t2929;
t2657 = -t2661 * t2887 + t2705 * t2883;
t2655 = -t2659 * t2887 + t2686 * t2883;
t2654 = t2658 * t2888 + t2662 * t2884;
t2653 = t2656 * t2888 + t2660 * t2884;
t2652 = -t2880 * t2657 + t2881 * t2930;
t2651 = t2881 * t2657 + t2880 * t2930;
t2650 = -t2880 * t2655 + t2881 * t2931;
t2649 = t2881 * t2655 + t2880 * t2931;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2864, -t2863, 0, -t2868 * t2885 + t2869 * t2889, 0, 0, 0, 0, 0, 0, -t2796 * t2885 + t2820 * t2889, -t2784 * t2885 + t2816 * t2889, -t2779 * t2885 + t2799 * t2889, -t2760 * t2885 + t2770 * t2889, 0, 0, 0, 0, 0, 0, -t2723 * t2885 + t2742 * t2889, -t2727 * t2885 + t2744 * t2889, -t2703 * t2885 + t2733 * t2889, -t2681 * t2885 + t2693 * t2889, 0, 0, 0, 0, 0, 0, -t2702 * t2885 + t2732 * t2889, -t2724 * t2885 + t2743 * t2889, -t2728 * t2885 + t2745 * t2889, -t2664 * t2885 + t2679 * t2889, 0, 0, 0, 0, 0, 0, t3022, t3029, t3023, -t2652 * t2885 + t2654 * t2889, 0, 0, 0, 0, 0, 0, t3022, t3023, -t3029, -t2650 * t2885 + t2653 * t2889; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2863, t2864, 0, t2868 * t2889 + t2869 * t2885, 0, 0, 0, 0, 0, 0, t2796 * t2889 + t2820 * t2885, t2784 * t2889 + t2816 * t2885, t2779 * t2889 + t2799 * t2885, t2760 * t2889 + t2770 * t2885, 0, 0, 0, 0, 0, 0, t2723 * t2889 + t2742 * t2885, t2727 * t2889 + t2744 * t2885, t2703 * t2889 + t2733 * t2885, t2681 * t2889 + t2693 * t2885, 0, 0, 0, 0, 0, 0, t2702 * t2889 + t2732 * t2885, t2724 * t2889 + t2743 * t2885, t2728 * t2889 + t2745 * t2885, t2664 * t2889 + t2679 * t2885, 0, 0, 0, 0, 0, 0, t3020, -t3028, t3021, t2652 * t2889 + t2654 * t2885, 0, 0, 0, 0, 0, 0, t3020, t3021, t3028, t2650 * t2889 + t2653 * t2885; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2795, t2783, t2778, t2759, 0, 0, 0, 0, 0, 0, t2721, t2725, t2701, t2680, 0, 0, 0, 0, 0, 0, t2700, t2722, t2726, t2663, 0, 0, 0, 0, 0, 0, t3013, -t2671, t3012, t2651, 0, 0, 0, 0, 0, 0, t3013, t3012, t2671, t2649; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2890, -qJDD(1), 0, t2869, 0, 0, 0, 0, 0, 0, t2820, t2816, t2799, t2770, 0, 0, 0, 0, 0, 0, t2742, t2744, t2733, t2693, 0, 0, 0, 0, 0, 0, t2732, t2743, t2745, t2679, 0, 0, 0, 0, 0, 0, t2996, -t2688, t2997, t2654, 0, 0, 0, 0, 0, 0, t2996, t2997, t2688, t2653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2890, 0, t2868, 0, 0, 0, 0, 0, 0, t2796, t2784, t2779, t2760, 0, 0, 0, 0, 0, 0, t2723, t2727, t2703, t2681, 0, 0, 0, 0, 0, 0, t2702, t2724, t2728, t2664, 0, 0, 0, 0, 0, 0, t3011, -t2673, t3010, t2652, 0, 0, 0, 0, 0, 0, t3011, t3010, t2673, t2650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2795, t2783, t2778, t2759, 0, 0, 0, 0, 0, 0, t2721, t2725, t2701, t2680, 0, 0, 0, 0, 0, 0, t2700, t2722, t2726, t2663, 0, 0, 0, 0, 0, 0, t3013, -t2671, t3012, t2651, 0, 0, 0, 0, 0, 0, t3013, t3012, t2671, t2649; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2854, t2851, t2828, t2815, 0, 0, 0, 0, 0, 0, t2766, -t2774, t2758, t2716, 0, 0, 0, 0, 0, 0, t2757, -t2766, t2774, t2692, 0, 0, 0, 0, 0, 0, t2987, -t2710, t2988, t2658, 0, 0, 0, 0, 0, 0, t2987, t2988, t2710, t2656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2850, t2838, t2827, t2814, 0, 0, 0, 0, 0, 0, -t2786, -t2791, -t2968, -t2780, 0, 0, 0, 0, 0, 0, -t2968, t2787, t2790, -t2731, 0, 0, 0, 0, 0, 0, -t2976, -t2740, -t2967, -t2662, 0, 0, 0, 0, 0, 0, -t2976, -t2967, t2740, -t2660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2829, t2826, t2855, t2832, 0, 0, 0, 0, 0, 0, -t2765, t2771, t2756, t2715, 0, 0, 0, 0, 0, 0, t2755, t2765, -t2771, t2691, 0, 0, 0, 0, 0, 0, t2986, t2712, t2989, t2657, 0, 0, 0, 0, 0, 0, t2986, t2989, -t2712, t2655; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2813, t2902, t2789, t2747, 0, 0, 0, 0, 0, 0, t2788, -t2813, -t2902, t2729, 0, 0, 0, 0, 0, 0, t2904, t2752, t2768, t2705, 0, 0, 0, 0, 0, 0, t2904, t2768, -t2752, t2686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2804, t2940, t2899, t2746, 0, 0, 0, 0, 0, 0, t2899, t2804, -t2940, -t2730, 0, 0, 0, 0, 0, 0, -t2975, t2734, -t2965, -t2661, 0, 0, 0, 0, 0, 0, -t2975, -t2965, -t2734, -t2659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2786, t2791, t2968, t2780, 0, 0, 0, 0, 0, 0, t2968, -t2787, -t2790, t2731, 0, 0, 0, 0, 0, 0, t2976, t2740, t2967, t2662, 0, 0, 0, 0, 0, 0, t2976, t2967, -t2740, t2660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2968, -t2787, -t2790, t2731, 0, 0, 0, 0, 0, 0, t2976, t2740, t2967, t2662, 0, 0, 0, 0, 0, 0, t2976, t2967, -t2740, t2660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2788, t2813, t2902, -t2729, 0, 0, 0, 0, 0, 0, -t2904, -t2752, -t2768, -t2705, 0, 0, 0, 0, 0, 0, -t2904, -t2768, t2752, -t2686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2899, -t2804, t2940, t2730, 0, 0, 0, 0, 0, 0, t2975, -t2734, t2965, t2661, 0, 0, 0, 0, 0, 0, t2975, t2965, t2734, t2659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2969, t2763, t2905, t2683, 0, 0, 0, 0, 0, 0, t2969, t2905, -t2763, t2669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2894, -t2944, t2900, t2926, 0, 0, 0, 0, 0, 0, -t2894, t2900, t2944, t2892; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2904, t2752, t2768, t2705, 0, 0, 0, 0, 0, 0, t2904, t2768, -t2752, t2686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2969, t2905, -t2763, t2669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2904, t2768, -t2752, t2686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2894, -t2900, -t2944, -t2892;];
f_new_reg  = t1;
