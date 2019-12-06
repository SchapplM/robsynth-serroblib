% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:34
% EndTime: 2019-12-05 16:37:42
% DurationCPUTime: 1.93s
% Computational Cost: add. (1413->269), mult. (3864->426), div. (0->0), fcn. (2861->10), ass. (0->138)
t101 = sin(pkin(10));
t107 = sin(qJ(2));
t109 = cos(qJ(3));
t102 = sin(pkin(5));
t159 = qJD(1) * t102;
t103 = cos(pkin(10));
t110 = cos(qJ(2));
t165 = t103 * t110;
t106 = sin(qJ(3));
t123 = pkin(3) * t106 - qJ(4) * t109;
t62 = t123 * qJD(3) - t106 * qJD(4);
t182 = t101 * t62 - (t101 * t107 + t109 * t165) * t159;
t142 = t110 * t159;
t143 = t107 * t159;
t170 = t101 * t109;
t181 = t142 * t170 + (-t143 + t62) * t103;
t104 = cos(pkin(5));
t158 = qJD(1) * t104;
t80 = qJD(2) * pkin(7) + t143;
t180 = -t106 * t80 + t109 * t158;
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t91 = t106 * t158;
t54 = t109 * t80 + t91;
t44 = qJD(3) * qJ(4) + t54;
t85 = -t109 * pkin(3) - t106 * qJ(4) - pkin(2);
t55 = t85 * qJD(2) - t142;
t15 = t101 * t55 + t103 * t44;
t155 = qJD(2) * t109;
t10 = -pkin(8) * t155 + t15;
t135 = -qJD(3) * pkin(3) + qJD(4);
t41 = t135 - t180;
t156 = qJD(2) * t106;
t97 = t103 * qJD(3);
t73 = t101 * t156 - t97;
t75 = t101 * qJD(3) + t103 * t156;
t12 = t73 * pkin(4) - t75 * pkin(8) + t41;
t122 = t105 * t10 - t108 * t12;
t126 = pkin(4) * t101 - pkin(8) * t103;
t116 = t126 * t155;
t129 = qJD(2) * t142;
t152 = qJD(3) * t109;
t27 = qJD(3) * t91 + t106 * t129 + t80 * t152;
t13 = qJD(3) * t116 + t27;
t149 = qJD(2) * qJD(3);
t139 = t106 * t149;
t25 = t109 * t129 + (qJD(4) + t180) * qJD(3);
t33 = (t62 + t143) * qJD(2);
t8 = t101 * t33 + t103 * t25;
t6 = pkin(8) * t139 + t8;
t1 = -t122 * qJD(5) + t105 * t13 + t108 * t6;
t146 = pkin(7) * t101 + pkin(4);
t154 = qJD(3) * t106;
t179 = -t146 * t154 - t181;
t148 = pkin(7) * t154;
t134 = t101 * t148;
t178 = t134 + t181;
t133 = t103 * t148;
t177 = -t133 + t182;
t77 = t123 * qJD(2);
t24 = t101 * t77 + t103 * t180;
t166 = t103 * t109;
t57 = pkin(7) * t166 + t101 * t85;
t176 = qJD(2) * pkin(2);
t174 = t103 * t85;
t173 = t27 * t101;
t172 = t27 * t103;
t171 = -t106 ^ 2 + t109 ^ 2;
t169 = t102 * t107;
t112 = qJD(2) ^ 2;
t168 = t102 * t112;
t167 = t103 * t108;
t164 = t105 * t106;
t163 = t108 * t109;
t111 = qJD(3) ^ 2;
t162 = t111 * t106;
t161 = t111 * t109;
t160 = -qJD(5) - t73;
t157 = qJD(2) * t102;
t153 = qJD(3) * t108;
t151 = qJD(5) * t105;
t150 = qJD(5) * t108;
t147 = t107 * t168;
t145 = t107 * t157;
t144 = t110 * t157;
t141 = t105 * t155;
t138 = t109 * t149;
t137 = t160 ^ 2;
t82 = -t103 * pkin(4) - t101 * pkin(8) - pkin(3);
t136 = pkin(8) * t156 - qJD(5) * t82 + t24;
t132 = t106 * t142;
t131 = t106 * t144;
t130 = t109 * t144;
t128 = t103 * t138;
t86 = t101 * t138;
t119 = pkin(7) + t126;
t59 = t119 * t106;
t127 = -qJD(5) * t59 - (-pkin(7) * t103 + pkin(8)) * t154 - t182;
t81 = -t142 - t176;
t125 = -t81 - t142;
t4 = t108 * t10 + t105 * t12;
t7 = -t101 * t25 + t103 * t33;
t14 = -t101 * t44 + t103 * t55;
t23 = -t101 * t180 + t103 * t77;
t65 = t104 * t106 + t109 * t169;
t38 = -t102 * t110 * t101 + t65 * t103;
t64 = -t104 * t109 + t106 * t169;
t121 = t64 * t105 + t38 * t108;
t120 = -t38 * t105 + t64 * t108;
t47 = t105 * t75 + t108 * t155;
t117 = t103 * t163 + t164;
t67 = t103 * t164 + t163;
t46 = -t109 * pkin(8) + t57;
t115 = qJD(5) * t46 - t119 * t152 + t132;
t114 = qJD(3) * (-t125 - t176);
t18 = -t47 * qJD(5) + t105 * t139 + t108 * t128;
t113 = -qJ(4) * t154 + (t135 - t41) * t109;
t2 = -t4 * qJD(5) - t105 * t6 + t108 * t13;
t19 = -qJD(5) * t141 + t105 * t128 - t108 * t139 + t75 * t150;
t98 = t101 ^ 2;
t68 = -t105 * t109 + t106 * t167;
t61 = t117 * qJD(2);
t60 = t103 * t141 - t108 * t156;
t56 = -pkin(7) * t170 + t174;
t49 = t108 * t75 - t141;
t45 = t146 * t109 - t174;
t40 = t65 * qJD(3) + t131;
t39 = -t64 * qJD(3) + t130;
t37 = t65 * t101 + t102 * t165;
t29 = -t106 * t153 - t109 * t151 + (t105 * t152 + t106 * t150) * t103;
t28 = t117 * qJD(3) - t67 * qJD(5);
t26 = t116 + t54;
t22 = t101 * t145 + t39 * t103;
t21 = t39 * t101 - t103 * t145;
t16 = -pkin(4) * t156 - t23;
t9 = pkin(4) * t155 - t14;
t5 = -pkin(4) * t139 - t7;
t3 = [0, 0, -t147, -t110 * t168, 0, 0, 0, 0, 0, -t109 * t147 + (-t40 - t131) * qJD(3), t106 * t147 + (-t39 - t130) * qJD(3), t40 * t73 + (t109 * t21 + (-t106 * t37 + t64 * t170) * qJD(3)) * qJD(2), t40 * t75 + (t109 * t22 + (-t106 * t38 + t64 * t166) * qJD(3)) * qJD(2), t21 * t75 - t22 * t73 + (-t38 * t101 + t103 * t37) * t138, -t14 * t21 + t15 * t22 + t27 * t64 - t7 * t37 + t8 * t38 + t41 * t40, 0, 0, 0, 0, 0, -(-t121 * qJD(5) - t22 * t105 + t40 * t108) * t160 + t120 * t86 + t21 * t47 + t37 * t19, (t120 * qJD(5) + t40 * t105 + t22 * t108) * t160 - t121 * t86 + t21 * t49 + t37 * t18; 0, 0, 0, 0, 0.2e1 * t106 * t138, 0.2e1 * t171 * t149, t161, -t162, 0, -pkin(7) * t161 + t106 * t114, pkin(7) * t162 + t109 * t114, (-t73 * t142 + t173 + (qJD(2) * t56 + t14) * qJD(3)) * t106 + (-t7 + (pkin(7) * t73 + t101 * t41) * qJD(3) + (t134 - t178) * qJD(2)) * t109, (-t75 * t142 + t172 + (-qJD(2) * t57 - t15) * qJD(3)) * t106 + (t8 + (pkin(7) * t75 + t103 * t41) * qJD(3) + (t133 + t177) * qJD(2)) * t109, -t178 * t75 - t177 * t73 + (-t101 * t8 - t103 * t7) * t106 + (-t101 * t15 - t103 * t14 + (-t101 * t57 - t103 * t56) * qJD(2)) * t152, -t41 * t132 + t7 * t56 + t8 * t57 + t177 * t15 + t178 * t14 + (t106 * t27 + t41 * t152) * pkin(7), t18 * t68 + t49 * t28, -t18 * t67 - t68 * t19 - t28 * t47 - t49 * t29, -t28 * t160 + (t106 * t18 + (qJD(2) * t68 + t49) * t152) * t101, t29 * t160 + (-t106 * t19 + (-qJD(2) * t67 - t47) * t152) * t101, (-t101 * t160 + t98 * t156) * t152, t45 * t19 + t9 * t29 + t5 * t67 - (t127 * t105 - t115 * t108) * t160 + t179 * t47 + (t2 * t106 + ((-t105 * t46 + t108 * t59) * qJD(2) - t122) * t152) * t101, t45 * t18 + t9 * t28 + t5 * t68 - (t115 * t105 + t127 * t108) * t160 + t179 * t49 + (-t1 * t106 + (-(t105 * t59 + t108 * t46) * qJD(2) - t4) * t152) * t101; 0, 0, 0, 0, -t106 * t112 * t109, -t171 * t112, 0, 0, 0, t54 * qJD(3) - t81 * t156 - t27, t125 * t155, -t172 - t54 * t73 + (t113 * t101 - t106 * t14 + t109 * t23) * qJD(2), t173 - t54 * t75 + (t113 * t103 + t106 * t15 - t109 * t24) * qJD(2), t23 * t75 + t24 * t73 + (-qJD(4) * t73 + t14 * t155 + t8) * t103 + (qJD(4) * t75 + t15 * t155 - t7) * t101, -t27 * pkin(3) - t14 * t23 - t15 * t24 - t41 * t54 + (-t101 * t14 + t103 * t15) * qJD(4) + (-t7 * t101 + t8 * t103) * qJ(4), -t49 * t61 + (t108 * t18 - t49 * t151) * t101, t61 * t47 + t49 * t60 + (-t105 * t18 - t108 * t19 + (t105 * t47 - t108 * t49) * qJD(5)) * t101, -t18 * t103 - (-t101 * t151 - t61) * t160 + (-t101 * t49 + t98 * t153) * t155, t19 * t103 - (-t101 * t150 + t60) * t160 + (-qJD(3) * t105 * t98 + t101 * t47) * t155, (t160 - t97) * t101 * t155, -t16 * t47 - t9 * t60 - (t136 * t105 - t108 * t26) * t160 + (-(-qJ(4) * t150 - qJD(4) * t105) * t160 - t2) * t103 + (t9 * t150 + qJ(4) * t19 + qJD(4) * t47 + t5 * t105 + ((-t105 * t103 * qJ(4) + t108 * t82) * qJD(3) + t122) * t155) * t101, -t16 * t49 - t9 * t61 - (t105 * t26 + t136 * t108) * t160 + ((-qJ(4) * t151 + qJD(4) * t108) * t160 + t1) * t103 + (-t9 * t151 + qJ(4) * t18 + qJD(4) * t49 + t5 * t108 + (-(qJ(4) * t167 + t105 * t82) * qJD(3) + t4) * t155) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t155 + t86, (t73 + t97) * t155, -t73 ^ 2 - t75 ^ 2, t14 * t75 + t15 * t73 + t27, 0, 0, 0, 0, 0, -t105 * t137 + t108 * t86 - t75 * t47, -t105 * t86 - t108 * t137 - t75 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t47, -t47 ^ 2 + t49 ^ 2, -t160 * t47 + t18, -t160 * t49 - t19, t86, -t160 * t4 - t9 * t49 + t2, t122 * t160 + t9 * t47 - t1;];
tauc_reg = t3;
