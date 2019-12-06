% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:27
% EndTime: 2019-12-05 16:56:34
% DurationCPUTime: 1.65s
% Computational Cost: add. (3977->231), mult. (7799->307), div. (0->0), fcn. (5460->10), ass. (0->160)
t141 = sin(qJ(4));
t144 = cos(qJ(4));
t142 = sin(qJ(3));
t172 = qJD(2) * t142;
t107 = -t144 * qJD(3) + t141 * t172;
t145 = cos(qJ(3));
t168 = qJD(2) * qJD(3);
t163 = t145 * t168;
t167 = t142 * qJDD(2);
t112 = t163 + t167;
t83 = -t107 * qJD(4) + t141 * qJDD(3) + t144 * t112;
t125 = t145 * qJD(2) - qJD(4);
t97 = t107 * t125;
t68 = t83 - t97;
t212 = qJ(5) * t68;
t128 = t142 * t168;
t166 = t145 * qJDD(2);
t113 = -t128 + t166;
t106 = -qJDD(4) + t113;
t109 = t141 * qJD(3) + t144 * t172;
t88 = t109 * t107;
t204 = -t106 - t88;
t211 = pkin(4) * t204;
t104 = t107 ^ 2;
t147 = qJD(2) ^ 2;
t194 = t145 * pkin(3);
t158 = -t142 * pkin(8) - t194;
t180 = sin(pkin(9));
t181 = cos(pkin(9));
t117 = -t181 * g(1) - t180 * g(2);
t143 = sin(qJ(2));
t146 = cos(qJ(2));
t137 = sin(pkin(5));
t138 = cos(pkin(5));
t152 = t180 * g(1) - t181 * g(2);
t151 = t138 * t152;
t173 = -g(3) + qJDD(1);
t208 = t137 * t173 + t151;
t78 = t146 * t117 + t208 * t143;
t71 = -t147 * pkin(2) + qJDD(2) * pkin(7) + t78;
t161 = t147 * t158 + t71;
t202 = qJD(3) ^ 2;
t149 = -t137 * t152 + t138 * t173;
t91 = t145 * t149;
t38 = -qJDD(3) * pkin(3) - t202 * pkin(8) + t161 * t142 - t91;
t160 = -t144 * qJDD(3) + t141 * t112;
t82 = -t109 * qJD(4) - t160;
t92 = -t125 * pkin(4) - t109 * qJ(5);
t20 = -t82 * pkin(4) - t104 * qJ(5) + t109 * t92 + qJDD(5) + t38;
t210 = t141 * t204;
t209 = t144 * t204;
t148 = t142 * t149;
t39 = -t202 * pkin(3) + qJDD(3) * pkin(8) + t161 * t145 + t148;
t156 = -t113 + t128;
t157 = t112 + t163;
t159 = t143 * t117 - t208 * t146;
t70 = -qJDD(2) * pkin(2) - t147 * pkin(7) + t159;
t44 = t156 * pkin(3) - t157 * pkin(8) + t70;
t19 = t141 * t44 + t144 * t39;
t154 = t82 * qJ(5) - 0.2e1 * qJD(5) * t107 + t125 * t92 + t19;
t105 = t109 ^ 2;
t123 = t125 ^ 2;
t85 = -t105 - t123;
t207 = -t154 + (t104 + t85) * pkin(4);
t205 = t83 + t97;
t64 = (qJD(4) + t125) * t109 + t160;
t84 = -t123 - t104;
t47 = t141 * t84 + t209;
t201 = pkin(3) * t47;
t74 = t106 - t88;
t185 = t141 * t74;
t52 = t144 * t85 + t185;
t200 = pkin(3) * t52;
t170 = qJD(5) * t109;
t100 = -0.2e1 * t170;
t18 = t141 * t39 - t144 * t44;
t153 = -t18 + t211 - t212;
t13 = t100 + t153;
t199 = pkin(4) * t13;
t198 = pkin(4) * t68;
t34 = -t141 * t64 - t144 * t68;
t197 = pkin(8) * t34;
t196 = pkin(8) * t47;
t195 = pkin(8) * t52;
t35 = t141 * t68 - t144 * t64;
t73 = -t104 - t105;
t22 = t142 * t73 + t145 * t35;
t193 = -pkin(2) * t34 + pkin(7) * t22;
t48 = t144 * t84 - t210;
t63 = (qJD(4) - t125) * t109 + t160;
t26 = t142 * t63 + t145 * t48;
t192 = -pkin(2) * t47 + pkin(7) * t26;
t183 = t144 * t74;
t53 = -t141 * t85 + t183;
t28 = t142 * t205 + t145 * t53;
t191 = -pkin(2) * t52 + pkin(7) * t28;
t190 = -pkin(3) * t73 + pkin(8) * t35;
t189 = -pkin(3) * t63 + pkin(8) * t48;
t188 = -pkin(3) * t205 + pkin(8) * t53;
t186 = t141 * t38;
t184 = t144 * t38;
t179 = qJ(5) * t141;
t178 = qJ(5) * t144;
t177 = t125 * t141;
t176 = t125 * t144;
t124 = t142 * t147 * t145;
t118 = qJDD(3) + t124;
t175 = t142 * t118;
t119 = qJDD(3) - t124;
t174 = t145 * t119;
t164 = t145 * t88;
t10 = t141 * t18 + t144 * t19;
t54 = t142 * t71 - t91;
t55 = t145 * t71 + t148;
t23 = t142 * t54 + t145 * t55;
t9 = t141 * t19 - t144 * t18;
t150 = t153 + t211;
t134 = t145 ^ 2;
t133 = t142 ^ 2;
t132 = t134 * t147;
t130 = t133 * t147;
t122 = -t132 - t202;
t121 = -t130 - t202;
t116 = t130 + t132;
t115 = (t133 + t134) * qJDD(2);
t114 = -0.2e1 * t128 + t166;
t111 = 0.2e1 * t163 + t167;
t101 = 0.2e1 * t170;
t94 = -t105 + t123;
t93 = t104 - t123;
t90 = -t142 * t121 - t174;
t89 = t145 * t122 - t175;
t86 = t105 - t104;
t69 = (t107 * t141 + t109 * t144) * t125;
t60 = -t109 * t176 + t141 * t83;
t59 = -t107 * t177 + t144 * t82;
t58 = t145 * t106 + t142 * (t107 * t144 - t109 * t141) * t125;
t57 = t141 * t93 - t183;
t56 = t144 * t94 + t210;
t41 = t142 * (t109 * t177 + t144 * t83) - t164;
t40 = t142 * (-t107 * t176 - t141 * t82) + t164;
t37 = -pkin(4) * t205 + qJ(5) * t74;
t33 = -t141 * t63 + t144 * t205;
t30 = t142 * (t144 * t93 + t185) + t145 * t64;
t29 = t142 * (-t141 * t94 + t209) - t145 * t68;
t24 = t142 * (-t141 * t205 - t144 * t63) - t145 * t86;
t16 = -qJ(5) * t85 + t20;
t15 = -pkin(4) * t63 + qJ(5) * t84 - t20;
t14 = -t104 * pkin(4) + t154;
t12 = t101 - t153 + t212;
t11 = t138 * (t142 * t53 - t145 * t205) + (t143 * t28 - t146 * t52) * t137;
t8 = t138 * (t142 * t48 - t145 * t63) + (t143 * t26 - t146 * t47) * t137;
t7 = -qJ(5) * t64 + (-t104 - t73) * pkin(4) + t154;
t6 = t138 * (t142 * t35 - t145 * t73) + (t143 * t22 - t146 * t34) * t137;
t5 = -pkin(4) * t20 + qJ(5) * t14;
t4 = t145 * t10 + t142 * t38;
t3 = -t141 * t13 + t144 * t14;
t2 = t144 * t13 + t141 * t14;
t1 = t142 * t20 + t145 * t3;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t173, 0, 0, 0, 0, 0, 0, (qJDD(2) * t146 - t143 * t147) * t137, (-qJDD(2) * t143 - t146 * t147) * t137, 0, t138 ^ 2 * t173 + (t143 * t78 - t146 * t159 - t151) * t137, 0, 0, 0, 0, 0, 0, t138 * (t145 * t118 + t142 * t122) + (t146 * t114 + t143 * t89) * t137, t138 * (-t142 * t119 + t145 * t121) + (-t146 * t111 + t143 * t90) * t137, (t115 * t143 + t116 * t146) * t137, t138 * (t142 * t55 - t145 * t54) + (t143 * t23 - t146 * t70) * t137, 0, 0, 0, 0, 0, 0, t8, t11, t6, t138 * (t142 * t10 - t145 * t38) + (t143 * t4 - t146 * t9) * t137, 0, 0, 0, 0, 0, 0, t8, t11, t6, t138 * (t142 * t3 - t145 * t20) + (t143 * t1 - t146 * t2) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t159, -t78, 0, 0, t157 * t142, t145 * t111 + t142 * t114, t175 + t145 * (-t130 + t202), -t156 * t145, t142 * (t132 - t202) + t174, 0, pkin(2) * t114 + pkin(7) * t89 - t145 * t70, -pkin(2) * t111 + pkin(7) * t90 + t142 * t70, pkin(2) * t116 + pkin(7) * t115 + t23, -pkin(2) * t70 + pkin(7) * t23, t41, t24, t29, t40, t30, t58, t142 * (t186 - t196) + t145 * (t18 - t201) + t192, t142 * (t184 - t195) + t145 * (t19 - t200) + t191, t142 * (-t9 - t197) - t34 * t194 + t193, pkin(7) * t4 + (-pkin(2) + t158) * t9, t41, t24, t29, t40, t30, t58, t142 * (-t141 * t15 - t178 * t204 - t196) + t145 * (t101 - t150 - t201) + t192, t142 * (-t141 * t37 + t144 * t16 - t195) + t145 * (-t200 - t207) + t191, t142 * (t144 * t12 - t141 * t7 - t197) + t145 * (-pkin(3) * t34 + t198) + t193, t142 * (-pkin(8) * t2 - t13 * t178 - t141 * t5) + t145 * (-pkin(3) * t2 - t199) - pkin(2) * t2 + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t130 - t132, t167, t124, t166, qJDD(3), -t54, -t55, 0, 0, t60, t33, t56, t59, t57, t69, -t184 + t189, t186 + t188, t10 + t190, -pkin(3) * t38 + pkin(8) * t10, t60, t33, t56, t59, t57, t69, t144 * t15 - t179 * t204 + t189, t141 * t16 + t144 * t37 + t188, t141 * t12 + t144 * t7 + t190, -pkin(3) * t20 + pkin(8) * t3 - t13 * t179 + t144 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t86, t68, -t88, -t64, -t106, -t18, -t19, 0, 0, t88, t86, t68, -t88, -t64, -t106, t100 + t150, t207, -t198, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t205, t73, t20;];
tauJ_reg = t17;
