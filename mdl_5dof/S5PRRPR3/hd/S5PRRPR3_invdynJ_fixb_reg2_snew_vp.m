% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:43
% EndTime: 2019-12-05 16:19:50
% DurationCPUTime: 2.13s
% Computational Cost: add. (7151->276), mult. (16141->398), div. (0->0), fcn. (11362->10), ass. (0->169)
t202 = -2 * qJD(4);
t143 = sin(pkin(9));
t144 = cos(pkin(9));
t150 = cos(qJ(3));
t147 = sin(qJ(3));
t173 = qJD(2) * t147;
t113 = -t144 * t150 * qJD(2) + t143 * t173;
t115 = (t150 * t143 + t147 * t144) * qJD(2);
t99 = t115 * t113;
t196 = qJDD(3) - t99;
t201 = t143 * t196;
t200 = t144 * t196;
t146 = sin(qJ(5));
t137 = qJDD(3) + qJDD(5);
t149 = cos(qJ(5));
t88 = t149 * t113 + t115 * t146;
t90 = -t113 * t146 + t115 * t149;
t63 = t90 * t88;
t197 = -t63 + t137;
t199 = t146 * t197;
t198 = t149 * t197;
t132 = t147 * qJDD(2);
t166 = qJD(2) * qJD(3);
t164 = t150 * t166;
t120 = t132 + t164;
t133 = t150 * qJDD(2);
t165 = t147 * t166;
t158 = t133 - t165;
t101 = t144 * t120 + t143 * t158;
t109 = qJD(3) * t113;
t78 = t109 + t101;
t152 = qJD(2) ^ 2;
t129 = t147 * t152 * t150;
t124 = qJDD(3) + t129;
t141 = -g(3) + qJDD(1);
t179 = sin(pkin(8));
t180 = cos(pkin(8));
t122 = -t180 * g(1) - t179 * g(2);
t148 = sin(qJ(2));
t156 = t179 * g(1) - t180 * g(2);
t194 = cos(qJ(2));
t154 = -t194 * t122 - t148 * t156;
t178 = qJDD(2) * pkin(6);
t97 = -t152 * pkin(2) - t154 + t178;
t84 = -t150 * t141 + t147 * t97;
t69 = (-t120 + t164) * qJ(4) + t124 * pkin(3) - t84;
t123 = qJD(3) * pkin(3) - qJ(4) * t173;
t140 = t150 ^ 2;
t135 = t140 * t152;
t85 = t147 * t141 + t150 * t97;
t70 = -pkin(3) * t135 + t158 * qJ(4) - qJD(3) * t123 + t85;
t159 = t115 * t202 - t143 * t70 + t144 * t69;
t31 = t113 * t202 + t143 * t69 + t144 * t70;
t195 = -t78 * pkin(7) + t159;
t86 = t88 ^ 2;
t87 = t90 ^ 2;
t111 = t113 ^ 2;
t112 = t115 ^ 2;
t138 = qJD(3) + qJD(5);
t136 = t138 ^ 2;
t153 = pkin(4) * t196 + t195;
t100 = -t120 * t143 + t144 * t158;
t102 = qJD(3) * pkin(4) - pkin(7) * t115;
t28 = -pkin(4) * t111 + pkin(7) * t100 - qJD(3) * t102 + t31;
t11 = t146 * t28 - t149 * t153;
t183 = t149 * t28;
t12 = t146 * t153 + t183;
t6 = -t11 * t149 + t12 * t146;
t193 = t143 * t6;
t192 = t144 * t6;
t160 = -t122 * t148 + t194 * t156;
t96 = -qJDD(2) * pkin(2) - t152 * pkin(6) - t160;
t71 = -t158 * pkin(3) - qJ(4) * t135 + t123 * t173 + qJDD(4) + t96;
t190 = t143 * t71;
t93 = qJDD(3) + t99;
t189 = t143 * t93;
t188 = t144 * t71;
t187 = t144 * t93;
t38 = -t100 * pkin(4) - t111 * pkin(7) + t102 * t115 + t71;
t186 = t146 * t38;
t60 = t63 + t137;
t185 = t146 * t60;
t15 = t143 * t31 + t144 * t159;
t184 = t147 * t15;
t182 = t149 * t38;
t181 = t149 * t60;
t177 = t138 * t146;
t176 = t138 * t149;
t175 = t147 * t124;
t125 = qJDD(3) - t129;
t174 = t150 * t125;
t172 = qJD(3) * t115;
t171 = qJD(3) * t143;
t170 = qJD(3) * t144;
t167 = qJD(5) + t138;
t7 = t11 * t146 + t149 * t12;
t16 = -t143 * t159 + t144 * t31;
t162 = t147 * t84 + t150 * t85;
t161 = -t149 * t100 + t101 * t146;
t121 = t133 - 0.2e1 * t165;
t157 = t100 * t146 + t101 * t149;
t76 = t100 + t172;
t155 = (-qJD(5) + t138) * t90 - t161;
t53 = -qJD(5) * t88 + t157;
t151 = qJD(3) ^ 2;
t139 = t147 ^ 2;
t134 = t139 * t152;
t128 = -t135 - t151;
t127 = -t134 - t151;
t119 = t132 + 0.2e1 * t164;
t105 = -t112 - t151;
t104 = -t112 + t151;
t103 = t111 - t151;
t91 = -t151 - t111;
t83 = t138 * t88;
t82 = -t87 + t136;
t81 = t86 - t136;
t79 = -t87 - t136;
t77 = -t109 + t101;
t75 = -t100 + t172;
t74 = -t111 - t112;
t73 = -t105 * t143 - t187;
t72 = t105 * t144 - t189;
t67 = t144 * t91 - t201;
t66 = t143 * t91 + t200;
t62 = t87 - t86;
t58 = -t136 - t86;
t57 = (t146 * t90 - t149 * t88) * t138;
t56 = (-t146 * t88 - t149 * t90) * t138;
t55 = t143 * t78 + t144 * t76;
t54 = t143 * t76 - t144 * t78;
t52 = -qJD(5) * t90 - t161;
t51 = -t86 - t87;
t50 = t149 * t81 - t185;
t49 = -t146 * t82 + t198;
t48 = t146 * t81 + t181;
t47 = t149 * t82 + t199;
t46 = -t146 * t79 - t181;
t45 = t149 * t79 - t185;
t44 = t53 + t83;
t43 = t53 - t83;
t42 = -t167 * t88 + t157;
t39 = t167 * t90 + t161;
t37 = t149 * t53 - t90 * t177;
t36 = t146 * t53 + t90 * t176;
t35 = -t146 * t52 + t88 * t176;
t34 = t149 * t52 + t88 * t177;
t33 = t149 * t58 - t199;
t32 = t146 * t58 + t198;
t26 = -t143 * t45 + t144 * t46;
t25 = t143 * t46 + t144 * t45;
t24 = -pkin(7) * t45 + t182;
t23 = t146 * t44 + t149 * t155;
t22 = -t146 * t43 - t149 * t39;
t21 = t146 * t155 - t149 * t44;
t20 = -t146 * t39 + t149 * t43;
t19 = -pkin(7) * t32 + t186;
t18 = -t143 * t32 + t144 * t33;
t17 = t143 * t33 + t144 * t32;
t14 = -pkin(4) * t42 + pkin(7) * t46 + t186;
t13 = -pkin(4) * t39 + pkin(7) * t33 - t182;
t9 = -t143 * t21 + t144 * t23;
t8 = t143 * t23 + t144 * t21;
t5 = -pkin(4) * t38 + pkin(7) * t7;
t4 = -pkin(7) * t21 - t6;
t3 = -pkin(4) * t51 + pkin(7) * t23 + t7;
t2 = t144 * t7 - t193;
t1 = t143 * t7 + t192;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, t124 * t150 + t128 * t147, -t125 * t147 + t127 * t150, 0, t147 * t85 - t150 * t84, 0, 0, 0, 0, 0, 0, t147 * t67 + t150 * t66, t147 * t73 + t150 * t72, t147 * t55 + t150 * t54, t147 * t16 + t15 * t150, 0, 0, 0, 0, 0, 0, t147 * t18 + t150 * t17, t147 * t26 + t150 * t25, t147 * t9 + t150 * t8, t1 * t150 + t147 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t160, t154, 0, 0, (t120 + t164) * t147, t119 * t150 + t121 * t147, t175 + t150 * (-t134 + t151), t121 * t150, t147 * (t135 - t151) + t174, 0, -t150 * t96 + pkin(2) * t121 + pkin(6) * (t128 * t150 - t175), t147 * t96 - pkin(2) * t119 + pkin(6) * (-t127 * t147 - t174), pkin(2) * (t134 + t135) + (t139 + t140) * t178 + t162, -pkin(2) * t96 + pkin(6) * t162, t147 * (t101 * t144 - t115 * t171) + t150 * (t101 * t143 + t115 * t170), t147 * (-t143 * t77 - t144 * t75) + t150 * (-t143 * t75 + t144 * t77), t147 * (-t104 * t143 + t200) + t150 * (t104 * t144 + t201), t147 * (-t100 * t143 + t113 * t170) + t150 * (t100 * t144 + t113 * t171), t147 * (t103 * t144 - t189) + t150 * (t103 * t143 + t187), (t147 * (-t113 * t144 + t115 * t143) + t150 * (-t113 * t143 - t115 * t144)) * qJD(3), t147 * (-qJ(4) * t66 + t190) + t150 * (-pkin(3) * t75 + qJ(4) * t67 - t188) - pkin(2) * t75 + pkin(6) * (-t147 * t66 + t150 * t67), t147 * (-qJ(4) * t72 + t188) + t150 * (-pkin(3) * t77 + qJ(4) * t73 + t190) - pkin(2) * t77 + pkin(6) * (-t147 * t72 + t150 * t73), t147 * (-qJ(4) * t54 - t15) + t150 * (-pkin(3) * t74 + qJ(4) * t55 + t16) - pkin(2) * t74 + pkin(6) * (-t147 * t54 + t150 * t55), -qJ(4) * t184 + t150 * (-pkin(3) * t71 + qJ(4) * t16) - pkin(2) * t71 + pkin(6) * (t150 * t16 - t184), t147 * (-t143 * t36 + t144 * t37) + t150 * (t143 * t37 + t144 * t36), t147 * (-t143 * t20 + t144 * t22) + t150 * (t143 * t22 + t144 * t20), t147 * (-t143 * t47 + t144 * t49) + t150 * (t143 * t49 + t144 * t47), t147 * (-t143 * t34 + t144 * t35) + t150 * (t143 * t35 + t144 * t34), t147 * (-t143 * t48 + t144 * t50) + t150 * (t143 * t50 + t144 * t48), t147 * (-t143 * t56 + t144 * t57) + t150 * (t143 * t57 + t144 * t56), t147 * (-qJ(4) * t17 - t13 * t143 + t144 * t19) + t150 * (-pkin(3) * t39 + qJ(4) * t18 + t13 * t144 + t143 * t19) - pkin(2) * t39 + pkin(6) * (-t147 * t17 + t150 * t18), t147 * (-qJ(4) * t25 - t14 * t143 + t144 * t24) + t150 * (-pkin(3) * t42 + qJ(4) * t26 + t14 * t144 + t143 * t24) - pkin(2) * t42 + pkin(6) * (-t147 * t25 + t150 * t26), t147 * (-qJ(4) * t8 - t143 * t3 + t144 * t4) + t150 * (-pkin(3) * t51 + qJ(4) * t9 + t143 * t4 + t144 * t3) - pkin(2) * t51 + pkin(6) * (-t147 * t8 + t150 * t9), t147 * (-pkin(7) * t192 - qJ(4) * t1 - t143 * t5) + t150 * (-pkin(3) * t38 - pkin(7) * t193 + qJ(4) * t2 + t144 * t5) - pkin(2) * t38 + pkin(6) * (-t1 * t147 + t150 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t134 - t135, t132, t129, t133, qJDD(3), -t84, -t85, 0, 0, t99, t112 - t111, t78, -t99, t76, qJDD(3), pkin(3) * t66 + t159, pkin(3) * t72 - t31, pkin(3) * t54, pkin(3) * t15, t63, t62, t44, -t63, t155, t137, pkin(3) * t17 + pkin(4) * t32 - t11, -t183 - t146 * t195 + pkin(3) * t25 + (-t146 * t196 + t45) * pkin(4), pkin(3) * t8 + pkin(4) * t21, pkin(3) * t1 + pkin(4) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t77, t74, t71, 0, 0, 0, 0, 0, 0, t39, t42, t51, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t62, t44, -t63, t155, t137, -t11, -t12, 0, 0;];
tauJ_reg = t10;
