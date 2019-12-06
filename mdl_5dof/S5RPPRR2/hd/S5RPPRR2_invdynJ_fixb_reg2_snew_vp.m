% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:57
% EndTime: 2019-12-05 17:40:06
% DurationCPUTime: 2.22s
% Computational Cost: add. (6540->253), mult. (15068->357), div. (0->0), fcn. (10658->8), ass. (0->166)
t141 = sin(qJ(4));
t138 = sin(pkin(8));
t139 = cos(pkin(8));
t144 = cos(qJ(4));
t156 = t138 * t144 + t139 * t141;
t119 = t156 * qJD(1);
t171 = t138 * t141;
t121 = (t139 * t144 - t171) * qJD(1);
t174 = t121 * t119;
t198 = qJDD(4) - t174;
t200 = t141 * t198;
t199 = t144 * t198;
t140 = sin(qJ(5));
t133 = qJDD(4) + qJDD(5);
t143 = cos(qJ(5));
t89 = t119 * t143 + t121 * t140;
t91 = -t119 * t140 + t121 * t143;
t70 = t91 * t89;
t194 = -t70 + t133;
t197 = t140 * t194;
t196 = t143 * t194;
t147 = qJD(1) ^ 2;
t142 = sin(qJ(1));
t145 = cos(qJ(1));
t163 = t142 * g(1) - g(2) * t145;
t158 = qJDD(2) - t163;
t154 = -t147 * qJ(2) + t158;
t186 = -qJ(3) - pkin(1);
t160 = -0.2e1 * qJD(1) * qJD(3) + qJDD(1) * t186 + t154;
t168 = t121 * qJD(4);
t78 = t156 * qJDD(1);
t103 = -t78 - t168;
t165 = t139 * qJDD(1);
t166 = t138 * qJDD(1);
t118 = -t141 * t166 + t144 * t165;
t169 = qJD(4) * t119;
t105 = t118 - t169;
t58 = -qJD(5) * t89 + t103 * t140 + t105 * t143;
t136 = qJD(4) + qJD(5);
t84 = t136 * t89;
t195 = -t84 + t58;
t134 = t138 ^ 2;
t135 = t139 ^ 2;
t170 = t134 + t135;
t193 = t170 * t147;
t190 = pkin(3) * t147;
t151 = (-pkin(6) * qJDD(1) - t138 * t190 + t160) * t139;
t189 = t138 * g(3);
t97 = -g(3) * t139 + t138 * t160;
t81 = -pkin(6) * t166 - t134 * t190 + t97;
t60 = t141 * t81 - t144 * (t151 + t189);
t191 = -t60 + (-t105 - t169) * pkin(7);
t87 = t89 ^ 2;
t88 = t91 ^ 2;
t116 = t119 ^ 2;
t117 = t121 ^ 2;
t132 = t136 ^ 2;
t150 = pkin(4) * t198 + t191;
t157 = qJD(4) * pkin(4) - pkin(7) * t121;
t61 = g(3) * t171 + t141 * t151 + t144 * t81;
t33 = -t116 * pkin(4) + t103 * pkin(7) - qJD(4) * t157 + t61;
t18 = t140 * t33 - t143 * t150;
t181 = t143 * t33;
t19 = t140 * t150 + t181;
t8 = t140 * t19 - t143 * t18;
t188 = t141 * t8;
t187 = t144 * t8;
t30 = t141 * t61 - t144 * t60;
t185 = t139 * t30;
t137 = qJDD(1) * qJ(2);
t159 = t145 * g(1) + t142 * g(2);
t155 = -t137 + t159;
t167 = qJD(2) * qJD(1);
t152 = -qJDD(3) + t155 - 0.2e1 * t167;
t85 = -pkin(3) * t166 + (pkin(6) * t170 - t186) * t147 + t152;
t54 = t103 * pkin(4) + t116 * pkin(7) - t121 * t157 + t85;
t184 = t140 * t54;
t66 = t70 + t133;
t183 = t140 * t66;
t182 = t141 * t85;
t180 = t143 * t54;
t179 = t143 * t66;
t178 = t144 * t85;
t177 = qJDD(1) * pkin(1);
t100 = qJDD(4) + t174;
t176 = t100 * t141;
t175 = t100 * t144;
t173 = t136 * t140;
t172 = t136 * t143;
t9 = t140 * t18 + t143 * t19;
t31 = t141 * t60 + t144 * t61;
t162 = -t103 * t143 + t140 * t105;
t108 = -t147 * t186 + t152;
t161 = -t108 + t137;
t69 = t138 * t97 + t139 * (t139 * t160 + t189);
t153 = (-qJD(5) + t136) * t91 - t162;
t146 = qJD(4) ^ 2;
t130 = 0.2e1 * t167;
t124 = t170 * qJDD(1);
t123 = t138 * t193;
t122 = t139 * t193;
t115 = -t154 + t177;
t112 = -t117 - t146;
t111 = -t117 + t146;
t110 = t116 - t146;
t104 = t118 - 0.2e1 * t169;
t102 = 0.2e1 * t168 + t78;
t98 = -t146 - t116;
t83 = -t88 + t132;
t82 = t87 - t132;
t80 = -t88 - t132;
t79 = -t116 - t117;
t76 = -t112 * t141 - t175;
t75 = t112 * t144 - t176;
t74 = t118 * t141 - t144 * t78;
t73 = -t118 * t144 - t141 * t78;
t72 = t144 * t98 - t200;
t71 = t141 * t98 + t199;
t68 = t88 - t87;
t64 = -t132 - t87;
t63 = (t140 * t91 - t143 * t89) * t136;
t62 = (-t140 * t89 - t143 * t91) * t136;
t57 = -qJD(5) * t91 - t162;
t56 = -t87 - t88;
t55 = t138 * t76 + t139 * t75;
t53 = t143 * t82 - t183;
t52 = -t140 * t83 + t196;
t51 = t140 * t82 + t179;
t50 = t143 * t83 + t197;
t49 = -t140 * t80 - t179;
t48 = t143 * t80 - t183;
t47 = t138 * t74 + t139 * t73;
t46 = t84 + t58;
t41 = (qJD(5) + t136) * t91 + t162;
t40 = t138 * t72 + t139 * t71;
t39 = t143 * t58 - t173 * t91;
t38 = t140 * t58 + t172 * t91;
t37 = -t140 * t57 + t172 * t89;
t36 = t143 * t57 + t173 * t89;
t35 = t143 * t64 - t197;
t34 = t140 * t64 + t196;
t29 = -pkin(7) * t48 - t180;
t28 = -t141 * t48 + t144 * t49;
t27 = t141 * t49 + t144 * t48;
t26 = -pkin(7) * t34 - t184;
t25 = t140 * t46 + t143 * t153;
t24 = -t140 * t195 - t143 * t41;
t23 = t140 * t153 - t143 * t46;
t22 = -t140 * t41 + t143 * t195;
t21 = -t141 * t34 + t144 * t35;
t20 = t141 * t35 + t144 * t34;
t16 = -pkin(4) * t195 + pkin(7) * t49 - t184;
t15 = t138 * t31 + t185;
t14 = -pkin(4) * t41 + pkin(7) * t35 + t180;
t13 = t138 * t28 + t139 * t27;
t12 = -t141 * t23 + t144 * t25;
t11 = t141 * t25 + t144 * t23;
t10 = t138 * t21 + t139 * t20;
t7 = pkin(4) * t54 + pkin(7) * t9;
t6 = -pkin(7) * t23 - t8;
t5 = -pkin(4) * t56 + pkin(7) * t25 + t9;
t4 = t11 * t139 + t12 * t138;
t3 = t144 * t9 - t188;
t2 = t141 * t9 + t187;
t1 = t138 * t3 + t139 * t2;
t17 = [0, 0, 0, 0, 0, qJDD(1), t163, t159, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t158 - 0.2e1 * t177, t130 + 0.2e1 * t137 - t159, pkin(1) * t115 + qJ(2) * (-t147 * pkin(1) + t130 - t155), t135 * qJDD(1), -0.2e1 * t138 * t165, 0, t134 * qJDD(1), 0, 0, -t123 * t186 + t138 * t161, -t122 * t186 + t139 * t161, -qJ(2) * t193 - t124 * t186 - t69, -qJ(2) * t108 + t186 * t69, t139 * (t105 * t144 - t141 * t168) - t138 * (t105 * t141 + t144 * t168), t139 * (-t102 * t144 - t104 * t141) - t138 * (-t102 * t141 + t104 * t144), t139 * (-t111 * t141 + t199) - t138 * (t111 * t144 + t200), t139 * (-t103 * t141 + t144 * t169) - t138 * (t103 * t144 + t141 * t169), t139 * (t110 * t144 - t176) - t138 * (t110 * t141 + t175), (t139 * (-t119 * t144 + t121 * t141) - t138 * (-t119 * t141 - t121 * t144)) * qJD(4), t139 * (-pkin(6) * t71 - t182) - t138 * (-pkin(3) * t102 + pkin(6) * t72 + t178) + qJ(2) * t102 + t186 * t40, t139 * (-pkin(6) * t75 - t178) - t138 * (-pkin(3) * t104 + pkin(6) * t76 - t182) + qJ(2) * t104 + t186 * t55, t139 * (-pkin(6) * t73 - t30) - t138 * (-pkin(3) * t79 + pkin(6) * t74 + t31) + qJ(2) * t79 + t186 * t47, -pkin(6) * t185 - t138 * (pkin(3) * t85 + pkin(6) * t31) - qJ(2) * t85 + t186 * t15, t139 * (-t141 * t38 + t144 * t39) - t138 * (t141 * t39 + t144 * t38), t139 * (-t141 * t22 + t144 * t24) - t138 * (t141 * t24 + t144 * t22), t139 * (-t141 * t50 + t144 * t52) - t138 * (t141 * t52 + t144 * t50), t139 * (-t141 * t36 + t144 * t37) - t138 * (t141 * t37 + t144 * t36), t139 * (-t141 * t51 + t144 * t53) - t138 * (t141 * t53 + t144 * t51), t139 * (-t141 * t62 + t144 * t63) - t138 * (t141 * t63 + t144 * t62), t139 * (-pkin(6) * t20 - t14 * t141 + t144 * t26) - t138 * (-pkin(3) * t41 + pkin(6) * t21 + t14 * t144 + t141 * t26) + qJ(2) * t41 + t186 * t10, t139 * (-pkin(6) * t27 - t141 * t16 + t144 * t29) - t138 * (-pkin(3) * t195 + pkin(6) * t28 + t141 * t29 + t144 * t16) + qJ(2) * t195 + t186 * t13, t139 * (-pkin(6) * t11 - t141 * t5 + t144 * t6) - t138 * (-pkin(3) * t56 + pkin(6) * t12 + t141 * t6 + t144 * t5) + qJ(2) * t56 + t186 * t4, t139 * (-pkin(6) * t2 - pkin(7) * t187 - t141 * t7) - t138 * (pkin(3) * t54 + pkin(6) * t3 - pkin(7) * t188 + t144 * t7) - qJ(2) * t54 + t186 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t147, -t115, 0, 0, 0, 0, 0, 0, -t123, -t122, -t124, t69, 0, 0, 0, 0, 0, 0, t40, t55, t47, t15, 0, 0, 0, 0, 0, 0, t10, t13, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t165, -t193, -t108, 0, 0, 0, 0, 0, 0, t102, t104, t79, -t85, 0, 0, 0, 0, 0, 0, t41, t195, t56, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t117 - t116, t118, -t174, -t78, qJDD(4), -t60, -t61, 0, 0, t70, t68, t46, -t70, t153, t133, pkin(4) * t34 - t18, -t181 - t140 * t191 + (-t140 * t198 + t48) * pkin(4), pkin(4) * t23, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t68, t46, -t70, t153, t133, -t18, -t19, 0, 0;];
tauJ_reg = t17;
