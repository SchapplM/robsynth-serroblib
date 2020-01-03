% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:18
% EndTime: 2019-12-31 18:04:24
% DurationCPUTime: 1.98s
% Computational Cost: add. (5566->264), mult. (13730->364), div. (0->0), fcn. (9616->8), ass. (0->166)
t139 = sin(qJ(5));
t135 = qJDD(4) + qJDD(5);
t137 = sin(pkin(8));
t138 = cos(pkin(8));
t140 = sin(qJ(4));
t143 = cos(qJ(4));
t159 = t137 * t140 + t138 * t143;
t111 = t159 * qJD(1);
t176 = qJD(1) * t138;
t177 = qJD(1) * t137;
t112 = -t140 * t176 + t143 * t177;
t142 = cos(qJ(5));
t80 = t142 * t111 + t139 * t112;
t82 = -t139 * t111 + t142 * t112;
t60 = t82 * t80;
t203 = -t60 + t135;
t209 = t139 * t203;
t94 = t112 * t111;
t201 = qJDD(4) - t94;
t208 = t140 * t201;
t207 = t142 * t203;
t206 = t143 * t201;
t141 = sin(qJ(1));
t197 = cos(qJ(1));
t156 = t197 * g(1) + t141 * g(2);
t198 = qJD(1) ^ 2;
t205 = -t198 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t156;
t173 = t112 * qJD(4);
t70 = t159 * qJDD(1);
t152 = t70 + t173;
t130 = t137 * qJDD(1);
t169 = t138 * qJDD(1);
t110 = t143 * t130 - t140 * t169;
t174 = t111 * qJD(4);
t92 = t110 - t174;
t52 = -t80 * qJD(5) - t139 * t152 + t142 * t92;
t136 = qJD(4) + qJD(5);
t77 = t136 * t80;
t204 = -t77 + t52;
t145 = t137 ^ 2;
t147 = t138 ^ 2;
t202 = t145 + t147;
t119 = t202 * t198;
t172 = t138 * t198;
t180 = t137 * qJ(3);
t196 = t138 * pkin(2);
t161 = -t180 - t196;
t117 = t161 * qJD(1);
t95 = t138 * g(3) + t205 * t137;
t71 = t117 * t177 + qJDD(3) + t95;
t150 = (-pkin(3) * t172 - pkin(6) * qJDD(1)) * t137 + t71;
t163 = -t137 * g(3) + t205 * t138;
t158 = t117 * t176 + t163;
t170 = t147 * t198;
t68 = -pkin(3) * t170 - pkin(6) * t169 + t158;
t48 = t140 * t68 - t143 * t150;
t200 = -t48 + (-t92 - t174) * pkin(7);
t78 = t80 ^ 2;
t79 = t82 ^ 2;
t199 = t111 ^ 2;
t109 = t112 ^ 2;
t134 = t136 ^ 2;
t149 = pkin(4) * t201 + t200;
t49 = t140 * t150 + t143 * t68;
t96 = qJD(4) * pkin(4) - t112 * pkin(7);
t28 = -t199 * pkin(4) - t152 * pkin(7) - qJD(4) * t96 + t49;
t13 = t139 * t28 - t142 * t149;
t189 = t142 * t28;
t14 = t139 * t149 + t189;
t6 = -t142 * t13 + t139 * t14;
t195 = t140 * t6;
t194 = t143 * t6;
t166 = t141 * g(1) - t197 * g(2);
t162 = -qJDD(2) + t166;
t151 = (qJ(2) * qJD(1) + 0.2e1 * qJD(3) * t137) * qJD(1) + t162;
t155 = t180 + t138 * (pkin(2) + pkin(3)) + pkin(1);
t171 = t145 * t198;
t69 = t155 * qJDD(1) + t151 + (-t170 - t171) * pkin(6);
t41 = t152 * pkin(4) - t199 * pkin(7) + t112 * t96 + t69;
t193 = t139 * t41;
t57 = t60 + t135;
t192 = t139 * t57;
t191 = t140 * t69;
t88 = qJDD(4) + t94;
t190 = t140 * t88;
t188 = t142 * t41;
t187 = t142 * t57;
t186 = t143 * t69;
t185 = t143 * t88;
t184 = qJ(2) * t137 * t119;
t183 = qJDD(1) * pkin(1);
t182 = t136 * t139;
t181 = t136 * t142;
t129 = t145 * qJDD(1);
t131 = t147 * qJDD(1);
t179 = qJ(2) * (t131 + t129) + pkin(1) * t119;
t175 = t198 * qJ(2);
t178 = -t202 * t138 * t175 + pkin(1) * t169;
t167 = t137 * t169;
t7 = t139 * t13 + t142 * t14;
t165 = t137 * t95 + t138 * t163;
t164 = t139 * t92 + t142 * t152;
t25 = t140 * t49 - t143 * t48;
t26 = t140 * t48 + t143 * t49;
t157 = pkin(1) - t161;
t154 = t157 * qJDD(1);
t153 = (-qJD(5) + t136) * t82 - t164;
t144 = qJD(4) ^ 2;
t108 = t162 + t175 + t183;
t99 = -t109 - t144;
t98 = -t109 + t144;
t97 = -t144 + t199;
t91 = t110 - 0.2e1 * t174;
t90 = 0.2e1 * t173 + t70;
t86 = -t144 - t199;
t84 = t154 + t151;
t75 = -t79 + t134;
t74 = t78 - t134;
t73 = -t79 - t134;
t72 = -t109 - t199;
t67 = -t140 * t99 - t185;
t66 = t143 * t99 - t190;
t64 = t140 * t110 - t143 * t70;
t63 = -t143 * t110 - t140 * t70;
t62 = t143 * t86 - t208;
t61 = t140 * t86 + t206;
t59 = t79 - t78;
t55 = -t134 - t78;
t54 = (t139 * t82 - t142 * t80) * t136;
t53 = (-t139 * t80 - t142 * t82) * t136;
t51 = -t82 * qJD(5) - t164;
t50 = -t78 - t79;
t47 = t142 * t74 - t192;
t46 = -t139 * t75 + t207;
t45 = t139 * t74 + t187;
t44 = t142 * t75 + t209;
t43 = -t139 * t73 - t187;
t42 = t142 * t73 - t192;
t40 = t77 + t52;
t35 = (qJD(5) + t136) * t82 + t164;
t34 = t142 * t52 - t82 * t182;
t33 = t139 * t52 + t82 * t181;
t32 = -t139 * t51 + t80 * t181;
t31 = t142 * t51 + t80 * t182;
t30 = t142 * t55 - t209;
t29 = t139 * t55 + t207;
t24 = -t140 * t42 + t143 * t43;
t23 = t140 * t43 + t143 * t42;
t22 = -pkin(7) * t42 + t188;
t21 = -pkin(7) * t29 + t193;
t20 = t139 * t40 + t142 * t153;
t19 = -t139 * t204 - t142 * t35;
t18 = t139 * t153 - t142 * t40;
t17 = -t139 * t35 + t142 * t204;
t16 = -t140 * t29 + t143 * t30;
t15 = t140 * t30 + t143 * t29;
t11 = -pkin(4) * t204 + pkin(7) * t43 + t193;
t10 = -pkin(4) * t35 + pkin(7) * t30 - t188;
t9 = -t140 * t18 + t143 * t20;
t8 = t140 * t20 + t143 * t18;
t5 = -pkin(4) * t41 + pkin(7) * t7;
t4 = -pkin(7) * t18 - t6;
t3 = -pkin(4) * t50 + pkin(7) * t20 + t7;
t2 = t143 * t7 - t195;
t1 = t140 * t7 + t194;
t12 = [0, 0, 0, 0, 0, qJDD(1), t166, t156, 0, 0, t129, 0.2e1 * t167, 0, t131, 0, 0, t138 * t108 + t178, t184 + (-t108 - t183) * t137, t165 + t179, pkin(1) * t108 + qJ(2) * t165, t129, 0, -0.2e1 * t167, 0, 0, t131, ((pkin(1) + 0.2e1 * t180 + 0.2e1 * t196) * qJDD(1) + t151) * t138 + t178, t137 * (qJ(3) * t119 + t71) + t138 * (pkin(2) * t119 + t158) + t179, -t184 + (t151 + 0.2e1 * t154) * t137, qJ(2) * (t137 * t71 + t138 * t158) + t157 * t84, t137 * (-t140 * t173 + t143 * t92) + t138 * (-t140 * t92 - t143 * t173), t137 * (-t140 * t91 - t143 * t90) + t138 * (t140 * t90 - t143 * t91), t137 * (-t140 * t98 + t206) + t138 * (-t143 * t98 - t208), t137 * (t140 * t152 + t143 * t174) + t138 * (-t140 * t174 + t143 * t152), t137 * (t143 * t97 - t190) + t138 * (-t140 * t97 - t185), (t137 * (-t143 * t111 + t140 * t112) + t138 * (t140 * t111 + t143 * t112)) * qJD(4), t137 * (-pkin(6) * t61 + t191) + t138 * (-pkin(6) * t62 + t186) + qJ(2) * (t137 * t61 + t138 * t62) + t155 * t90, t137 * (-pkin(6) * t66 + t186) + t138 * (-pkin(6) * t67 - t191) + qJ(2) * (t137 * t66 + t138 * t67) + t155 * t91, t137 * (-pkin(6) * t63 - t25) + t138 * (-pkin(6) * t64 - t26) + qJ(2) * (t137 * t63 + t138 * t64) + t155 * t72, t155 * t69 + (-pkin(6) + qJ(2)) * (t137 * t25 + t138 * t26), t137 * (-t140 * t33 + t143 * t34) + t138 * (-t140 * t34 - t143 * t33), t137 * (-t140 * t17 + t143 * t19) + t138 * (-t140 * t19 - t143 * t17), t137 * (-t140 * t44 + t143 * t46) + t138 * (-t140 * t46 - t143 * t44), t137 * (-t140 * t31 + t143 * t32) + t138 * (-t140 * t32 - t143 * t31), t137 * (-t140 * t45 + t143 * t47) + t138 * (-t140 * t47 - t143 * t45), t137 * (-t140 * t53 + t143 * t54) + t138 * (-t140 * t54 - t143 * t53), t137 * (-pkin(6) * t15 - t140 * t10 + t143 * t21) + t138 * (-pkin(6) * t16 - t143 * t10 - t140 * t21) + qJ(2) * (t137 * t15 + t138 * t16) + t155 * t35, t137 * (-pkin(6) * t23 - t140 * t11 + t143 * t22) + t138 * (-pkin(6) * t24 - t143 * t11 - t140 * t22) + qJ(2) * (t137 * t23 + t138 * t24) + t155 * t204, t137 * (-pkin(6) * t8 - t140 * t3 + t143 * t4) + t138 * (-pkin(6) * t9 - t140 * t4 - t143 * t3) + qJ(2) * (t137 * t8 + t138 * t9) + t155 * t50, t137 * (-pkin(6) * t1 - pkin(7) * t194 - t140 * t5) + t138 * (-pkin(6) * t2 + pkin(7) * t195 - t143 * t5) + qJ(2) * (t137 * t1 + t138 * t2) + t155 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, t130, -t119, -t108, 0, 0, 0, 0, 0, 0, -t169, -t119, -t130, -t84, 0, 0, 0, 0, 0, 0, -t90, -t91, -t72, -t69, 0, 0, 0, 0, 0, 0, -t35, -t204, -t50, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137 * t172, t130, -t171, t71, 0, 0, 0, 0, 0, 0, t61, t66, t63, t25, 0, 0, 0, 0, 0, 0, t15, t23, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t109 - t199, t110, -t94, -t70, qJDD(4), -t48, -t49, 0, 0, t60, t59, t40, -t60, t153, t135, pkin(4) * t29 - t13, -t189 - t139 * t200 + (-t139 * t201 + t42) * pkin(4), pkin(4) * t18, pkin(4) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t59, t40, -t60, t153, t135, -t13, -t14, 0, 0;];
tauJ_reg = t12;
