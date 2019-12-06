% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:08
% DurationCPUTime: 2.27s
% Computational Cost: add. (8492->289), mult. (17967->407), div. (0->0), fcn. (12794->10), ass. (0->182)
t166 = sin(qJ(5));
t159 = qJDD(3) + qJDD(4);
t154 = qJDD(5) + t159;
t167 = sin(qJ(4));
t171 = cos(qJ(4));
t172 = cos(qJ(3));
t168 = sin(qJ(3));
t189 = qJD(2) * t168;
t132 = -t171 * t172 * qJD(2) + t167 * t189;
t134 = (t172 * t167 + t168 * t171) * qJD(2);
t170 = cos(qJ(5));
t107 = t170 * t132 + t166 * t134;
t109 = -t166 * t132 + t170 * t134;
t80 = t109 * t107;
t213 = -t80 + t154;
t219 = t166 * t213;
t115 = t134 * t132;
t211 = -t115 + t159;
t218 = t167 * t211;
t217 = t170 * t213;
t216 = t171 * t211;
t160 = qJD(3) + qJD(4);
t155 = qJD(5) + t160;
t101 = t155 * t107;
t187 = qJD(2) * qJD(3);
t182 = t172 * t187;
t186 = t168 * qJDD(2);
t139 = t182 + t186;
t183 = t168 * t187;
t185 = t172 * qJDD(2);
t177 = -t183 + t185;
t180 = t167 * t139 - t171 * t177;
t95 = -t134 * qJD(4) - t180;
t96 = -t132 * qJD(4) + t171 * t139 + t167 * t177;
t57 = -t107 * qJD(5) + t166 * t95 + t170 * t96;
t215 = -t101 + t57;
t127 = t160 * t132;
t214 = -t127 + t96;
t212 = -t96 - t127;
t105 = t107 ^ 2;
t106 = t109 ^ 2;
t130 = t132 ^ 2;
t131 = t134 ^ 2;
t153 = t155 ^ 2;
t158 = t160 ^ 2;
t210 = t172 ^ 2;
t209 = qJD(2) ^ 2;
t163 = sin(pkin(9));
t164 = cos(pkin(9));
t144 = -t164 * g(1) - t163 * g(2);
t169 = sin(qJ(2));
t173 = cos(qJ(2));
t190 = -g(3) + qJDD(1);
t125 = t173 * t144 + t169 * t190;
t117 = -t209 * pkin(2) + qJDD(2) * pkin(6) + t125;
t143 = -t163 * g(1) + t164 * g(2);
t103 = t168 * t117 - t172 * t143;
t188 = t172 * t209;
t148 = t168 * t188;
t145 = qJDD(3) + t148;
t81 = (-t139 + t182) * pkin(7) + t145 * pkin(3) - t103;
t184 = pkin(7) * t189;
t147 = qJD(3) * pkin(3) - t184;
t194 = t168 * t143;
t82 = t194 + (-t147 - t184) * qJD(3) + (-pkin(3) * t188 + qJDD(2) * pkin(7) + t117) * t172;
t54 = t167 * t82 - t171 * t81;
t36 = t211 * pkin(4) + t212 * pkin(8) - t54;
t121 = t160 * pkin(4) - t134 * pkin(8);
t55 = t167 * t81 + t171 * t82;
t37 = -t130 * pkin(4) + t95 * pkin(8) - t160 * t121 + t55;
t18 = t166 * t37 - t170 * t36;
t19 = t166 * t36 + t170 * t37;
t9 = t166 * t19 - t170 * t18;
t208 = t167 * t9;
t207 = t171 * t9;
t124 = -t169 * t144 + t173 * t190;
t116 = -qJDD(2) * pkin(2) - t209 * pkin(6) - t124;
t157 = t210 * t209;
t91 = -t177 * pkin(3) - pkin(7) * t157 + t147 * t189 + t116;
t52 = -t95 * pkin(4) - t130 * pkin(8) + t134 * t121 + t91;
t206 = t166 * t52;
t74 = t80 + t154;
t205 = t166 * t74;
t204 = t167 * t91;
t29 = t167 * t55 - t171 * t54;
t203 = t168 * t29;
t202 = t170 * t52;
t201 = t170 * t74;
t200 = t171 * t91;
t199 = t155 * t166;
t198 = t155 * t170;
t197 = t160 * t167;
t196 = t160 * t171;
t112 = t115 + t159;
t195 = t167 * t112;
t193 = t168 * t145;
t192 = t171 * t112;
t191 = t172 * (qJDD(3) - t148);
t10 = t166 * t18 + t170 * t19;
t181 = t166 * t96 - t170 * t95;
t30 = t167 * t54 + t171 * t55;
t104 = t172 * t117 + t194;
t69 = t168 * t103 + t172 * t104;
t72 = -t153 - t105;
t49 = t166 * t72 + t217;
t179 = pkin(4) * t49 - t18;
t93 = -t106 - t153;
t62 = t170 * t93 - t205;
t178 = pkin(4) * t62 - t19;
t140 = -0.2e1 * t183 + t185;
t176 = (-qJD(5) + t155) * t109 - t181;
t175 = (-qJD(4) + t160) * t134 - t180;
t174 = qJD(3) ^ 2;
t161 = t168 ^ 2;
t156 = t161 * t209;
t142 = t156 + t157;
t141 = (t161 + t210) * qJDD(2);
t138 = 0.2e1 * t182 + t186;
t123 = -t131 + t158;
t122 = t130 - t158;
t120 = -t131 - t158;
t119 = -t191 - t168 * (-t156 - t174);
t118 = t172 * (-t157 - t174) - t193;
t114 = t131 - t130;
t110 = -t158 - t130;
t100 = -t106 + t153;
t99 = t105 - t153;
t97 = -t130 - t131;
t90 = -t167 * t120 - t192;
t89 = t171 * t120 - t195;
t83 = (qJD(4) + t160) * t134 + t180;
t79 = t106 - t105;
t77 = t171 * t110 - t218;
t76 = t167 * t110 + t216;
t71 = (-t107 * t170 + t109 * t166) * t155;
t70 = (-t107 * t166 - t109 * t170) * t155;
t68 = -t105 - t106;
t67 = t170 * t99 - t205;
t66 = -t166 * t100 + t217;
t65 = t166 * t99 + t201;
t64 = t170 * t100 + t219;
t63 = -t166 * t93 - t201;
t60 = -t168 * t89 + t172 * t90;
t59 = -t167 * t212 + t171 * t175;
t58 = t167 * t175 + t171 * t212;
t56 = -t109 * qJD(5) - t181;
t51 = -t168 * t76 + t172 * t77;
t50 = t170 * t72 - t219;
t46 = t101 + t57;
t42 = (qJD(5) + t155) * t109 + t181;
t41 = -t109 * t199 + t170 * t57;
t40 = t109 * t198 + t166 * t57;
t39 = t107 * t198 - t166 * t56;
t38 = t107 * t199 + t170 * t56;
t34 = -t167 * t62 + t171 * t63;
t33 = t167 * t63 + t171 * t62;
t32 = -t168 * t58 + t172 * t59;
t31 = -pkin(8) * t62 + t202;
t28 = -pkin(8) * t49 + t206;
t27 = -t167 * t49 + t171 * t50;
t26 = t167 * t50 + t171 * t49;
t25 = -t166 * t215 - t170 * t42;
t24 = t166 * t46 + t170 * t176;
t23 = -t166 * t42 + t170 * t215;
t22 = t166 * t176 - t170 * t46;
t21 = pkin(4) * t22;
t20 = -pkin(4) * t215 + pkin(8) * t63 + t206;
t16 = -pkin(4) * t42 + pkin(8) * t50 - t202;
t15 = -t168 * t33 + t172 * t34;
t14 = t172 * t30 - t203;
t13 = -t168 * t26 + t172 * t27;
t12 = -t167 * t22 + t171 * t24;
t11 = t167 * t24 + t171 * t22;
t8 = pkin(4) * t9;
t7 = -pkin(4) * t52 + pkin(8) * t10;
t6 = -pkin(8) * t22 - t9;
t5 = -pkin(4) * t68 + pkin(8) * t24 + t10;
t4 = -t168 * t11 + t172 * t12;
t3 = t171 * t10 - t208;
t2 = t167 * t10 + t207;
t1 = -t168 * t2 + t172 * t3;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t190, 0, 0, 0, 0, 0, 0, t173 * qJDD(2) - t169 * t209, -t169 * qJDD(2) - t173 * t209, 0, t173 * t124 + t169 * t125, 0, 0, 0, 0, 0, 0, t169 * t118 + t173 * t140, t169 * t119 - t173 * t138, t169 * t141 + t173 * t142, -t173 * t116 + t169 * t69, 0, 0, 0, 0, 0, 0, t169 * t51 - t173 * t83, t169 * t60 - t173 * t214, t169 * t32 - t173 * t97, t169 * t14 - t173 * t91, 0, 0, 0, 0, 0, 0, t169 * t13 - t173 * t42, t169 * t15 - t173 * t215, t169 * t4 - t173 * t68, t169 * t1 - t173 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t124, -t125, 0, 0, (t139 + t182) * t168, t172 * t138 + t168 * t140, t193 + t172 * (-t156 + t174), t140 * t172, t168 * (t157 - t174) + t191, 0, pkin(2) * t140 + pkin(6) * t118 - t172 * t116, -pkin(2) * t138 + pkin(6) * t119 + t168 * t116, pkin(2) * t142 + pkin(6) * t141 + t69, -pkin(2) * t116 + pkin(6) * t69, t168 * (-t134 * t197 + t171 * t96) + t172 * (t134 * t196 + t167 * t96), t168 * (-t167 * t214 - t171 * t83) + t172 * (-t167 * t83 + t171 * t214), t168 * (-t167 * t123 + t216) + t172 * (t171 * t123 + t218), t168 * (t132 * t196 - t167 * t95) + t172 * (t132 * t197 + t171 * t95), t168 * (t171 * t122 - t195) + t172 * (t167 * t122 + t192), (t168 * (-t132 * t171 + t134 * t167) + t172 * (-t132 * t167 - t134 * t171)) * t160, t168 * (-pkin(7) * t76 + t204) + t172 * (-pkin(3) * t83 + pkin(7) * t77 - t200) - pkin(2) * t83 + pkin(6) * t51, t168 * (-pkin(7) * t89 + t200) + t172 * (-pkin(3) * t214 + pkin(7) * t90 + t204) - pkin(2) * t214 + pkin(6) * t60, t168 * (-pkin(7) * t58 - t29) + t172 * (-pkin(3) * t97 + pkin(7) * t59 + t30) - pkin(2) * t97 + pkin(6) * t32, -pkin(7) * t203 + t172 * (-pkin(3) * t91 + pkin(7) * t30) - pkin(2) * t91 + pkin(6) * t14, t168 * (-t167 * t40 + t171 * t41) + t172 * (t167 * t41 + t171 * t40), t168 * (-t167 * t23 + t171 * t25) + t172 * (t167 * t25 + t171 * t23), t168 * (-t167 * t64 + t171 * t66) + t172 * (t167 * t66 + t171 * t64), t168 * (-t167 * t38 + t171 * t39) + t172 * (t167 * t39 + t171 * t38), t168 * (-t167 * t65 + t171 * t67) + t172 * (t167 * t67 + t171 * t65), t168 * (-t167 * t70 + t171 * t71) + t172 * (t167 * t71 + t171 * t70), t168 * (-pkin(7) * t26 - t167 * t16 + t171 * t28) + t172 * (-pkin(3) * t42 + pkin(7) * t27 + t171 * t16 + t167 * t28) - pkin(2) * t42 + pkin(6) * t13, t168 * (-pkin(7) * t33 - t167 * t20 + t171 * t31) + t172 * (-pkin(3) * t215 + pkin(7) * t34 + t167 * t31 + t171 * t20) - pkin(2) * t215 + pkin(6) * t15, t168 * (-pkin(7) * t11 - t167 * t5 + t171 * t6) + t172 * (-pkin(3) * t68 + pkin(7) * t12 + t167 * t6 + t171 * t5) - pkin(2) * t68 + pkin(6) * t4, t168 * (-pkin(7) * t2 - pkin(8) * t207 - t167 * t7) + t172 * (-pkin(3) * t52 + pkin(7) * t3 - pkin(8) * t208 + t171 * t7) - pkin(2) * t52 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t156 - t157, t186, t148, t185, qJDD(3), -t103, -t104, 0, 0, t115, t114, -t212, -t115, t175, t159, pkin(3) * t76 - t54, pkin(3) * t89 - t55, pkin(3) * t58, pkin(3) * t29, t80, t79, t46, -t80, t176, t154, pkin(3) * t26 + t179, pkin(3) * t33 + t178, pkin(3) * t11 + t21, pkin(3) * t2 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t114, -t212, -t115, t175, t159, -t54, -t55, 0, 0, t80, t79, t46, -t80, t176, t154, t179, t178, t21, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t79, t46, -t80, t176, t154, -t18, -t19, 0, 0;];
tauJ_reg = t17;
