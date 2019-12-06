% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:59
% EndTime: 2019-12-05 17:54:07
% DurationCPUTime: 2.40s
% Computational Cost: add. (8975->309), mult. (19759->447), div. (0->0), fcn. (13240->10), ass. (0->182)
t223 = -2 * qJD(4);
t156 = sin(pkin(9));
t158 = cos(pkin(9));
t164 = cos(qJ(3));
t162 = sin(qJ(3));
t190 = qJD(1) * t162;
t124 = -qJD(1) * t158 * t164 + t156 * t190;
t126 = (t156 * t164 + t158 * t162) * qJD(1);
t109 = t126 * t124;
t214 = qJDD(3) - t109;
t222 = t156 * t214;
t221 = t158 * t214;
t161 = sin(qJ(5));
t151 = qJDD(3) + qJDD(5);
t163 = cos(qJ(5));
t100 = -t124 * t161 + t126 * t163;
t98 = t124 * t163 + t126 * t161;
t77 = t100 * t98;
t215 = -t77 + t151;
t220 = t161 * t215;
t219 = t163 * t215;
t182 = qJD(1) * qJD(3);
t177 = t164 * t182;
t181 = t162 * qJDD(1);
t134 = t177 + t181;
t212 = qJD(1) ^ 2;
t183 = t164 * t212;
t144 = t162 * t183;
t139 = qJDD(3) + t144;
t210 = sin(qJ(1));
t211 = cos(qJ(1));
t169 = -g(2) * t210 + g(3) * t211;
t132 = -pkin(1) * t212 - t169;
t157 = sin(pkin(8));
t159 = cos(pkin(8));
t170 = g(2) * t211 + g(3) * t210;
t167 = qJDD(1) * pkin(1) + t170;
t191 = t132 * t159 + t157 * t167;
t102 = -pkin(2) * t212 + qJDD(1) * pkin(6) + t191;
t154 = -g(1) + qJDD(2);
t89 = t162 * t102 - t154 * t164;
t74 = (-t134 + t177) * qJ(4) + t139 * pkin(3) - t89;
t179 = qJ(4) * t190;
t138 = qJD(3) * pkin(3) - t179;
t193 = t162 * t154;
t76 = t193 + (-t138 - t179) * qJD(3) + (-pkin(3) * t183 + qJ(4) * qJDD(1) + t102) * t164;
t172 = t126 * t223 - t156 * t76 + t158 * t74;
t178 = t162 * t182;
t180 = t164 * qJDD(1);
t171 = -t178 + t180;
t111 = t158 * t134 + t156 * t171;
t189 = qJD(3) * t124;
t173 = -t111 - t189;
t218 = t173 * pkin(7) + t172;
t110 = -t134 * t156 + t158 * t171;
t61 = -qJD(5) * t98 + t110 * t161 + t111 * t163;
t152 = qJD(3) + qJD(5);
t95 = t152 * t98;
t216 = t61 - t95;
t37 = t124 * t223 + t156 * t74 + t158 * t76;
t96 = t98 ^ 2;
t97 = t100 ^ 2;
t122 = t124 ^ 2;
t123 = t126 ^ 2;
t150 = t152 ^ 2;
t213 = t164 ^ 2;
t166 = pkin(4) * t214 + t218;
t114 = qJD(3) * pkin(4) - pkin(7) * t126;
t33 = -pkin(4) * t122 + pkin(7) * t110 - qJD(3) * t114 + t37;
t16 = t161 * t33 - t163 * t166;
t201 = t163 * t33;
t17 = t161 * t166 + t201;
t8 = -t16 * t163 + t161 * t17;
t209 = t156 * t8;
t208 = t158 * t8;
t174 = -t157 * t132 + t159 * t167;
t101 = -qJDD(1) * pkin(2) - pkin(6) * t212 - t174;
t149 = t213 * t212;
t80 = -pkin(3) * t171 - qJ(4) * t149 + t138 * t190 + qJDD(4) + t101;
t206 = t156 * t80;
t205 = t158 * t80;
t41 = -pkin(4) * t110 - pkin(7) * t122 + t114 * t126 + t80;
t204 = t161 * t41;
t71 = t77 + t151;
t203 = t161 * t71;
t20 = t156 * t37 + t158 * t172;
t202 = t162 * t20;
t200 = t163 * t41;
t199 = t163 * t71;
t198 = t152 * t161;
t197 = t152 * t163;
t105 = qJDD(3) + t109;
t196 = t156 * t105;
t195 = t158 * t105;
t194 = t162 * t139;
t140 = qJDD(3) - t144;
t192 = t164 * t140;
t188 = qJD(3) * t126;
t187 = qJD(3) * t156;
t186 = qJD(3) * t158;
t9 = t16 * t161 + t163 * t17;
t21 = -t156 * t172 + t158 * t37;
t90 = t164 * t102 + t193;
t64 = t162 * t89 + t164 * t90;
t175 = -t110 * t163 + t161 * t111;
t86 = t110 + t188;
t135 = -0.2e1 * t178 + t180;
t168 = (-qJD(5) + t152) * t100 - t175;
t165 = qJD(3) ^ 2;
t153 = t162 ^ 2;
t148 = t153 * t212;
t143 = -t149 - t165;
t142 = -t148 - t165;
t137 = t148 + t149;
t136 = (t153 + t213) * qJDD(1);
t133 = 0.2e1 * t177 + t181;
t117 = -t123 - t165;
t116 = -t123 + t165;
t115 = t122 - t165;
t113 = -t142 * t162 - t192;
t112 = t143 * t164 - t194;
t103 = -t165 - t122;
t93 = -t97 + t150;
t92 = t96 - t150;
t91 = -t97 - t150;
t87 = t111 - t189;
t85 = -t110 + t188;
t84 = -t122 - t123;
t82 = -t117 * t156 - t195;
t81 = t117 * t158 - t196;
t79 = t103 * t158 - t222;
t78 = t103 * t156 + t221;
t75 = t97 - t96;
t69 = -t150 - t96;
t66 = (t100 * t161 - t163 * t98) * t152;
t65 = (-t100 * t163 - t161 * t98) * t152;
t63 = -t156 * t173 + t158 * t86;
t62 = t156 * t86 + t158 * t173;
t60 = -qJD(5) * t100 - t175;
t59 = -t96 - t97;
t58 = -t162 * t81 + t164 * t82;
t57 = t163 * t92 - t203;
t56 = -t161 * t93 + t219;
t55 = t161 * t92 + t199;
t54 = t163 * t93 + t220;
t53 = -t161 * t91 - t199;
t52 = t163 * t91 - t203;
t51 = t61 + t95;
t46 = (qJD(5) + t152) * t100 + t175;
t45 = -t100 * t198 + t163 * t61;
t44 = t100 * t197 + t161 * t61;
t43 = -t161 * t60 + t197 * t98;
t42 = t163 * t60 + t198 * t98;
t40 = -t162 * t78 + t164 * t79;
t39 = t163 * t69 - t220;
t38 = t161 * t69 + t219;
t34 = -t162 * t62 + t164 * t63;
t31 = -t156 * t52 + t158 * t53;
t30 = t156 * t53 + t158 * t52;
t29 = -pkin(7) * t52 + t200;
t28 = t161 * t51 + t163 * t168;
t27 = -t161 * t216 - t163 * t46;
t26 = t161 * t168 - t163 * t51;
t25 = -t161 * t46 + t163 * t216;
t24 = -pkin(7) * t38 + t204;
t23 = -t156 * t38 + t158 * t39;
t22 = t156 * t39 + t158 * t38;
t19 = -pkin(4) * t216 + pkin(7) * t53 + t204;
t18 = -pkin(4) * t46 + pkin(7) * t39 - t200;
t14 = -t162 * t30 + t164 * t31;
t13 = -t156 * t26 + t158 * t28;
t12 = t156 * t28 + t158 * t26;
t11 = -t162 * t22 + t164 * t23;
t10 = t164 * t21 - t202;
t7 = -pkin(4) * t41 + pkin(7) * t9;
t6 = -t12 * t162 + t13 * t164;
t5 = -pkin(7) * t26 - t8;
t4 = -pkin(4) * t59 + pkin(7) * t28 + t9;
t3 = t158 * t9 - t209;
t2 = t156 * t9 + t208;
t1 = -t162 * t2 + t164 * t3;
t15 = [0, 0, 0, 0, 0, qJDD(1), t170, t169, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t159 - t157 * t212) + t174, pkin(1) * (-qJDD(1) * t157 - t159 * t212) - t191, 0, pkin(1) * (t157 * t191 + t159 * t174), (t134 + t177) * t162, t133 * t164 + t135 * t162, t194 + t164 * (-t148 + t165), t135 * t164, t162 * (t149 - t165) + t192, 0, -t164 * t101 + pkin(2) * t135 + pkin(6) * t112 + pkin(1) * (t112 * t157 + t135 * t159), t162 * t101 - pkin(2) * t133 + pkin(6) * t113 + pkin(1) * (t113 * t157 - t133 * t159), pkin(2) * t137 + pkin(6) * t136 + pkin(1) * (t136 * t157 + t137 * t159) + t64, -pkin(2) * t101 + pkin(6) * t64 + pkin(1) * (-t101 * t159 + t157 * t64), t162 * (t111 * t158 - t126 * t187) + t164 * (t111 * t156 + t126 * t186), t162 * (-t156 * t87 - t158 * t85) + t164 * (-t156 * t85 + t158 * t87), t162 * (-t116 * t156 + t221) + t164 * (t116 * t158 + t222), t162 * (-t110 * t156 + t124 * t186) + t164 * (t110 * t158 + t124 * t187), t162 * (t115 * t158 - t196) + t164 * (t115 * t156 + t195), (t162 * (-t124 * t158 + t126 * t156) + t164 * (-t124 * t156 - t126 * t158)) * qJD(3), t162 * (-qJ(4) * t78 + t206) + t164 * (-pkin(3) * t85 + qJ(4) * t79 - t205) - pkin(2) * t85 + pkin(6) * t40 + pkin(1) * (t157 * t40 - t159 * t85), t162 * (-qJ(4) * t81 + t205) + t164 * (-pkin(3) * t87 + qJ(4) * t82 + t206) - pkin(2) * t87 + pkin(6) * t58 + pkin(1) * (t157 * t58 - t159 * t87), t162 * (-qJ(4) * t62 - t20) + t164 * (-pkin(3) * t84 + qJ(4) * t63 + t21) - pkin(2) * t84 + pkin(6) * t34 + pkin(1) * (t157 * t34 - t159 * t84), -qJ(4) * t202 + t164 * (-pkin(3) * t80 + qJ(4) * t21) - pkin(2) * t80 + pkin(6) * t10 + pkin(1) * (t10 * t157 - t159 * t80), t162 * (-t156 * t44 + t158 * t45) + t164 * (t156 * t45 + t158 * t44), t162 * (-t156 * t25 + t158 * t27) + t164 * (t156 * t27 + t158 * t25), t162 * (-t156 * t54 + t158 * t56) + t164 * (t156 * t56 + t158 * t54), t162 * (-t156 * t42 + t158 * t43) + t164 * (t156 * t43 + t158 * t42), t162 * (-t156 * t55 + t158 * t57) + t164 * (t156 * t57 + t158 * t55), t162 * (-t156 * t65 + t158 * t66) + t164 * (t156 * t66 + t158 * t65), t162 * (-qJ(4) * t22 - t156 * t18 + t158 * t24) + t164 * (-pkin(3) * t46 + qJ(4) * t23 + t156 * t24 + t158 * t18) - pkin(2) * t46 + pkin(6) * t11 + pkin(1) * (t11 * t157 - t159 * t46), t162 * (-qJ(4) * t30 - t156 * t19 + t158 * t29) + t164 * (-pkin(3) * t216 + qJ(4) * t31 + t156 * t29 + t158 * t19) - pkin(2) * t216 + pkin(6) * t14 + pkin(1) * (t14 * t157 - t159 * t216), t162 * (-qJ(4) * t12 - t156 * t4 + t158 * t5) + t164 * (-pkin(3) * t59 + qJ(4) * t13 + t156 * t5 + t158 * t4) - pkin(2) * t59 + pkin(6) * t6 + pkin(1) * (t157 * t6 - t159 * t59), t162 * (-pkin(7) * t208 - qJ(4) * t2 - t156 * t7) + t164 * (-pkin(3) * t41 - pkin(7) * t209 + qJ(4) * t3 + t158 * t7) - pkin(2) * t41 + pkin(6) * t1 + pkin(1) * (t1 * t157 - t159 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, 0, 0, 0, 0, 0, 0, t139 * t164 + t143 * t162, -t140 * t162 + t142 * t164, 0, t162 * t90 - t164 * t89, 0, 0, 0, 0, 0, 0, t162 * t79 + t164 * t78, t162 * t82 + t164 * t81, t162 * t63 + t164 * t62, t162 * t21 + t164 * t20, 0, 0, 0, 0, 0, 0, t162 * t23 + t164 * t22, t162 * t31 + t164 * t30, t12 * t164 + t13 * t162, t162 * t3 + t164 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t148 - t149, t181, t144, t180, qJDD(3), -t89, -t90, 0, 0, t109, t123 - t122, -t173, -t109, t86, qJDD(3), pkin(3) * t78 + t172, pkin(3) * t81 - t37, pkin(3) * t62, pkin(3) * t20, t77, t75, t51, -t77, t168, t151, pkin(3) * t22 + pkin(4) * t38 - t16, -t201 - t161 * t218 + pkin(3) * t30 + (-t161 * t214 + t52) * pkin(4), pkin(3) * t12 + pkin(4) * t26, pkin(3) * t2 + pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t87, t84, t80, 0, 0, 0, 0, 0, 0, t46, t216, t59, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t75, t51, -t77, t168, t151, -t16, -t17, 0, 0;];
tauJ_reg = t15;
