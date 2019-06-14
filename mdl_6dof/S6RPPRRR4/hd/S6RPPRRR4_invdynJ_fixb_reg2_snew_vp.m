% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:41:35
% EndTime: 2019-05-05 15:41:45
% DurationCPUTime: 3.29s
% Computational Cost: add. (16193->352), mult. (29858->485), div. (0->0), fcn. (17344->10), ass. (0->216)
t186 = sin(qJ(6));
t188 = sin(qJ(4));
t218 = qJD(1) * qJD(4);
t168 = t188 * t218;
t192 = cos(qJ(4));
t171 = t192 * qJDD(1);
t153 = -t171 + t168;
t146 = qJDD(5) - t153;
t143 = qJDD(6) + t146;
t187 = sin(qJ(5));
t191 = cos(qJ(5));
t222 = qJD(1) * t188;
t147 = t191 * qJD(4) + t187 * t222;
t148 = -t187 * qJD(4) + t191 * t222;
t190 = cos(qJ(6));
t123 = -t190 * t147 - t186 * t148;
t125 = t186 * t147 - t190 * t148;
t95 = t125 * t123;
t248 = -t95 + t143;
t253 = t186 * t248;
t129 = t147 * t148;
t247 = -t129 + t146;
t252 = t187 * t247;
t251 = t190 * t248;
t250 = t191 * t247;
t166 = t192 * qJD(1) + qJD(5);
t161 = qJD(6) + t166;
t111 = t161 * t123;
t214 = t192 * t218;
t217 = t188 * qJDD(1);
t151 = -t214 - t217;
t204 = -t187 * qJDD(4) - t191 * t151;
t119 = t147 * qJD(5) - t204;
t209 = -t191 * qJDD(4) + t187 * t151;
t199 = t148 * qJD(5) - t209;
t73 = -t123 * qJD(6) + t190 * t119 + t186 * t199;
t249 = -t111 + t73;
t139 = t166 * t147;
t102 = t119 - t139;
t211 = t186 * t119 - t190 * t199;
t56 = (qJD(6) - t161) * t125 + t211;
t99 = (qJD(5) - t166) * t148 - t209;
t121 = t123 ^ 2;
t122 = t125 ^ 2;
t246 = t147 ^ 2;
t145 = t148 ^ 2;
t160 = t161 ^ 2;
t164 = t166 ^ 2;
t194 = qJD(1) ^ 2;
t245 = qJD(4) ^ 2;
t205 = -t153 - t168;
t206 = -t151 + t214;
t178 = qJDD(1) * qJ(2);
t189 = sin(qJ(1));
t193 = cos(qJ(1));
t207 = t193 * g(1) + t189 * g(2);
t202 = 0.2e1 * qJD(2) * qJD(1) - t207;
t198 = t178 + t202;
t243 = pkin(1) + pkin(2);
t136 = -t243 * t194 + t198;
t182 = sin(pkin(10));
t183 = cos(pkin(10));
t223 = t189 * g(1) - t193 * g(2);
t213 = -qJDD(2) + t223;
t201 = -t194 * qJ(2) - t213;
t195 = -t243 * qJDD(1) + t201;
t210 = t182 * t136 - t183 * t195;
t96 = qJDD(1) * pkin(3) - t194 * pkin(7) + t210;
t78 = t205 * pkin(4) + t206 * pkin(8) + t96;
t208 = pkin(4) * t192 + pkin(8) * t188;
t200 = t194 * t208;
t224 = g(3) + qJDD(3);
t106 = t183 * t136 + t182 * t195;
t97 = -t194 * pkin(3) - qJDD(1) * pkin(7) + t106;
t94 = t188 * t224 + t192 * t97;
t81 = -t245 * pkin(4) + qJDD(4) * pkin(8) - t192 * t200 + t94;
t45 = t187 * t81 - t191 * t78;
t36 = t247 * pkin(5) - t102 * pkin(9) - t45;
t135 = t166 * pkin(5) + t148 * pkin(9);
t46 = t187 * t78 + t191 * t81;
t40 = -t246 * pkin(5) + t199 * pkin(9) - t166 * t135 + t46;
t18 = t186 * t40 - t190 * t36;
t19 = t186 * t36 + t190 * t40;
t9 = -t190 * t18 + t186 * t19;
t244 = pkin(5) * t9;
t59 = t111 + t73;
t32 = -t186 * t56 - t190 * t59;
t242 = pkin(5) * t32;
t241 = t187 * t9;
t240 = t191 * t9;
t93 = t188 * t97 - t192 * t224;
t80 = -qJDD(4) * pkin(4) - t245 * pkin(8) - t188 * t200 + t93;
t47 = -t199 * pkin(5) - t246 * pkin(9) - t148 * t135 + t80;
t239 = t186 * t47;
t85 = t95 + t143;
t238 = t186 * t85;
t237 = t187 * t80;
t236 = t190 * t47;
t235 = t190 * t85;
t234 = t191 * t80;
t233 = qJDD(1) * pkin(1);
t232 = t161 * t186;
t231 = t161 * t190;
t230 = t166 * t187;
t229 = t166 * t191;
t115 = t129 + t146;
t228 = t187 * t115;
t165 = t192 * t194 * t188;
t158 = qJDD(4) + t165;
t227 = t188 * t158;
t226 = t191 * t115;
t159 = qJDD(4) - t165;
t225 = t192 * t159;
t220 = qJD(5) + t166;
t216 = t192 * t95;
t215 = t192 * t129;
t212 = qJ(2) * t183 - pkin(7);
t10 = t186 * t18 + t190 * t19;
t26 = t187 * t45 + t191 * t46;
t25 = t187 * t46 - t191 * t45;
t62 = t188 * t93 + t192 * t94;
t89 = -t160 - t121;
t49 = t186 * t89 + t251;
t203 = pkin(5) * t49 - t18;
t104 = -t122 - t160;
t63 = t190 * t104 - t238;
t197 = pkin(5) * t63 - t19;
t196 = qJ(2) * t182 + pkin(3) + t208;
t180 = t192 ^ 2;
t179 = t188 ^ 2;
t174 = t180 * t194;
t173 = t179 * t194;
t163 = -t174 - t245;
t162 = -t173 - t245;
t157 = t173 + t174;
t156 = (-t179 - t180) * qJDD(1);
t155 = t183 * qJDD(1) + t182 * t194;
t154 = -t182 * qJDD(1) + t183 * t194;
t152 = -t171 + 0.2e1 * t168;
t150 = 0.2e1 * t214 + t217;
t140 = -t201 + t233;
t138 = -t145 + t164;
t137 = -t164 + t246;
t132 = -t188 * t162 - t225;
t131 = t192 * t163 - t227;
t128 = t145 - t246;
t127 = -t145 - t164;
t126 = t182 * t156 + t183 * t157;
t120 = -t164 - t246;
t113 = t145 + t246;
t110 = -t122 + t160;
t109 = t121 - t160;
t108 = t182 * t132 + t183 * t150;
t107 = t182 * t131 + t183 * t152;
t103 = -t220 * t147 + t204;
t101 = t119 + t139;
t100 = t220 * t148 - t209;
t92 = t122 - t121;
t91 = -t187 * t127 - t226;
t90 = t191 * t127 - t228;
t88 = t191 * t120 - t252;
t87 = t187 * t120 + t250;
t83 = (-t123 * t190 + t125 * t186) * t161;
t82 = (-t123 * t186 - t125 * t190) * t161;
t75 = -t121 - t122;
t74 = t182 * t106 - t183 * t210;
t72 = -t125 * qJD(6) - t211;
t71 = t187 * t102 + t191 * t99;
t70 = -t191 * t102 + t187 * t99;
t69 = t190 * t109 - t238;
t68 = -t186 * t110 + t251;
t67 = t186 * t109 + t235;
t66 = t190 * t110 + t253;
t65 = -t188 * t103 + t192 * t91;
t64 = -t186 * t104 - t235;
t61 = -t188 * t100 + t192 * t88;
t55 = (qJD(6) + t161) * t125 + t211;
t54 = -t125 * t232 + t190 * t73;
t53 = t125 * t231 + t186 * t73;
t52 = t123 * t231 - t186 * t72;
t51 = t123 * t232 + t190 * t72;
t50 = t190 * t89 - t253;
t48 = -t188 * t113 + t192 * t71;
t43 = t182 * t62 - t183 * t96;
t42 = t182 * t65 - t183 * t90;
t41 = t182 * t61 - t183 * t87;
t39 = -t187 * t63 + t191 * t64;
t38 = t187 * t64 + t191 * t63;
t37 = t182 * t48 - t183 * t70;
t34 = t186 * t59 - t190 * t56;
t33 = -t186 * t249 - t190 * t55;
t31 = -t186 * t55 + t190 * t249;
t30 = -pkin(9) * t63 + t236;
t29 = -t187 * t49 + t191 * t50;
t28 = t187 * t50 + t191 * t49;
t27 = -pkin(9) * t49 + t239;
t24 = t188 * t249 + t192 * t39;
t23 = -pkin(5) * t249 + pkin(9) * t64 + t239;
t22 = t188 * t80 + t192 * t26;
t21 = t188 * t55 + t192 * t29;
t20 = -pkin(5) * t55 + pkin(9) * t50 - t236;
t16 = -t187 * t32 + t191 * t34;
t15 = t187 * t34 + t191 * t32;
t14 = t182 * t24 - t183 * t38;
t13 = t192 * t16 + t188 * t75;
t12 = t182 * t21 - t183 * t28;
t11 = t182 * t22 - t183 * t25;
t8 = -pkin(5) * t47 + pkin(9) * t10;
t7 = t182 * t13 - t183 * t15;
t6 = -pkin(9) * t32 - t9;
t5 = -pkin(5) * t75 + pkin(9) * t34 + t10;
t4 = t191 * t10 - t241;
t3 = t187 * t10 + t240;
t2 = t188 * t47 + t192 * t4;
t1 = t182 * t2 - t183 * t3;
t17 = [0, 0, 0, 0, 0, qJDD(1), t223, t207, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t213 + 0.2e1 * t233, 0, 0.2e1 * t178 + t202, pkin(1) * t140 + qJ(2) * (-t194 * pkin(1) + t198), 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t154 + t243 * t155 + t210, qJ(2) * t155 + t243 * t154 + t106, 0, qJ(2) * (t183 * t106 + t182 * t210) - t243 * t74, t206 * t188, t192 * t150 - t188 * t152, -t192 * (-t173 + t245) - t227, t205 * t192, -t225 - t188 * (t174 - t245), 0, qJ(2) * (t183 * t131 - t182 * t152) + t192 * t96 - pkin(3) * t152 - pkin(7) * t131 - t243 * t107, qJ(2) * (t183 * t132 - t182 * t150) - pkin(3) * t150 - pkin(7) * t132 - t188 * t96 - t243 * t108, qJ(2) * (t183 * t156 - t182 * t157) - pkin(3) * t157 - pkin(7) * t156 - t243 * t126 - t62, qJ(2) * (t182 * t96 + t183 * t62) + pkin(3) * t96 - pkin(7) * t62 - t243 * t43, t215 - t188 * (t191 * t119 + t148 * t230), t192 * t128 - t188 * (t191 * t100 - t187 * t101), t192 * t102 - t188 * (-t187 * t138 + t250), -t215 - t188 * (-t147 * t229 - t187 * t199), t192 * t99 - t188 * (t191 * t137 - t228), t192 * t146 - t188 * (t147 * t191 - t148 * t187) * t166, qJ(2) * (t182 * t87 + t183 * t61) - t192 * (-pkin(4) * t87 + t45) + pkin(3) * t87 - pkin(7) * t61 - t188 * (-pkin(8) * t87 + t237) - t243 * t41, qJ(2) * (t182 * t90 + t183 * t65) - t192 * (-pkin(4) * t90 + t46) + pkin(3) * t90 - pkin(7) * t65 - t188 * (-pkin(8) * t90 + t234) - t243 * t42, t188 * t25 + t196 * t70 + t212 * t48 - t243 * t37, -t243 * t11 + t196 * t25 + t212 * t22, t216 - t188 * (-t187 * t53 + t191 * t54), t192 * t92 - t188 * (-t187 * t31 + t191 * t33), t192 * t59 - t188 * (-t187 * t66 + t191 * t68), -t216 - t188 * (-t187 * t51 + t191 * t52), -t192 * t56 - t188 * (-t187 * t67 + t191 * t69), t192 * t143 - t188 * (-t187 * t82 + t191 * t83), qJ(2) * (t182 * t28 + t183 * t21) - t192 * (-pkin(4) * t28 - t203) + pkin(3) * t28 - pkin(7) * t21 - t188 * (-pkin(8) * t28 - t187 * t20 + t191 * t27) - t243 * t12, qJ(2) * (t182 * t38 + t183 * t24) - t192 * (-pkin(4) * t38 - t197) + pkin(3) * t38 - pkin(7) * t24 - t188 * (-pkin(8) * t38 - t187 * t23 + t191 * t30) - t243 * t14, qJ(2) * (t183 * t13 + t182 * t15) - t192 * (-pkin(4) * t15 - t242) + pkin(3) * t15 - pkin(7) * t13 - t188 * (-pkin(8) * t15 - t187 * t5 + t191 * t6) - t243 * t7, qJ(2) * (t182 * t3 + t183 * t2) - t192 * (-pkin(4) * t3 - t244) + pkin(3) * t3 - pkin(7) * t2 - t188 * (-pkin(8) * t3 - pkin(9) * t240 - t187 * t8) - t243 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t194, -t140, 0, 0, 0, 0, 0, 0, -t155, -t154, 0, t74, 0, 0, 0, 0, 0, 0, t107, t108, t126, t43, 0, 0, 0, 0, 0, 0, t41, t42, t37, t11, 0, 0, 0, 0, 0, 0, t12, t14, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, 0, 0, 0, 0, 0, 0, t192 * t158 + t188 * t163, -t188 * t159 + t192 * t162, 0, t188 * t94 - t192 * t93, 0, 0, 0, 0, 0, 0, t192 * t100 + t188 * t88, t192 * t103 + t188 * t91, t192 * t113 + t188 * t71, t188 * t26 - t192 * t80, 0, 0, 0, 0, 0, 0, t188 * t29 - t192 * t55, t188 * t39 - t192 * t249, t188 * t16 - t192 * t75, t188 * t4 - t192 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t173 - t174, -t217, t165, -t171, qJDD(4), -t93, -t94, 0, 0, t187 * t119 - t148 * t229, t187 * t100 + t191 * t101, t191 * t138 + t252, -t147 * t230 + t191 * t199, t187 * t137 + t226, (t147 * t187 + t148 * t191) * t166, pkin(4) * t100 + pkin(8) * t88 - t234, pkin(4) * t103 + pkin(8) * t91 + t237, pkin(4) * t113 + pkin(8) * t71 + t26, -pkin(4) * t80 + pkin(8) * t26, t187 * t54 + t191 * t53, t187 * t33 + t191 * t31, t187 * t68 + t191 * t66, t187 * t52 + t191 * t51, t187 * t69 + t191 * t67, t187 * t83 + t191 * t82, -pkin(4) * t55 + pkin(8) * t29 + t187 * t27 + t191 * t20, -pkin(4) * t249 + pkin(8) * t39 + t187 * t30 + t191 * t23, -pkin(4) * t75 + pkin(8) * t16 + t187 * t6 + t191 * t5, -pkin(4) * t47 + pkin(8) * t4 - pkin(9) * t241 + t191 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t128, t102, -t129, t99, t146, -t45, -t46, 0, 0, t95, t92, t59, -t95, -t56, t143, t203, t197, t242, t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t92, t59, -t95, -t56, t143, -t18, -t19, 0, 0;];
tauJ_reg  = t17;
