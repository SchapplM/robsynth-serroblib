% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR13_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:53
% EndTime: 2019-12-31 18:33:03
% DurationCPUTime: 4.43s
% Computational Cost: add. (6850->326), mult. (16820->429), div. (0->0), fcn. (11914->8), ass. (0->189)
t161 = sin(pkin(8));
t162 = cos(pkin(8));
t169 = qJD(3) ^ 2;
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t211 = t161 * t168;
t184 = t162 * t165 + t211;
t145 = t184 * qJD(1);
t239 = t145 ^ 2;
t124 = t239 + t169;
t212 = t161 * t165;
t143 = (-t162 * t168 + t212) * qJD(1);
t214 = t145 * t143;
t258 = qJDD(3) + t214;
t272 = t165 * t258;
t71 = t124 * t168 + t272;
t270 = t168 * t258;
t73 = -t124 * t165 + t270;
t283 = qJ(2) * (t161 * t71 - t162 * t73);
t282 = pkin(6) * t71;
t281 = pkin(6) * t73;
t125 = t239 - t169;
t259 = qJDD(3) - t214;
t269 = t168 * t259;
t271 = t165 * t259;
t279 = t161 * (t125 * t165 + t269) - t162 * (t125 * t168 - t271);
t240 = t143 ^ 2;
t120 = t240 - t169;
t276 = t161 * (-t120 * t168 + t272) - t162 * (t120 * t165 + t270);
t97 = -t169 - t240;
t59 = t165 * t97 + t269;
t62 = -t168 * t97 + t271;
t273 = qJ(2) * (t161 * t59 + t162 * t62);
t268 = pkin(6) * t59;
t267 = pkin(6) * t62;
t131 = t143 * qJD(3);
t142 = t184 * qJDD(1);
t108 = t142 - t131;
t79 = -t131 + t108;
t248 = qJ(4) * t79;
t208 = t145 * qJD(3);
t170 = qJD(1) ^ 2;
t166 = sin(qJ(1));
t236 = cos(qJ(1));
t195 = g(1) * t166 - t236 * g(2);
t187 = -qJDD(2) + t195;
t157 = t161 ^ 2;
t158 = t162 ^ 2;
t210 = t157 + t158;
t233 = t162 * pkin(2);
t98 = t170 * (pkin(6) * t210 + qJ(2)) + (pkin(1) + t233) * qJDD(1) + t187;
t266 = pkin(3) * t208 - 0.2e1 * qJD(4) * t145 - t98;
t182 = g(1) * t236 + t166 * g(2);
t251 = -t170 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t182;
t244 = -t240 - t239;
t256 = pkin(1) * t244;
t255 = pkin(2) * t244;
t234 = g(3) * t162;
t85 = -t234 + (-pkin(6) * qJDD(1) + t170 * t233 - t251) * t161;
t188 = -g(3) * t161 + t251 * t162;
t203 = t162 * qJDD(1);
t213 = t158 * t170;
t88 = -pkin(2) * t213 + pkin(6) * t203 + t188;
t53 = t165 * t88 - t168 * t85;
t94 = pkin(3) * t143 - qJ(4) * t145;
t39 = -qJDD(3) * pkin(3) - t169 * qJ(4) + t145 * t94 + qJDD(4) + t53;
t254 = -pkin(7) * t259 + t39;
t164 = sin(qJ(5));
t167 = cos(qJ(5));
t116 = qJD(3) * t164 - t167 * t143;
t118 = qJD(3) * t167 + t143 * t164;
t83 = t118 * t116;
t95 = qJDD(5) + t108;
t249 = -t83 + t95;
t253 = t164 * t249;
t252 = t167 * t249;
t246 = t131 + t108;
t245 = t239 - t240;
t219 = qJDD(1) * pkin(1);
t220 = qJ(2) * t170;
t136 = t187 + t219 + t220;
t243 = qJ(2) * t213 + t157 * t220 - t136 - t219;
t114 = t116 ^ 2;
t115 = t118 ^ 2;
t135 = qJD(5) + t145;
t133 = t135 ^ 2;
t238 = 0.2e1 * qJD(4);
t237 = pkin(3) + pkin(7);
t235 = pkin(3) * t168;
t54 = t165 * t85 + t168 * t88;
t33 = t165 * t54 - t168 * t53;
t232 = t161 * t33;
t204 = t161 * qJDD(1);
t186 = t165 * t204 - t168 * t203;
t106 = t186 + t208;
t119 = pkin(4) * t145 - qJD(3) * pkin(7);
t180 = -t169 * pkin(3) - t143 * t94 + t54;
t205 = qJDD(3) * qJ(4);
t25 = t205 - t240 * pkin(7) - t106 * pkin(4) + (t238 + t119) * qJD(3) + t180;
t231 = t164 * t25;
t57 = t83 + t95;
t230 = t164 * t57;
t228 = t165 * t98;
t173 = -t248 + t266;
t23 = -pkin(4) * t240 + t106 * t237 - t119 * t145 + t173;
t226 = t167 * t23;
t225 = t167 * t25;
t224 = t167 * t57;
t222 = t168 * t98;
t216 = t135 * t164;
t215 = t135 * t167;
t207 = qJD(5) + t135;
t199 = t165 * t83;
t198 = t168 * t83;
t196 = qJ(4) * t165 + pkin(2);
t175 = pkin(4) * t246 + t254;
t12 = t164 * t23 - t167 * t175;
t34 = t165 * t53 + t168 * t54;
t193 = t161 * (t251 * t161 + t234) + t162 * t188;
t190 = qJDD(3) * t164 - t167 * t106;
t13 = t164 * t175 + t226;
t6 = -t12 * t167 + t13 * t164;
t7 = t12 * t164 + t13 * t167;
t183 = qJDD(3) * t167 + t106 * t164;
t181 = (-qJD(5) + t135) * t118 - t190;
t69 = -qJD(5) * t116 + t183;
t179 = qJD(3) * t238 + t180;
t178 = t161 * (t168 * t108 - t165 * t208) + t162 * (t165 * t108 + t168 * t208);
t38 = t179 + t205;
t177 = t161 * (t106 * t165 + t131 * t168) + t162 * (-t168 * t106 + t131 * t165);
t176 = (t161 * (-t143 * t168 + t145 * t165) + t162 * (-t143 * t165 - t145 * t168)) * qJD(3);
t174 = -pkin(3) * t106 - t266;
t154 = t158 * qJDD(1);
t153 = t157 * qJDD(1);
t148 = t210 * t170;
t107 = t142 - 0.2e1 * t131;
t105 = t186 + 0.2e1 * t208;
t89 = t135 * t116;
t87 = -t115 + t133;
t86 = t114 - t133;
t80 = t115 - t114;
t77 = t106 + t208;
t76 = t106 - t208;
t75 = -t115 - t133;
t70 = -t133 - t114;
t68 = -qJD(5) * t118 - t190;
t67 = -t114 - t115;
t66 = t165 * t246 - t168 * t186;
t65 = t142 * t165 - t168 * t76;
t64 = -t165 * t186 - t168 * t246;
t63 = -t142 * t168 - t165 * t76;
t55 = (t116 * t164 + t118 * t167) * t135;
t51 = -t116 * t207 + t183;
t50 = t69 + t89;
t49 = t69 - t89;
t46 = t118 * t207 + t190;
t45 = -t118 * t215 - t164 * t69;
t44 = -t116 * t216 - t167 * t68;
t43 = -t164 * t86 - t224;
t42 = -t167 * t87 - t253;
t41 = -t164 * t75 - t224;
t40 = t167 * t75 - t230;
t37 = t167 * t70 - t253;
t36 = t164 * t70 + t252;
t35 = t174 + t248;
t32 = -qJ(4) * t244 + t39;
t31 = -pkin(3) * t244 + t38;
t30 = (t106 + t77) * pkin(3) + t173;
t29 = t174 + 0.2e1 * t248;
t28 = t164 * t50 + t167 * t181;
t27 = t164 * t181 - t167 * t50;
t26 = t164 * t46 - t167 * t49;
t22 = t165 * t40 + t168 * t51;
t21 = t165 * t51 - t168 * t40;
t20 = t165 * t36 + t168 * t46;
t19 = t165 * t46 - t168 * t36;
t18 = t165 * t27 + t168 * t67;
t17 = t165 * t67 - t168 * t27;
t14 = pkin(4) * t27 - qJ(4) * t28;
t11 = pkin(4) * t51 - t237 * t41 - t231;
t10 = pkin(4) * t46 - t237 * t37 + t225;
t9 = -t226 - t164 * t254 - qJ(4) * t41 + (-t164 * t246 + t40) * pkin(4);
t8 = pkin(4) * t36 - qJ(4) * t37 - t12;
t5 = t165 * t6 + t168 * t25;
t4 = t165 * t25 - t168 * t6;
t3 = pkin(4) * t67 - t237 * t28 - t7;
t2 = pkin(4) * t6 - qJ(4) * t7;
t1 = pkin(4) * t25 - t237 * t7;
t15 = [0, 0, 0, 0, 0, qJDD(1), t195, t182, 0, 0, t153, 0.2e1 * t161 * t203, 0, t154, 0, 0, -t243 * t162, t243 * t161, pkin(1) * t148 + qJ(2) * (t154 + t153) + t193, pkin(1) * t136 + qJ(2) * t193, t178, t161 * (-t105 * t168 - t107 * t165) + t162 * (-t105 * t165 + t107 * t168), t279, t177, -t276, t176, t161 * (-t228 - t268) + t162 * (-pkin(2) * t105 + t222 - t267) - pkin(1) * t105 - t273, t161 * (-t222 + t282) + t162 * (-pkin(2) * t107 - t228 - t281) - pkin(1) * t107 + t283, t161 * (-pkin(6) * t63 - t33) + t162 * (pkin(6) * t65 - t255 + t34) - t256 + qJ(2) * (-t161 * t63 + t162 * t65), -pkin(6) * t232 + t162 * (pkin(2) * t98 + pkin(6) * t34) + pkin(1) * t98 + qJ(2) * (t162 * t34 - t232), t176, -t279, t276, t178, t161 * (-t165 * t79 - t168 * t77) + t162 * (-t165 * t77 + t168 * t79), t177, t161 * (-pkin(6) * t64 - t165 * t31 + t168 * t32) + t162 * (pkin(6) * t66 + t165 * t32 + t168 * t31 - t255) - t256 + qJ(2) * (-t161 * t64 + t162 * t66), t161 * (-t165 * t30 + t268) + t162 * (t168 * t30 + t267) + t273 + (qJ(4) * t211 + t162 * t196 + pkin(1)) * t77, t161 * (t168 * t29 - t282) + t162 * (t165 * t29 + t281) - t283 + (-pkin(3) * t212 + t162 * (pkin(2) + t235) + pkin(1)) * t79, (t161 * (-pkin(3) * t165 + qJ(4) * t168) + t162 * (t196 + t235) + pkin(1)) * t35 + (qJ(2) + pkin(6)) * (-(t165 * t38 - t168 * t39) * t161 + (t165 * t39 + t168 * t38) * t162), t161 * (-t165 * t45 + t198) + t162 * (t168 * t45 + t199), t161 * (-t165 * t26 + t168 * t80) + t162 * (t165 * t80 + t168 * t26), t161 * (-t165 * t42 + t168 * t50) + t162 * (t165 * t50 + t168 * t42), t161 * (-t165 * t44 - t198) + t162 * (t168 * t44 - t199), t161 * (-t165 * t43 + t168 * t181) + t162 * (t165 * t181 + t168 * t43), t161 * (-t165 * t55 + t168 * t95) + t162 * (t165 * t95 + t168 * t55), t161 * (-pkin(6) * t19 - t10 * t165 + t168 * t8) + t162 * (-pkin(2) * t37 + pkin(6) * t20 + t10 * t168 + t165 * t8) - pkin(1) * t37 + qJ(2) * (-t161 * t19 + t162 * t20), t161 * (-pkin(6) * t21 - t11 * t165 + t168 * t9) + t162 * (-pkin(2) * t41 + pkin(6) * t22 + t11 * t168 + t165 * t9) - pkin(1) * t41 + qJ(2) * (-t161 * t21 + t162 * t22), t161 * (-pkin(6) * t17 + t14 * t168 - t165 * t3) + t162 * (-pkin(2) * t28 + pkin(6) * t18 + t14 * t165 + t168 * t3) - pkin(1) * t28 + qJ(2) * (-t161 * t17 + t162 * t18), t161 * (-pkin(6) * t4 - t1 * t165 + t168 * t2) + t162 * (-pkin(2) * t7 + pkin(6) * t5 + t1 * t168 + t165 * t2) - pkin(1) * t7 + qJ(2) * (-t161 * t4 + t162 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203, t204, -t148, -t136, 0, 0, 0, 0, 0, 0, t105, t107, t244, -t98, 0, 0, 0, 0, 0, 0, t244, -t77, -t79, -t35, 0, 0, 0, 0, 0, 0, t37, t41, t28, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, t245, t142, -t214, -t186, qJDD(3), -t53, -t54, 0, 0, qJDD(3), -t246, t76, t214, t245, -t214, -pkin(3) * t246 - qJ(4) * t186, -pkin(3) * t259 - qJ(4) * t97 + t39, pkin(3) * t124 + (qJDD(3) + t258) * qJ(4) + t179, -pkin(3) * t39 + qJ(4) * t38, -t118 * t216 + t167 * t69, -t164 * t49 - t167 * t46, -t164 * t87 + t252, t116 * t215 - t164 * t68, t167 * t86 - t230, (-t116 * t167 + t118 * t164) * t135, qJ(4) * t46 - t237 * t36 + t231, qJ(4) * t51 - t237 * t40 + t225, qJ(4) * t67 - t237 * t27 - t6, qJ(4) * t25 - t237 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, t259, -t124, t39, 0, 0, 0, 0, 0, 0, t36, t40, t27, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t80, t50, -t83, t181, t95, -t12, -t13, 0, 0;];
tauJ_reg = t15;
