% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:15
% EndTime: 2022-01-20 11:43:21
% DurationCPUTime: 2.08s
% Computational Cost: add. (3101->307), mult. (4777->395), div. (0->0), fcn. (3381->16), ass. (0->199)
t185 = sin(qJ(5));
t189 = cos(qJ(5));
t177 = qJD(1) + qJD(2);
t183 = cos(pkin(9));
t190 = cos(qJ(3));
t252 = t183 * t190;
t233 = t177 * t252;
t182 = sin(pkin(9));
t186 = sin(qJ(3));
t253 = t182 * t186;
t101 = t177 * t253 - t233;
t124 = t182 * t190 + t183 * t186;
t103 = t124 * t177;
t212 = t101 * t185 - t103 * t189;
t116 = t124 * qJD(3);
t175 = qJDD(1) + qJDD(2);
t250 = t190 * t175;
t251 = t186 * t175;
t217 = t182 * t251 - t183 * t250;
t55 = t116 * t177 + t217;
t241 = qJD(3) * t186;
t231 = t177 * t241;
t199 = t124 * t175 - t182 * t231;
t240 = qJD(3) * t190;
t230 = t177 * t240;
t56 = t183 * t230 + t199;
t195 = qJD(5) * t212 - t185 * t56 - t189 * t55;
t176 = qJD(3) + qJD(5);
t259 = t212 * t176;
t286 = t195 - t259;
t285 = qJ(4) + pkin(7);
t89 = t189 * t101;
t50 = -t103 * t185 - t89;
t261 = t176 * t50;
t239 = qJD(5) * t185;
t8 = -qJD(5) * t89 - t103 * t239 - t185 * t55 + t189 * t56;
t284 = t8 - t261;
t283 = t212 * t50;
t187 = sin(qJ(2));
t243 = qJD(2) * t187;
t165 = pkin(1) * t243;
t191 = cos(qJ(2));
t274 = pkin(1) * t191;
t246 = -qJD(1) * t165 + qJDD(1) * t274;
t181 = qJ(1) + qJ(2);
t171 = cos(t181);
t265 = g(2) * t171;
t273 = pkin(2) * t175;
t282 = -t246 - t273 + t265;
t170 = sin(t181);
t158 = g(1) * t170;
t277 = t265 - t158;
t281 = t212 ^ 2 - t50 ^ 2;
t178 = qJ(3) + pkin(9);
t168 = qJ(5) + t178;
t154 = sin(t168);
t269 = pkin(8) * t101;
t275 = pkin(1) * t187;
t237 = qJD(1) * t275;
t226 = t177 * t285 + t237;
t92 = t226 * t190;
t260 = t183 * t92;
t91 = t226 * t186;
t85 = qJD(3) * pkin(3) - t91;
t40 = t182 * t85 + t260;
t23 = t40 - t269;
t155 = cos(t168);
t255 = t155 * t171;
t256 = t155 * t170;
t161 = pkin(3) * t190 + pkin(2);
t244 = qJD(1) * t191;
t235 = pkin(1) * t244;
t100 = -t161 * t177 + qJD(4) - t235;
t57 = t101 * pkin(4) + t100;
t280 = g(1) * t255 + g(2) * t256 + g(3) * t154 + t23 * t239 - t57 * t50;
t257 = t154 * t171;
t258 = t154 * t170;
t238 = qJDD(1) * t187;
t242 = qJD(2) * t191;
t107 = t175 * pkin(7) + (qJD(1) * t242 + t238) * pkin(1);
t202 = qJ(4) * t175 + qJD(4) * t177 + t107;
t210 = qJD(3) * t226;
t33 = qJDD(3) * pkin(3) - t186 * t202 - t190 * t210;
t36 = -t186 * t210 + t190 * t202;
t10 = -t182 * t36 + t183 * t33;
t4 = qJDD(3) * pkin(4) - pkin(8) * t56 + t10;
t11 = t182 * t33 + t183 * t36;
t5 = -pkin(8) * t55 + t11;
t279 = g(1) * t257 + g(2) * t258 - g(3) * t155 - t185 * t5 + t189 * t4 + t57 * t212;
t169 = t190 * qJD(4);
t227 = qJD(3) * t285;
t113 = -t186 * t227 + t169;
t114 = -t186 * qJD(4) - t190 * t227;
t263 = -t113 * t182 + t114 * t183 + t124 * t235;
t123 = -t252 + t253;
t262 = t113 * t183 + t114 * t182 + t123 * t235;
t278 = g(1) * t171 + g(2) * t170;
t276 = qJD(5) - t176;
t272 = pkin(2) * t177;
t271 = pkin(3) * t182;
t270 = pkin(3) * t186;
t268 = pkin(8) * t103;
t117 = t123 * qJD(3);
t267 = pkin(8) * t117;
t266 = pkin(8) * t124;
t264 = g(3) * t190;
t160 = pkin(7) + t275;
t249 = -qJ(4) - t160;
t224 = qJD(3) * t249;
t236 = pkin(1) * t242;
t73 = t186 * t224 + t190 * t236 + t169;
t74 = (-qJD(4) - t236) * t186 + t190 * t224;
t35 = t182 * t74 + t183 * t73;
t79 = t182 * t92;
t42 = -t183 * t91 - t79;
t254 = t177 * t186;
t121 = t249 * t186;
t172 = t190 * qJ(4);
t122 = t160 * t190 + t172;
t66 = t121 * t182 + t122 * t183;
t136 = -t235 - t272;
t248 = t136 * t241 + t158 * t190;
t146 = t285 * t186;
t147 = pkin(7) * t190 + t172;
t84 = -t146 * t182 + t147 * t183;
t179 = t186 ^ 2;
t245 = -t190 ^ 2 + t179;
t164 = pkin(3) * t241;
t234 = t136 * t240 + t186 * t282;
t232 = t177 * t243;
t88 = pkin(4) * t116 + t164;
t34 = -t182 * t73 + t183 * t74;
t39 = t183 * t85 - t79;
t41 = t182 * t91 - t260;
t65 = t121 * t183 - t122 * t182;
t83 = -t146 * t183 - t147 * t182;
t225 = t161 * t171 + t170 * t285;
t223 = t177 * t237;
t221 = t88 - t237;
t120 = t123 * pkin(8);
t61 = -t120 + t84;
t220 = qJD(5) * t61 - t263 - t267;
t112 = t116 * pkin(8);
t60 = t83 - t266;
t219 = -qJD(5) * t60 + t112 - t262;
t21 = qJD(3) * pkin(4) - t268 + t39;
t216 = -t185 * t21 - t189 * t23;
t43 = t65 - t266;
t44 = -t120 + t66;
t215 = -t185 * t44 + t189 * t43;
t214 = t185 * t43 + t189 * t44;
t213 = -t10 * t124 - t11 * t123 - t116 * t40 + t117 * t39 - t278;
t67 = t123 * t189 + t124 * t185;
t68 = -t123 * t185 + t124 * t189;
t211 = -t161 * t170 + t171 * t285;
t96 = pkin(4) * t123 - t161;
t209 = -t246 + t277;
t156 = pkin(3) * t183 + pkin(4);
t208 = t156 * t185 + t189 * t271;
t207 = t156 * t189 - t185 * t271;
t62 = pkin(3) * t231 - t161 * t175 + qJDD(4) - t246;
t22 = t55 * pkin(4) + t62;
t26 = -qJD(5) * t67 - t185 * t116 - t189 * t117;
t206 = -g(1) * t258 + g(2) * t257 + t22 * t68 + t26 * t57;
t166 = sin(t178);
t205 = -t100 * t117 + t62 * t124 + t166 * t277;
t27 = qJD(5) * t68 + t116 * t189 - t185 * t117;
t204 = g(1) * t256 - g(2) * t255 + t22 * t67 + t27 * t57;
t167 = cos(t178);
t203 = t100 * t116 + t62 * t123 - t167 * t277;
t200 = -t136 * t177 - t107 + t278;
t193 = qJD(3) ^ 2;
t198 = pkin(7) * t193 - t223 - t273;
t162 = -pkin(2) - t274;
t197 = pkin(1) * t232 + t160 * t193 + t162 * t175;
t196 = -pkin(7) * qJDD(3) + (t235 - t272) * qJD(3);
t194 = -qJDD(3) * t160 + (t162 * t177 - t236) * qJD(3);
t192 = cos(qJ(1));
t188 = sin(qJ(1));
t174 = qJDD(3) + qJDD(5);
t173 = t177 ^ 2;
t144 = -t161 - t274;
t142 = qJDD(3) * t190 - t186 * t193;
t141 = qJDD(3) * t186 + t190 * t193;
t134 = t165 + t164;
t108 = t175 * t179 + 0.2e1 * t186 * t230;
t87 = t96 - t274;
t78 = -0.2e1 * qJD(3) * t177 * t245 + 0.2e1 * t186 * t250;
t77 = t165 + t88;
t72 = pkin(3) * t254 + pkin(4) * t103;
t25 = t42 - t268;
t24 = t41 + t269;
t19 = -t112 + t35;
t18 = t34 + t267;
t13 = -t174 * t67 - t176 * t27;
t12 = t174 * t68 + t176 * t26;
t2 = -t212 * t26 + t68 * t8;
t1 = t195 * t68 + t212 * t27 + t26 * t50 - t67 * t8;
t3 = [qJDD(1), g(1) * t188 - g(2) * t192, g(1) * t192 + g(2) * t188, t175, (t175 * t191 - t232) * pkin(1) - t209, ((-qJDD(1) - t175) * t187 + (-qJD(1) - t177) * t242) * pkin(1) + t278, t108, t78, t141, t142, 0, t194 * t186 + (-t197 - t282) * t190 + t248, t194 * t190 + (t197 - t158) * t186 + t234, qJD(3) * t34 + qJDD(3) * t65 + t101 * t134 + t144 * t55 + t203, -qJD(3) * t35 - qJDD(3) * t66 + t103 * t134 + t144 * t56 + t205, -t101 * t35 - t103 * t34 - t55 * t66 - t56 * t65 + t213, t11 * t66 + t40 * t35 + t10 * t65 + t39 * t34 + t62 * t144 + t100 * t134 - g(1) * (-pkin(1) * t188 + t211) - g(2) * (pkin(1) * t192 + t225), t2, t1, t12, t13, 0, -t77 * t50 - t87 * t195 + (-qJD(5) * t214 + t189 * t18 - t185 * t19) * t176 + t215 * t174 + t204, -t77 * t212 + t87 * t8 - (qJD(5) * t215 + t185 * t18 + t189 * t19) * t176 - t214 * t174 + t206; 0, 0, 0, t175, -t209 + t223, (-t238 + (-qJD(2) + t177) * t244) * pkin(1) + t278, t108, t78, t141, t142, 0, t196 * t186 + (-t198 - t282) * t190 + t248, t196 * t190 + (t198 - t158) * t186 + t234, -t101 * t237 + qJDD(3) * t83 - t161 * t55 + (t101 * t270 + t263) * qJD(3) + t203, -t103 * t237 - qJDD(3) * t84 - t161 * t56 + (t103 * t270 - t262) * qJD(3) + t205, -t101 * t262 - t103 * t263 - t55 * t84 - t56 * t83 + t213, t11 * t84 + t10 * t83 - t62 * t161 - g(1) * t211 - g(2) * t225 + t262 * t40 + t263 * t39 + (t164 - t237) * t100, t2, t1, t12, t13, 0, -t96 * t195 + (-t185 * t61 + t189 * t60) * t174 - t221 * t50 + (t185 * t219 - t189 * t220) * t176 + t204, t96 * t8 - (t185 * t60 + t189 * t61) * t174 - t221 * t212 + (t185 * t220 + t189 * t219) * t176 + t206; 0, 0, 0, 0, 0, 0, -t186 * t173 * t190, t245 * t173, t251, t250, qJDD(3), t186 * t200 - t264, g(3) * t186 + t190 * t200, -g(3) * t167 - qJD(3) * t41 - t100 * t103 + t278 * t166 + (qJDD(3) * t183 - t101 * t254) * pkin(3) + t10, g(3) * t166 + t42 * qJD(3) + t100 * t101 + t278 * t167 + (-qJDD(3) * t182 - t103 * t254) * pkin(3) - t11, (t40 + t41) * t103 + (-t39 + t42) * t101 + (-t182 * t55 - t183 * t56) * pkin(3), -t39 * t41 - t40 * t42 + (-t264 + t10 * t183 + t11 * t182 + (-t100 * t177 + t278) * t186) * pkin(3), t283, t281, t284, t286, t174, t207 * t174 + t72 * t50 - (-t185 * t25 + t189 * t24) * t176 + (-t176 * t208 + t216) * qJD(5) + t279, -t208 * t174 - t189 * t5 - t185 * t4 + t72 * t212 + (t185 * t24 + t189 * t25) * t176 + (-t176 * t207 - t189 * t21) * qJD(5) + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103 * qJD(3) + t217, (-t101 + t233) * qJD(3) + t199, -t101 ^ 2 - t103 ^ 2, t101 * t40 + t103 * t39 + t277 + t62, 0, 0, 0, 0, 0, -t195 - t259, t8 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, t281, t284, t286, t174, t216 * t276 + t279, (-t176 * t23 - t4) * t185 + (-t21 * t276 - t5) * t189 + t280;];
tau_reg = t3;
