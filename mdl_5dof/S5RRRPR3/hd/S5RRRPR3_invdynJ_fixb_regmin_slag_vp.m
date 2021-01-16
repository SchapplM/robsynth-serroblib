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
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 22:59:36
% EndTime: 2021-01-15 22:59:46
% DurationCPUTime: 2.06s
% Computational Cost: add. (3101->305), mult. (4777->390), div. (0->0), fcn. (3381->16), ass. (0->196)
t186 = sin(qJ(2));
t240 = qJD(2) * t186;
t164 = pkin(1) * t240;
t190 = cos(qJ(2));
t271 = pkin(1) * t190;
t243 = -qJD(1) * t164 + qJDD(1) * t271;
t174 = qJDD(1) + qJDD(2);
t270 = pkin(2) * t174;
t106 = -t243 - t270;
t180 = qJ(1) + qJ(2);
t169 = sin(t180);
t170 = cos(t180);
t215 = g(2) * t170 + g(3) * t169;
t283 = t106 + t215;
t184 = sin(qJ(5));
t188 = cos(qJ(5));
t176 = qJD(1) + qJD(2);
t182 = cos(pkin(9));
t189 = cos(qJ(3));
t249 = t182 * t189;
t231 = t176 * t249;
t181 = sin(pkin(9));
t185 = sin(qJ(3));
t250 = t181 * t185;
t101 = t176 * t250 - t231;
t124 = t181 * t189 + t182 * t185;
t103 = t124 * t176;
t208 = t101 * t184 - t188 * t103;
t116 = t124 * qJD(3);
t247 = t189 * t174;
t248 = t185 * t174;
t213 = t181 * t248 - t182 * t247;
t55 = t176 * t116 + t213;
t238 = qJD(3) * t185;
t229 = t176 * t238;
t201 = t124 * t174 - t181 * t229;
t237 = qJD(3) * t189;
t228 = t176 * t237;
t56 = t182 * t228 + t201;
t194 = t208 * qJD(5) - t184 * t56 - t188 * t55;
t175 = qJD(3) + qJD(5);
t254 = t208 * t175;
t282 = t194 - t254;
t281 = qJ(4) + pkin(7);
t89 = t188 * t101;
t50 = -t103 * t184 - t89;
t256 = t175 * t50;
t236 = qJD(5) * t184;
t8 = -qJD(5) * t89 - t103 * t236 - t184 * t55 + t188 * t56;
t280 = t8 - t256;
t279 = t208 * t50;
t275 = g(2) * t169 - g(3) * t170;
t278 = t208 ^ 2 - t50 ^ 2;
t177 = qJ(3) + pkin(9);
t167 = qJ(5) + t177;
t154 = sin(t167);
t155 = cos(t167);
t265 = pkin(8) * t101;
t272 = pkin(1) * t186;
t234 = qJD(1) * t272;
t225 = t281 * t176 + t234;
t92 = t225 * t189;
t255 = t182 * t92;
t91 = t225 * t185;
t85 = qJD(3) * pkin(3) - t91;
t40 = t181 * t85 + t255;
t23 = t40 - t265;
t160 = pkin(3) * t189 + pkin(2);
t241 = qJD(1) * t190;
t232 = pkin(1) * t241;
t100 = -t160 * t176 + qJD(4) - t232;
t57 = t101 * pkin(4) + t100;
t277 = g(1) * t154 + t275 * t155 + t23 * t236 - t57 * t50;
t252 = t154 * t170;
t253 = t154 * t169;
t235 = qJDD(1) * t186;
t239 = qJD(2) * t190;
t107 = t174 * pkin(7) + (qJD(1) * t239 + t235) * pkin(1);
t203 = qJ(4) * t174 + qJD(4) * t176 + t107;
t207 = qJD(3) * t225;
t33 = qJDD(3) * pkin(3) - t203 * t185 - t189 * t207;
t36 = -t185 * t207 + t203 * t189;
t10 = -t181 * t36 + t182 * t33;
t4 = qJDD(3) * pkin(4) - pkin(8) * t56 + t10;
t11 = t181 * t33 + t182 * t36;
t5 = -pkin(8) * t55 + t11;
t276 = -g(1) * t155 + g(2) * t253 - g(3) * t252 - t184 * t5 + t188 * t4 + t57 * t208;
t168 = t189 * qJD(4);
t226 = qJD(3) * t281;
t113 = -t185 * t226 + t168;
t114 = -t185 * qJD(4) - t189 * t226;
t258 = -t113 * t181 + t182 * t114 + t124 * t232;
t123 = -t249 + t250;
t257 = t182 * t113 + t181 * t114 + t123 * t232;
t192 = qJD(3) ^ 2;
t274 = pkin(7) * t192 - t270;
t273 = qJD(5) - t175;
t269 = pkin(2) * t176;
t268 = pkin(3) * t181;
t267 = pkin(3) * t185;
t264 = pkin(8) * t103;
t117 = t123 * qJD(3);
t263 = pkin(8) * t117;
t262 = pkin(8) * t124;
t261 = g(1) * t189;
t159 = pkin(7) + t272;
t246 = -qJ(4) - t159;
t223 = qJD(3) * t246;
t233 = pkin(1) * t239;
t73 = t185 * t223 + t189 * t233 + t168;
t74 = (-qJD(4) - t233) * t185 + t189 * t223;
t35 = t181 * t74 + t182 * t73;
t79 = t181 * t92;
t42 = -t182 * t91 - t79;
t251 = t176 * t185;
t121 = t246 * t185;
t171 = t189 * qJ(4);
t122 = t159 * t189 + t171;
t66 = t181 * t121 + t182 * t122;
t145 = t281 * t185;
t146 = pkin(7) * t189 + t171;
t84 = -t181 * t145 + t182 * t146;
t245 = t169 * t160 - t170 * t281;
t178 = t185 ^ 2;
t242 = -t189 ^ 2 + t178;
t163 = pkin(3) * t238;
t230 = t176 * t240;
t88 = pkin(4) * t116 + t163;
t34 = -t181 * t73 + t182 * t74;
t39 = t182 * t85 - t79;
t41 = t181 * t91 - t255;
t65 = t182 * t121 - t122 * t181;
t83 = -t182 * t145 - t146 * t181;
t224 = t170 * t160 + t169 * t281;
t62 = pkin(3) * t229 - t160 * t174 + qJDD(4) - t243;
t22 = t55 * pkin(4) + t62;
t67 = t188 * t123 + t124 * t184;
t26 = -t67 * qJD(5) - t184 * t116 - t188 * t117;
t68 = -t123 * t184 + t124 * t188;
t222 = g(2) * t252 + g(3) * t253 + t22 * t68 + t57 * t26;
t165 = sin(t177);
t221 = -t100 * t117 + t62 * t124 + t215 * t165;
t220 = t176 * t234;
t134 = -t232 - t269;
t219 = t134 * t237 + t283 * t185;
t218 = t88 - t234;
t120 = t123 * pkin(8);
t61 = -t120 + t84;
t217 = qJD(5) * t61 - t258 - t263;
t112 = t116 * pkin(8);
t60 = t83 - t262;
t216 = -qJD(5) * t60 + t112 - t257;
t21 = qJD(3) * pkin(4) - t264 + t39;
t212 = -t184 * t21 - t188 * t23;
t43 = t65 - t262;
t44 = -t120 + t66;
t211 = -t184 * t44 + t188 * t43;
t210 = t184 * t43 + t188 * t44;
t209 = -t10 * t124 - t11 * t123 - t40 * t116 + t39 * t117 - t275;
t96 = pkin(4) * t123 - t160;
t156 = pkin(3) * t182 + pkin(4);
t206 = t156 * t184 + t188 * t268;
t205 = t156 * t188 - t184 * t268;
t202 = -t134 * t176 - t107 + t275;
t27 = t68 * qJD(5) + t188 * t116 - t184 * t117;
t200 = -t215 * t155 + t22 * t67 + t57 * t27;
t166 = cos(t177);
t199 = t100 * t116 + t62 * t123 - t215 * t166;
t198 = -t215 + t220;
t161 = -pkin(2) - t271;
t197 = pkin(1) * t230 + t159 * t192 + t161 * t174;
t195 = -pkin(7) * qJDD(3) + (t232 - t269) * qJD(3);
t193 = -qJDD(3) * t159 + (t161 * t176 - t233) * qJD(3);
t191 = cos(qJ(1));
t187 = sin(qJ(1));
t173 = qJDD(3) + qJDD(5);
t172 = t176 ^ 2;
t143 = -t160 - t271;
t141 = qJDD(3) * t189 - t185 * t192;
t140 = qJDD(3) * t185 + t189 * t192;
t132 = t164 + t163;
t118 = t134 * t238;
t108 = t174 * t178 + 0.2e1 * t185 * t228;
t87 = t96 - t271;
t78 = -0.2e1 * t242 * t176 * qJD(3) + 0.2e1 * t185 * t247;
t77 = t164 + t88;
t72 = pkin(3) * t251 + pkin(4) * t103;
t25 = t42 - t264;
t24 = t41 + t265;
t19 = -t112 + t35;
t18 = t34 + t263;
t13 = -t173 * t67 - t175 * t27;
t12 = t173 * t68 + t175 * t26;
t2 = -t208 * t26 + t68 * t8;
t1 = t194 * t68 + t208 * t27 + t26 * t50 - t67 * t8;
t3 = [qJDD(1), -g(2) * t191 - g(3) * t187, g(2) * t187 - g(3) * t191, t174, (t174 * t190 - t230) * pkin(1) - t215 + t243, ((-qJDD(1) - t174) * t186 + (-qJD(1) - t176) * t239) * pkin(1) + t275, t108, t78, t140, t141, 0, t118 + t193 * t185 + (-t197 - t283) * t189, t197 * t185 + t193 * t189 + t219, qJD(3) * t34 + qJDD(3) * t65 + t101 * t132 + t143 * t55 + t199, -qJD(3) * t35 - qJDD(3) * t66 + t103 * t132 + t143 * t56 + t221, -t101 * t35 - t103 * t34 - t55 * t66 - t56 * t65 + t209, t11 * t66 + t40 * t35 + t10 * t65 + t39 * t34 + t62 * t143 + t100 * t132 - g(2) * (pkin(1) * t191 + t224) - g(3) * (pkin(1) * t187 + t245), t2, t1, t12, t13, 0, -t77 * t50 - t87 * t194 + (-qJD(5) * t210 + t188 * t18 - t184 * t19) * t175 + t211 * t173 + t200, -t77 * t208 + t87 * t8 - (qJD(5) * t211 + t184 * t18 + t188 * t19) * t175 - t210 * t173 + t222; 0, 0, 0, t174, t198 + t243, (-t235 + (-qJD(2) + t176) * t241) * pkin(1) + t275, t108, t78, t140, t141, 0, t118 + t195 * t185 + (-t106 + t198 - t274) * t189, t195 * t189 + (-t220 + t274) * t185 + t219, -t101 * t234 + qJDD(3) * t83 - t160 * t55 + (t101 * t267 + t258) * qJD(3) + t199, -t103 * t234 - qJDD(3) * t84 - t160 * t56 + (t103 * t267 - t257) * qJD(3) + t221, -t257 * t101 - t258 * t103 - t55 * t84 - t56 * t83 + t209, t11 * t84 + t10 * t83 - t62 * t160 - g(2) * t224 - g(3) * t245 + t257 * t40 + t258 * t39 + (t163 - t234) * t100, t2, t1, t12, t13, 0, -t96 * t194 + (-t184 * t61 + t188 * t60) * t173 - t218 * t50 + (t184 * t216 - t188 * t217) * t175 + t200, t96 * t8 - (t184 * t60 + t188 * t61) * t173 - t218 * t208 + (t184 * t217 + t188 * t216) * t175 + t222; 0, 0, 0, 0, 0, 0, -t185 * t172 * t189, t242 * t172, t248, t247, qJDD(3), t202 * t185 - t261, g(1) * t185 + t189 * t202, -g(1) * t166 - qJD(3) * t41 - t100 * t103 + t275 * t165 + (qJDD(3) * t182 - t101 * t251) * pkin(3) + t10, g(1) * t165 + qJD(3) * t42 + t100 * t101 + t275 * t166 + (-qJDD(3) * t181 - t103 * t251) * pkin(3) - t11, (t40 + t41) * t103 + (-t39 + t42) * t101 + (-t181 * t55 - t182 * t56) * pkin(3), -t39 * t41 - t40 * t42 + (-t261 + t10 * t182 + t11 * t181 + (-t100 * t176 + t275) * t185) * pkin(3), t279, t278, t280, t282, t173, t205 * t173 + t72 * t50 - (-t184 * t25 + t188 * t24) * t175 + (-t175 * t206 + t212) * qJD(5) + t276, -t206 * t173 - t188 * t5 - t184 * t4 + t72 * t208 + (t184 * t24 + t188 * t25) * t175 + (-t175 * t205 - t188 * t21) * qJD(5) + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103 * qJD(3) + t213, (-t101 + t231) * qJD(3) + t201, -t101 ^ 2 - t103 ^ 2, t101 * t40 + t103 * t39 + t215 + t62, 0, 0, 0, 0, 0, -t194 - t254, t8 + t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t279, t278, t280, t282, t173, t212 * t273 + t276, (-t175 * t23 - t4) * t184 + (-t21 * t273 - t5) * t188 + t277;];
tau_reg = t3;
