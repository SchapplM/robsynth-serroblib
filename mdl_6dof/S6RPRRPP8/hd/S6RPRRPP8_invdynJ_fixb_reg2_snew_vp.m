% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:52:54
% EndTime: 2019-05-05 21:53:06
% DurationCPUTime: 4.02s
% Computational Cost: add. (6291->308), mult. (12162->332), div. (0->0), fcn. (7364->6), ass. (0->202)
t271 = pkin(7) + pkin(1);
t185 = sin(qJ(3));
t188 = cos(qJ(3));
t184 = sin(qJ(4));
t187 = cos(qJ(4));
t235 = qJD(1) * t188;
t160 = qJD(3) * t184 + t187 * t235;
t177 = t188 * qJDD(1);
t232 = qJD(1) * qJD(3);
t224 = t185 * t232;
t164 = t177 - t224;
t106 = qJD(4) * t160 - t187 * qJDD(3) + t164 * t184;
t175 = qJD(1) * t185 + qJD(4);
t246 = t175 * t160;
t283 = t106 + t246;
t223 = t188 * t232;
t230 = t185 * qJDD(1);
t163 = -t223 - t230;
t157 = qJDD(4) - t163;
t158 = -t187 * qJD(3) + t184 * t235;
t247 = t160 * t158;
t97 = t247 - t157;
t260 = t184 * t97;
t174 = t175 ^ 2;
t274 = t158 ^ 2;
t282 = -t174 - t274;
t301 = t282 * t187 + t260;
t313 = t185 * t301 - t188 * t283;
t87 = t187 * t97;
t48 = -t282 * t184 + t87;
t323 = -qJ(2) * t48 - t271 * t313;
t273 = t160 ^ 2;
t280 = -t273 - t174;
t281 = t247 + t157;
t291 = t184 * t281;
t316 = -t187 * t280 + t291;
t289 = t187 * t281;
t317 = t184 * t280 + t289;
t209 = t184 * qJDD(3) + t187 * t164;
t80 = (qJD(4) + t175) * t158 - t209;
t325 = t185 * t317 - t188 * t80;
t332 = qJ(2) * t316 - t271 * t325;
t308 = pkin(3) * t48;
t307 = pkin(8) * t48;
t335 = pkin(8) * t301;
t136 = -t273 + t174;
t314 = t188 * (t136 * t184 + t87);
t75 = (-qJD(4) + t175) * t158 + t209;
t333 = t185 * t75 - t314;
t331 = pkin(3) * t316;
t330 = pkin(8) * t316;
t329 = pkin(8) * t317;
t234 = qJD(5) * t175;
t327 = -0.2e1 * t234 - t331;
t326 = pkin(3) * t80 - t329;
t135 = t274 - t174;
t284 = t106 - t246;
t324 = -t185 * t284 + t188 * (t135 * t187 - t291);
t67 = t187 * t284;
t107 = -qJD(4) * t158 + t209;
t142 = t158 * t175;
t79 = t107 + t142;
t295 = t184 * t79 - t67;
t91 = t273 + t274;
t304 = t188 * t91;
t297 = t185 * t295 + t304;
t321 = t271 * t297;
t320 = t135 * t184 + t289;
t319 = t187 * t136 - t260;
t315 = -pkin(3) * t283 + t335;
t257 = t187 * t283;
t117 = t274 - t273;
t302 = t185 * t117;
t312 = t188 * (-t184 * t80 + t257) + t302;
t309 = pkin(3) * t91;
t306 = (t106 + t283) * pkin(4);
t305 = qJ(5) * t91;
t265 = t184 * t283;
t300 = t187 * t80 + t265;
t298 = pkin(8) * t295 + t309;
t66 = t184 * t284;
t40 = -t187 * t79 - t66;
t293 = qJ(5) * t282;
t191 = qJD(1) ^ 2;
t285 = t271 * t191;
t78 = t107 - t142;
t131 = pkin(5) * t160 - qJ(6) * t175;
t272 = -2 * qJD(6);
t279 = -t106 * qJ(6) + t160 * t131 + t158 * t272;
t267 = pkin(4) + qJ(6);
t69 = qJ(5) * t284;
t278 = -t267 * t79 - t69;
t88 = qJ(5) * t281;
t277 = -t267 * t280 + t88;
t111 = pkin(4) * t158 - qJ(5) * t160;
t231 = qJD(2) * qJD(1);
t179 = 0.2e1 * t231;
t181 = qJDD(1) * qJ(2);
t186 = sin(qJ(1));
t189 = cos(qJ(1));
t214 = t189 * g(1) + t186 * g(2);
t207 = -t181 + t214;
t205 = t179 - t207;
t211 = -t164 + t224;
t212 = -t163 + t223;
t63 = t212 * pkin(3) + t211 * pkin(8) + t205 - t285;
t221 = t186 * g(1) - t189 * g(2);
t213 = qJDD(2) - t221;
t238 = t191 * qJ(2);
t203 = t213 - t238;
t140 = -t271 * qJDD(1) + t203;
t113 = t188 * g(3) - t185 * t140;
t190 = qJD(3) ^ 2;
t216 = pkin(3) * t185 - pkin(8) * t188;
t206 = t191 * t216;
t83 = -t190 * pkin(3) + qJDD(3) * pkin(8) - t185 * t206 - t113;
t37 = t184 * t83 - t187 * t63;
t24 = -t157 * pkin(4) - t174 * qJ(5) + t160 * t111 + qJDD(5) + t37;
t198 = t107 * pkin(5) + qJ(6) * t97 + t24;
t14 = (pkin(5) * t158 + t272) * t175 + t198;
t171 = 0.2e1 * t234;
t38 = t184 * t63 + t187 * t83;
t217 = t174 * pkin(4) - t157 * qJ(5) + t158 * t111 - t38;
t201 = -t106 * pkin(5) - qJ(6) * t274 + t175 * t131 + qJDD(6) - t217;
t16 = t171 + t201;
t276 = qJ(5) * t16 - t267 * t14;
t275 = t267 * t97 - t293;
t112 = t185 * g(3) + t188 * t140;
t82 = qJDD(3) * pkin(3) + t190 * pkin(8) - t188 * t206 + t112;
t196 = -pkin(4) * t246 + 0.2e1 * qJD(5) * t160 + t82;
t193 = qJ(5) * t78 + t196;
t245 = t175 * t184;
t133 = t160 * t245;
t227 = t185 * t247;
t218 = t188 * (t187 * t107 - t133) + t227;
t268 = t106 * pkin(4);
t19 = (t78 - t80) * qJ(5) + t196 - t268;
t270 = pkin(4) * t184;
t269 = pkin(4) * t187;
t262 = t184 * t82;
t255 = t187 * t82;
t252 = qJ(5) * t187;
t251 = qJDD(1) * pkin(1);
t244 = t175 * t187;
t182 = t185 ^ 2;
t243 = t182 * t191;
t183 = t188 ^ 2;
t242 = t183 * t191;
t241 = t185 * t157;
t226 = t185 * t191 * t188;
t240 = t185 * (qJDD(3) + t226);
t239 = t188 * (qJDD(3) - t226);
t236 = t182 + t183;
t228 = t158 * t244;
t225 = -pkin(4) * t75 - t69;
t222 = qJ(5) * t184 + pkin(3);
t18 = t184 * t37 + t187 * t38;
t219 = t188 * (t106 * t184 + t228) - t227;
t65 = t184 * t107 + t160 * t244;
t23 = t171 - t217;
t215 = -pkin(4) * t24 + qJ(5) * t23;
t64 = -t187 * t106 + t158 * t245;
t210 = t184 * t38 - t187 * t37;
t62 = t188 * t112 - t185 * t113;
t208 = qJ(2) + t216;
t204 = (-t158 * t184 - t160 * t187) * t175;
t202 = -pkin(4) * t280 - t217 + t88;
t200 = t188 * (t133 - t228) + t241;
t197 = pkin(4) * t97 + t24 - t293;
t195 = t175 * t272 + t198;
t192 = t193 + t279;
t166 = t236 * qJDD(1);
t165 = t177 - 0.2e1 * t224;
t162 = 0.2e1 * t223 + t230;
t145 = -t203 + t251;
t134 = t207 - 0.2e1 * t231 + t285;
t130 = -t240 + t188 * (-t190 - t242);
t129 = t185 * (-t190 - t243) + t239;
t45 = pkin(5) * t97 - qJ(5) * t283;
t42 = t184 * t75 - t67;
t39 = -t187 * t75 - t66;
t29 = pkin(5) * t281 - t267 * t80;
t26 = t185 * t42 + t304;
t25 = t193 - t268;
t22 = t24 + t305;
t21 = pkin(4) * t91 + t23;
t20 = -t193 + t306;
t15 = pkin(5) * t274 + t192 - t268;
t13 = t18 * t185 + t188 * t82;
t12 = (t280 + t274) * pkin(5) + t19 + t279;
t11 = t184 * t24 + t187 * t23;
t10 = t184 * t23 - t187 * t24;
t9 = t305 + (t79 + t142) * pkin(5) + t195;
t8 = -pkin(5) * t284 + t267 * t91 + t16;
t7 = (t282 + t274) * pkin(5) - t306 + t192 - qJ(6) * t283;
t6 = t11 * t185 + t188 * t25;
t5 = t14 * t184 + t16 * t187;
t4 = -t14 * t187 + t16 * t184;
t3 = pkin(5) * t14 + qJ(5) * t15;
t2 = pkin(5) * t16 + t267 * t15;
t1 = t15 * t188 + t185 * t5;
t17 = [0, 0, 0, 0, 0, qJDD(1), t221, t214, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t213 - 0.2e1 * t251, t179 + 0.2e1 * t181 - t214, pkin(1) * t145 + qJ(2) * (-t191 * pkin(1) + t205), -t211 * t188, -t162 * t188 - t165 * t185, t239 - t185 * (t190 - t242), t212 * t185, t188 * (-t190 + t243) - t240, 0, qJ(2) * t162 - t271 * t129 - t185 * t134, qJ(2) * t165 - t271 * t130 - t188 * t134, t271 * t166 - t236 * t238 - t62, -qJ(2) * t134 - t271 * t62, t218, t188 * (-t184 * t78 - t257) - t302, t185 * t79 - t314, t219, t324, t200, t188 * (-t262 + t307) - t185 * (t37 + t308) + t323, t188 * (-t255 + t330) - t185 * (t38 + t331) - t332, -t188 * t210 + t208 * t40 - t321, -t271 * t13 + t208 * t210, t241 + t188 * (-t158 * t187 + t160 * t184) * t175, -t333, -t324, t218, -t312, t219, t188 * (-pkin(8) * t39 - t184 * t21 + t187 * t22) - t185 * (-pkin(3) * t39 - t225) + qJ(2) * t39 - t271 * t26, t188 * (-t184 * t20 + t252 * t283 - t307) - t185 * (-t197 - t308) - t323, t188 * (t187 * t19 + t270 * t80 - t330) - t185 * (-t202 + t327) + t332, t188 * (-pkin(8) * t10 + (t252 - t270) * t25) - t185 * (-pkin(3) * t10 - t215) + qJ(2) * t10 - t271 * t6, t200, -t324, t333, t219, t312, t218, t188 * (-pkin(8) * t40 - t184 * t8 + t187 * t9) - t185 * (-pkin(3) * t40 - t278) + qJ(2) * t40 - t321, t188 * (t12 * t187 - t184 * t29 - t330) - t185 * (-t201 - t277 + t327) + t332, t188 * (-t184 * t7 + t187 * t45 + t307) - t185 * (t14 + t275 + t308) + t323, t188 * (-pkin(8) * t4 - t184 * t2 + t187 * t3) - t185 * (-pkin(3) * t4 - t276) + qJ(2) * t4 - t271 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t191, -t145, 0, 0, 0, 0, 0, 0, t129, t130, -t166, t62, 0, 0, 0, 0, 0, 0, t313, -t325, t297, t13, 0, 0, 0, 0, 0, 0, t26, -t313, t325, t6, 0, 0, 0, 0, 0, 0, t297, t325, t313, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, (-t182 + t183) * t191, t177, -t226, -t230, qJDD(3), t112, t113, 0, 0, t65, t187 * t78 - t265, t319, t64, t320, t204, t255 + t315, -t262 + t326, t18 + t298, pkin(3) * t82 + pkin(8) * t18, t204, -t319, -t320, t65, -t300, t64, pkin(8) * t42 + t184 * t22 + t187 * t21 + t309, t187 * t20 + t222 * t283 - t335, t329 + t184 * t19 - (pkin(3) + t269) * t80, pkin(8) * t11 + (t222 + t269) * t25, t204, -t320, t319, t64, t300, t65, t184 * t9 + t187 * t8 + t298, t12 * t184 + t187 * t29 - t326, t184 * t45 + t187 * t7 + t315, pkin(3) * t15 + pkin(8) * t5 + t184 * t3 + t187 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, -t117, t79, -t247, -t284, t157, -t37, -t38, 0, 0, t157, -t75, t284, t247, -t117, -t247, t225, t197, t171 + t202, t215, t157, t284, t75, -t247, t117, t247, t278, t16 + t277, -pkin(5) * t142 - t195 - t275, t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t97, t280, t24, 0, 0, 0, 0, 0, 0, t79, t280, t97, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t284, t281, t282, t16;];
tauJ_reg  = t17;
