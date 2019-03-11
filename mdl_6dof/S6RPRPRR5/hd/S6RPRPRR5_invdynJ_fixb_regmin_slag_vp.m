% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:39
% EndTime: 2019-03-09 03:49:48
% DurationCPUTime: 3.83s
% Computational Cost: add. (4809->395), mult. (11356->480), div. (0->0), fcn. (9058->12), ass. (0->215)
t175 = qJDD(3) - qJDD(5);
t179 = qJD(3) - qJD(5);
t188 = sin(qJ(6));
t192 = cos(qJ(6));
t264 = qJD(6) * t192;
t265 = qJD(6) * t188;
t185 = cos(pkin(10));
t303 = cos(qJ(3));
t256 = t303 * t185;
t235 = qJD(1) * t256;
t184 = sin(pkin(10));
t190 = sin(qJ(3));
t275 = t190 * t184;
t255 = qJD(1) * t275;
t125 = -t235 + t255;
t138 = t303 * t184 + t190 * t185;
t127 = t138 * qJD(1);
t189 = sin(qJ(5));
t193 = cos(qJ(5));
t266 = qJD(5) * t193;
t267 = qJD(5) * t189;
t247 = qJDD(1) * t303;
t261 = t185 * qJDD(1);
t257 = qJD(3) * t235 + t184 * t247 + t190 * t261;
t80 = qJD(3) * t255 - t257;
t130 = t138 * qJD(3);
t262 = t184 * qJDD(1);
t225 = -t185 * t247 + t190 * t262;
t81 = qJD(1) * t130 + t225;
t27 = t125 * t266 - t127 * t267 + t189 * t81 - t193 * t80;
t313 = t189 * t125 + t193 * t127;
t10 = -t188 * t175 - t179 * t264 + t192 * t27 - t265 * t313;
t223 = t188 * t179 - t192 * t313;
t11 = -qJD(6) * t223 + t192 * t175 + t188 * t27;
t315 = -t193 * t125 + t189 * t127;
t319 = qJD(6) + t315;
t296 = t223 * t319;
t56 = t192 * t179 + t188 * t313;
t297 = t56 * t319;
t336 = t192 * (t10 - t297) + (-t11 + t296) * t188;
t295 = t223 * t313;
t28 = qJD(5) * t313 - t189 * t80 - t193 * t81;
t24 = qJDD(6) + t28;
t291 = t188 * t24;
t329 = t192 * t319;
t327 = t319 * t329 + t291;
t333 = t295 + t327;
t292 = t10 * t188;
t332 = t223 * t329 - t292;
t288 = t313 * t179;
t330 = t28 + t288;
t293 = pkin(7) + qJ(2);
t148 = t293 * t185;
t140 = qJD(1) * t148;
t122 = t190 * t140;
t147 = t293 * t184;
t139 = qJD(1) * t147;
t84 = -t303 * t139 - t122;
t272 = qJD(4) - t84;
t248 = -t192 * t24 + t265 * t319;
t290 = t188 * t319;
t236 = t290 * t315 + t248;
t294 = t56 * t313;
t328 = t236 - t294;
t326 = t313 * pkin(5);
t166 = t185 * pkin(2) + pkin(1);
t144 = -t166 * qJD(1) + qJD(2);
t59 = t125 * pkin(3) - t127 * qJ(4) + t144;
t41 = -t125 * pkin(4) - t59;
t12 = pkin(5) * t315 - pkin(9) * t313 + t41;
t195 = -pkin(3) - pkin(4);
t320 = -t127 * pkin(8) + t272;
t47 = t195 * qJD(3) + t320;
t182 = qJD(3) * qJ(4);
t85 = -t190 * t139 + t303 * t140;
t55 = t125 * pkin(8) + t85;
t48 = t182 + t55;
t26 = t189 * t47 + t193 * t48;
t21 = -t179 * pkin(9) + t26;
t5 = t192 * t12 - t188 * t21;
t325 = t313 * t5;
t6 = t188 * t12 + t192 * t21;
t324 = t313 * t6;
t323 = t319 * t313;
t287 = t315 * t179;
t322 = -t27 + t287;
t321 = t313 * t315;
t183 = qJDD(1) * pkin(1);
t170 = qJDD(2) - t183;
t194 = cos(qJ(1));
t174 = g(2) * t194;
t191 = sin(qJ(1));
t314 = g(1) * t191 - t174;
t215 = -t170 + t314;
t318 = t313 ^ 2 - t315 ^ 2;
t263 = qJD(1) * qJD(2);
t308 = t293 * qJDD(1) + t263;
t100 = t308 * t184;
t101 = t308 * t185;
t254 = qJD(3) * t303;
t268 = qJD(3) * t190;
t239 = t303 * t100 + t190 * t101 - t139 * t268 + t140 * t254;
t224 = qJDD(4) + t239;
t18 = t80 * pkin(8) + t195 * qJDD(3) + t224;
t180 = qJDD(3) * qJ(4);
t181 = qJD(3) * qJD(4);
t260 = -t190 * t100 + t303 * t101 - t139 * t254;
t34 = -t140 * t268 + t180 + t181 + t260;
t19 = t81 * pkin(8) + t34;
t213 = t189 * t18 + t193 * t19 + t47 * t266 - t267 * t48;
t178 = pkin(10) + qJ(3);
t172 = cos(t178);
t171 = sin(t178);
t281 = t171 * t193;
t217 = t171 * t189 + t172 * t193;
t95 = t217 * t191;
t232 = g(2) * t95 + g(3) * (-t172 * t189 + t281);
t97 = t217 * t194;
t317 = g(1) * t97 + t315 * t41 - t213 + t232;
t279 = t172 * t191;
t94 = t189 * t279 - t191 * t281;
t278 = t172 * t194;
t280 = t171 * t194;
t96 = t189 * t278 - t193 * t280;
t212 = g(1) * t96 + g(2) * t94 + g(3) * t217;
t221 = -t193 * t18 + t189 * t19 + t48 * t266 + t47 * t267;
t316 = -t313 * t41 + t212 - t221;
t89 = -t190 * t147 + t303 * t148;
t270 = t193 * qJ(4) + t189 * t195;
t300 = g(2) * t191;
t301 = g(1) * t194;
t231 = t300 + t301;
t312 = qJ(2) * qJDD(1);
t142 = -pkin(9) + t270;
t2 = t175 * pkin(5) + t221;
t211 = -t2 + t212;
t77 = t127 * pkin(3) + t125 * qJ(4);
t49 = -t127 * pkin(4) - t77;
t311 = (-pkin(9) * t315 + qJD(6) * t142 - t326 + t49) * t319 + t211;
t310 = (t319 * pkin(9) + t326) * t319 - t211;
t60 = -t147 * t254 + qJD(2) * t256 + (-qJD(2) * t184 - qJD(3) * t148) * t190;
t309 = -t60 * qJD(3) - t89 * qJDD(3) - t171 * t314;
t25 = -t189 * t48 + t193 * t47;
t20 = t179 * pkin(5) - t25;
t137 = -t256 + t275;
t218 = t193 * t137 - t189 * t138;
t78 = t137 * pkin(3) - t138 * qJ(4) - t166;
t53 = -t137 * pkin(4) - t78;
t83 = t189 * t137 + t193 * t138;
t22 = -pkin(5) * t218 - t83 * pkin(9) + t53;
t251 = -t175 * pkin(9) + qJD(6) * t12 + t213;
t88 = t303 * t147 + t190 * t148;
t63 = -t138 * pkin(8) + t88;
t64 = t137 * pkin(8) + t89;
t33 = t189 * t63 + t193 * t64;
t129 = t184 * t268 - t185 * t254;
t38 = qJD(5) * t218 - t193 * t129 + t189 * t130;
t32 = t189 * t64 - t193 * t63;
t42 = t130 * pkin(8) + t60;
t61 = t138 * qJD(2) + qJD(3) * t89;
t43 = t129 * pkin(8) + t61;
t7 = -qJD(5) * t32 + t189 * t43 + t193 * t42;
t307 = t2 * t83 + t20 * t38 - t33 * t24 - (qJD(6) * t22 + t7) * t319 + t251 * t218 + t301;
t306 = g(3) * t171 + qJD(3) * (t84 + t122) + t172 * t231 - t260;
t305 = t127 ^ 2;
t304 = g(1) * t95;
t298 = t20 * t83;
t220 = -t189 * qJ(4) + t193 * t195;
t286 = qJD(5) * t220 - t189 * t55 + t320 * t193;
t285 = t270 * qJD(5) + t320 * t189 + t193 * t55;
t284 = qJD(6) * t21;
t283 = qJDD(3) * pkin(3);
t282 = t127 * t125;
t273 = t85 * qJD(3);
t269 = t184 ^ 2 + t185 ^ 2;
t259 = t83 * t265;
t258 = t319 * t264;
t52 = t130 * pkin(3) + t129 * qJ(4) - t138 * qJD(4);
t143 = -pkin(2) * t261 + t170;
t242 = t269 * qJD(1) ^ 2;
t240 = t179 * t319;
t238 = t179 ^ 2;
t237 = 0.2e1 * t269;
t233 = t22 * t24 + t304;
t229 = t24 * t83 + t319 * t38;
t228 = t172 * pkin(3) + t171 * qJ(4);
t226 = -t284 - t174;
t40 = -t130 * pkin(4) - t52;
t214 = t166 + t228;
t31 = t81 * pkin(3) + t80 * qJ(4) - t127 * qJD(4) + t143;
t210 = t232 - t251;
t209 = g(1) * t280 - g(3) * t172 + t171 * t300 - t239;
t208 = t215 + t183;
t13 = -t81 * pkin(4) - t31;
t207 = -pkin(9) * t24 + (t20 + t25) * t319;
t205 = g(1) * t279 - g(2) * t278 - t61 * qJD(3) - t88 * qJDD(3);
t204 = -t142 * t24 + (-t20 - t286) * t319;
t202 = t237 * t263 - t231;
t200 = t59 * t127 + qJDD(4) - t209;
t198 = 0.2e1 * t127 * qJD(3) + t225;
t141 = pkin(5) - t220;
t119 = t125 ^ 2;
t87 = -t191 * t188 + t97 * t192;
t86 = -t97 * t188 - t191 * t192;
t75 = t182 + t85;
t74 = -qJD(3) * pkin(3) + t272;
t51 = (t125 - t255) * qJD(3) + t257;
t50 = (t125 + t255) * qJD(3) - t257;
t39 = qJD(5) * t83 - t189 * t129 - t193 * t130;
t35 = t224 - t283;
t9 = t39 * pkin(5) - t38 * pkin(9) + t40;
t8 = qJD(5) * t33 + t189 * t42 - t193 * t43;
t4 = t28 * pkin(5) - t27 * pkin(9) + t13;
t3 = t192 * t4;
t1 = [qJDD(1), t314, t231, t208 * t185, -t208 * t184, t237 * t312 + t202, pkin(1) * t215 + (t269 * t312 + t202) * qJ(2), -t127 * t129 - t80 * t138, t129 * t125 - t127 * t130 + t80 * t137 - t138 * t81, -t129 * qJD(3) + t138 * qJDD(3), -t130 * qJD(3) - t137 * qJDD(3), 0, t144 * t130 + t143 * t137 - t166 * t81 + t205, -t144 * t129 + t143 * t138 + t166 * t80 + t309, t52 * t125 + t59 * t130 + t31 * t137 + t78 * t81 + t205, -t60 * t125 + t61 * t127 - t74 * t129 - t75 * t130 - t34 * t137 + t35 * t138 - t88 * t80 - t89 * t81 - t231, -t52 * t127 + t59 * t129 - t31 * t138 + t78 * t80 - t309, t31 * t78 + t34 * t89 + t35 * t88 + t59 * t52 + t75 * t60 + t74 * t61 + (-g(1) * t293 - g(2) * t214) * t194 + (g(1) * t214 - g(2) * t293) * t191, t27 * t83 + t313 * t38, t218 * t27 - t83 * t28 - t313 * t39 - t315 * t38, -t83 * t175 - t38 * t179, -t175 * t218 + t39 * t179, 0, -g(2) * t97 - t13 * t218 + t32 * t175 + t8 * t179 + t53 * t28 + t315 * t40 + t41 * t39 + t304, -g(1) * t94 + g(2) * t96 + t13 * t83 + t33 * t175 + t7 * t179 + t53 * t27 + t313 * t40 + t41 * t38, t223 * t259 + (t10 * t83 - t223 * t38) * t192 (t188 * t223 - t192 * t56) * t38 + (-t292 - t11 * t192 + (t188 * t56 + t192 * t223) * qJD(6)) * t83, -t10 * t218 + t192 * t229 - t223 * t39 - t259 * t319, t11 * t218 - t188 * t229 - t258 * t83 - t56 * t39, -t218 * t24 + t319 * t39, -g(2) * t87 + t32 * t11 - t3 * t218 + t5 * t39 + t8 * t56 + (t9 * t319 + (t21 * t218 - t319 * t33 + t298) * qJD(6) + t233) * t192 + t307 * t188, -g(2) * t86 + t32 * t10 - t6 * t39 - t8 * t223 + (-(-qJD(6) * t33 + t9) * t319 + (t4 - t284) * t218 - qJD(6) * t298 - t233) * t188 + t307 * t192; 0, 0, 0, -t261, t262, -t242, -qJ(2) * t242 - t215, 0, 0, 0, 0, 0, t198, -t50, t198, -t119 - t305, t50, t75 * t125 - t74 * t127 + t31 - t314, 0, 0, 0, 0, 0, -t28 + t288, -t27 - t287, 0, 0, 0, 0, 0, t236 + t294, -t295 + t327; 0, 0, 0, 0, 0, 0, 0, t282, -t119 + t305, t51, -t225, qJDD(3), -t144 * t127 + t209 + t273, t144 * t125 + t306, -t77 * t125 - t200 + t273 + 0.2e1 * t283, pkin(3) * t80 - t81 * qJ(4) + (t75 - t85) * t127 + (t74 - t272) * t125, -t59 * t125 + t77 * t127 + 0.2e1 * t180 + 0.2e1 * t181 - t306, -t35 * pkin(3) - g(3) * t228 + t34 * qJ(4) + t272 * t75 - t59 * t77 - t74 * t85 + t231 * (pkin(3) * t171 - qJ(4) * t172) -t321, -t318, t322, t330, t175, -t175 * t220 + t179 * t285 - t315 * t49 - t316, t175 * t270 + t179 * t286 - t313 * t49 - t317, t332, -t336, -t333, t328, t323, t141 * t11 + t204 * t188 - t192 * t311 + t285 * t56 + t325, t141 * t10 + t188 * t311 + t204 * t192 - t223 * t285 - t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t282, t51, -qJD(3) ^ 2 - t305, -t75 * qJD(3) + t200 - t283, 0, 0, 0, 0, 0, -t127 * t315 - t193 * t175 - t189 * t238, -t127 * t313 + t189 * t175 - t193 * t238, 0, 0, 0, 0, 0, -t127 * t329 + (t188 * t240 - t11) * t193 + (-t179 * t56 - t258 - t291) * t189, t127 * t290 + (t192 * t240 - t10) * t193 + (t179 * t223 + t248) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t321, t318, -t322, -t330, -t175, -t26 * t179 + t316, -t25 * t179 + t317, -t332, t336, t333, -t328, -t323, -pkin(5) * t11 + t207 * t188 - t192 * t310 - t26 * t56 - t325, -pkin(5) * t10 + t188 * t310 + t207 * t192 + t223 * t26 + t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223 * t56, t223 ^ 2 - t56 ^ 2, t10 + t297, -t11 - t296, t24, -g(1) * t86 + t188 * t210 + t192 * t226 + t20 * t223 + t319 * t6 + t3, g(1) * t87 + t20 * t56 + t5 * t319 + (-t226 - t4) * t188 + t210 * t192;];
tau_reg  = t1;
