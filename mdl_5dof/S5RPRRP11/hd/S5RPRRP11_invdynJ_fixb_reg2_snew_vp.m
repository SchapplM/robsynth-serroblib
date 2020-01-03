% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:33
% EndTime: 2019-12-31 18:54:44
% DurationCPUTime: 4.74s
% Computational Cost: add. (8868->343), mult. (21181->445), div. (0->0), fcn. (15427->8), ass. (0->221)
t185 = sin(pkin(8));
t186 = cos(pkin(8));
t189 = sin(qJ(3));
t192 = cos(qJ(3));
t251 = t185 * t189;
t170 = (-t186 * t192 + t251) * qJD(1);
t165 = qJD(4) + t170;
t276 = t165 ^ 2;
t250 = t186 * t189;
t210 = t185 * t192 + t250;
t172 = t210 * qJD(1);
t188 = sin(qJ(4));
t191 = cos(qJ(4));
t154 = -t191 * qJD(3) + t188 * t172;
t277 = t154 ^ 2;
t124 = t277 - t276;
t234 = t186 * qJDD(1);
t235 = t185 * qJDD(1);
t212 = t189 * t235 - t192 * t234;
t237 = t172 * qJD(3);
t144 = -t212 - t237;
t135 = qJDD(4) - t144;
t156 = t188 * qJD(3) + t191 * t172;
t254 = t156 * t154;
t94 = -t254 - t135;
t266 = t188 * t94;
t68 = -t191 * t124 - t266;
t129 = t165 * t156;
t169 = t210 * qJDD(1);
t238 = t170 * qJD(3);
t146 = t169 - t238;
t221 = t191 * qJDD(3) - t188 * t146;
t205 = t156 * qJD(4) - t221;
t77 = -t129 + t205;
t342 = t185 * (t189 * t77 + t192 * t68) + t186 * (t189 * t68 - t192 * t77);
t153 = t156 ^ 2;
t291 = -t153 - t276;
t57 = t191 * t291 + t266;
t341 = pkin(1) * t57;
t340 = pkin(2) * t57;
t339 = pkin(3) * t57;
t338 = pkin(7) * t57;
t260 = t191 * t94;
t59 = -t188 * t291 + t260;
t337 = pkin(7) * t59;
t336 = t189 * t59;
t335 = t192 * t59;
t290 = t153 - t277;
t209 = -t188 * qJDD(3) - t191 * t146;
t203 = -t154 * qJD(4) - t209;
t255 = t154 * t165;
t297 = t255 - t203;
t267 = t188 * t297;
t292 = t129 + t205;
t46 = t191 * t292 - t267;
t332 = t185 * (-t189 * t290 + t192 * t46) + t186 * (t189 * t46 + t192 * t290);
t288 = -t254 + t135;
t265 = t188 * t288;
t285 = -t276 - t277;
t295 = t191 * t285 - t265;
t313 = t189 * t295 - t192 * t292;
t329 = pkin(6) * t313;
t328 = qJ(5) * t297;
t327 = -t188 * t124 + t260;
t259 = t191 * t288;
t296 = t188 * t285 + t259;
t312 = t189 * t292 + t192 * t295;
t326 = -pkin(2) * t296 + pkin(6) * t312;
t325 = qJ(2) * (-t185 * t313 + t186 * t312) - pkin(1) * t296;
t286 = t255 + t203;
t125 = -t153 + t276;
t314 = -t188 * t125 + t259;
t324 = t185 * (t189 * t286 + t192 * t314) + t186 * (t189 * t314 - t192 * t286);
t321 = pkin(3) * t296;
t320 = pkin(7) * t295;
t319 = pkin(7) * t296;
t316 = t191 * t125 + t265;
t182 = t185 ^ 2;
t183 = t186 ^ 2;
t240 = t182 + t183;
t190 = sin(qJ(1));
t274 = cos(qJ(1));
t208 = t274 * g(1) + t190 * g(2);
t315 = 0.2e1 * qJD(2) * qJD(1) - t208;
t289 = t153 + t277;
t311 = pkin(3) * t289;
t148 = t172 * t170;
t284 = qJDD(3) - t148;
t309 = t189 * t284;
t307 = t189 * t289;
t303 = t192 * t284;
t301 = t192 * t289;
t193 = qJD(1) ^ 2;
t298 = -t193 * pkin(1) + qJDD(1) * qJ(2) + t315;
t294 = -t188 * t292 - t191 * t297;
t224 = t190 * g(1) - t274 * g(2);
t213 = -qJDD(2) + t224;
t227 = pkin(2) * t186 + pkin(1);
t139 = t227 * qJDD(1) + (t240 * pkin(6) + qJ(2)) * t193 + t213;
t293 = pkin(6) + qJ(2);
t242 = t193 * qJ(2);
t256 = qJDD(1) * pkin(1);
t166 = t213 + t242 + t256;
t283 = t240 * t242 - t166 - t256;
t112 = t154 * pkin(4) - t156 * qJ(5);
t136 = t170 * pkin(3) - t172 * pkin(7);
t197 = (-t293 * qJDD(1) + t227 * t193 - t315) * t185;
t272 = t186 * g(3);
t196 = t197 - t272;
t219 = -t185 * g(3) + t186 * t298;
t127 = -t183 * t193 * pkin(2) + pkin(6) * t234 + t219;
t245 = t192 * t127;
t275 = qJD(3) ^ 2;
t62 = -t275 * pkin(3) + qJDD(3) * pkin(7) - t170 * t136 + t189 * t196 + t245;
t65 = (-t146 + t238) * pkin(7) + (-t144 + t237) * pkin(3) - t139;
t35 = t188 * t65 + t191 * t62;
t223 = -t135 * qJ(5) + t154 * t112 - t35;
t282 = -(t291 + t276) * pkin(4) - qJ(5) * t94 - t223;
t252 = t165 * t191;
t229 = t154 * t252;
t206 = t188 * t205 + t229;
t228 = t192 * t254;
t230 = t189 * t254;
t279 = t185 * (t192 * t206 - t230) + t186 * (t189 * t206 + t228);
t253 = t165 * t188;
t122 = t156 * t253;
t214 = t122 - t229;
t278 = t185 * (t189 * t135 + t192 * t214) + t186 * (-t192 * t135 + t189 * t214);
t167 = t170 ^ 2;
t168 = t172 ^ 2;
t273 = pkin(4) * t191;
t86 = t189 * t127 - t192 * t196;
t87 = -g(3) * t250 + t189 * t197 + t245;
t50 = t189 * t87 - t192 * t86;
t271 = t185 * t50;
t61 = -qJDD(3) * pkin(3) - t275 * pkin(7) + t172 * t136 + t86;
t270 = t188 * t61;
t268 = t188 * t286;
t263 = t191 * t61;
t261 = t191 * t286;
t257 = qJ(5) * t191;
t247 = t189 * t139;
t141 = qJDD(3) + t148;
t246 = t189 * t141;
t244 = t192 * t139;
t243 = t192 * t141;
t239 = qJD(5) * t165;
t226 = -pkin(3) * t192 - pkin(2);
t225 = -qJ(5) * t188 - pkin(3);
t34 = t188 * t62 - t191 * t65;
t16 = t188 * t34 + t191 * t35;
t51 = t189 * t86 + t192 * t87;
t220 = t185 * (t185 * t298 + t272) + t186 * t219;
t157 = 0.2e1 * t239;
t218 = t157 - t223;
t24 = -pkin(4) * t276 + t218;
t25 = -t135 * pkin(4) - qJ(5) * t276 + t156 * t112 + qJDD(5) + t34;
t217 = -pkin(4) * t25 + qJ(5) * t24;
t216 = -pkin(4) * t286 - qJ(5) * t77;
t215 = t154 * t253 - t191 * t205;
t15 = t188 * t35 - t191 * t34;
t204 = (-t154 * t188 - t156 * t191) * t165;
t201 = t205 * pkin(4) + t328 + t61;
t200 = 0.2e1 * qJD(5) * t156 - t201;
t75 = t191 * t203 - t122;
t199 = t185 * (t192 * t75 + t230) + t186 * (t189 * t75 - t228);
t198 = pkin(4) * t288 + qJ(5) * t285 - t25;
t178 = t183 * qJDD(1);
t177 = t182 * qJDD(1);
t173 = t240 * t193;
t160 = -t168 - t275;
t159 = -t168 + t275;
t158 = t167 - t275;
t145 = t169 - 0.2e1 * t238;
t143 = t212 + 0.2e1 * t237;
t137 = -t275 - t167;
t116 = -t167 - t168;
t109 = -t189 * t160 - t243;
t108 = t192 * t160 - t246;
t99 = t189 * t169 - t192 * t212;
t98 = -t192 * t169 - t189 * t212;
t97 = t192 * t137 - t309;
t96 = t189 * t137 + t303;
t84 = (qJD(4) + t165) * t154 + t209;
t79 = (-qJD(4) + t165) * t156 + t221;
t74 = t156 * t252 + t188 * t203;
t49 = t191 * t79 + t268;
t47 = -t191 * t77 + t268;
t45 = t188 * t79 - t261;
t44 = -t188 * t77 - t261;
t42 = -t189 * t84 + t335;
t40 = t192 * t84 + t336;
t38 = t189 * t297 - t335;
t36 = -t192 * t297 - t336;
t32 = t263 - t338;
t31 = t192 * t49 - t307;
t30 = t192 * t47 - t307;
t29 = t189 * t49 + t301;
t28 = t189 * t47 + t301;
t27 = t270 - t319;
t26 = (pkin(4) * t165 - 0.2e1 * qJD(5)) * t156 + t201;
t23 = t35 - t339;
t22 = -pkin(3) * t44 - t216;
t21 = t34 - t321;
t20 = qJ(5) * t289 + t25;
t19 = (-t292 - t129) * pkin(4) + t200;
t18 = -pkin(4) * t129 + t200 - t328;
t17 = (t289 - t276) * pkin(4) + t218;
t14 = -t198 - t321;
t11 = -0.2e1 * t239 - t282 + t339;
t10 = -t188 * t19 - t257 * t292 - t319;
t9 = pkin(4) * t267 + t191 * t18 + t338;
t8 = -pkin(7) * t45 - t15;
t7 = t188 * t25 + t191 * t24;
t6 = t188 * t24 - t191 * t25;
t5 = -pkin(7) * t44 - t188 * t17 + t191 * t20;
t4 = t189 * t26 + t192 * t7;
t3 = t189 * t7 - t192 * t26;
t2 = -pkin(7) * t6 + (pkin(4) * t188 - t257) * t26;
t1 = -pkin(3) * t6 - t217;
t12 = [0, 0, 0, 0, 0, qJDD(1), t224, t208, 0, 0, t177, 0.2e1 * t185 * t234, 0, t178, 0, 0, -t283 * t186, t283 * t185, pkin(1) * t173 + qJ(2) * (t178 + t177) + t220, pkin(1) * t166 + qJ(2) * t220, t185 * (t192 * t146 - t189 * t237) + t186 * (t189 * t146 + t192 * t237), t185 * (-t192 * t143 - t189 * t145) + t186 * (-t189 * t143 + t192 * t145), t185 * (-t189 * t159 + t303) + t186 * (t192 * t159 + t309), t185 * (-t189 * t144 + t192 * t238) + t186 * (t192 * t144 + t189 * t238), t185 * (t192 * t158 - t246) + t186 * (t189 * t158 + t243), (t185 * (-t170 * t192 + t172 * t189) + t186 * (-t170 * t189 - t172 * t192)) * qJD(3), t185 * (-pkin(6) * t96 - t247) + t186 * (-pkin(2) * t143 + pkin(6) * t97 + t244) - pkin(1) * t143 + qJ(2) * (-t185 * t96 + t186 * t97), t185 * (-pkin(6) * t108 - t244) + t186 * (-pkin(2) * t145 + pkin(6) * t109 - t247) - pkin(1) * t145 + qJ(2) * (-t185 * t108 + t186 * t109), t185 * (-pkin(6) * t98 - t50) + t186 * (-pkin(2) * t116 + pkin(6) * t99 + t51) - pkin(1) * t116 + qJ(2) * (-t185 * t98 + t186 * t99), -pkin(6) * t271 + t186 * (pkin(2) * t139 + pkin(6) * t51) + pkin(1) * t139 + qJ(2) * (t186 * t51 - t271), t199, -t332, t324, t279, -t342, t278, t185 * (-t189 * t21 + t192 * t27 - t329) + t186 * (t189 * t27 + t192 * t21 + t326) + t325, t185 * (-pkin(6) * t40 - t189 * t23 + t192 * t32) + t186 * (pkin(6) * t42 + t189 * t32 + t192 * t23 - t340) - t341 + qJ(2) * (-t185 * t40 + t186 * t42), t185 * (-pkin(6) * t29 + t192 * t8) + t186 * (pkin(6) * t31 + t189 * t8) + qJ(2) * (-t185 * t29 + t186 * t31) + (pkin(3) * t251 + t186 * t226 - pkin(1)) * t45, (t185 * (pkin(3) * t189 - pkin(7) * t192) + t186 * (-pkin(7) * t189 + t226) - pkin(1)) * t15 + t293 * (-t185 * (t189 * t16 - t192 * t61) + t186 * (t192 * t16 + t189 * t61)), t199, t324, t332, t278, t342, t279, t185 * (t192 * t10 - t189 * t14 - t329) + t186 * (t189 * t10 + t192 * t14 + t326) + t325, t185 * (-pkin(6) * t28 - t189 * t22 + t192 * t5) + t186 * (-pkin(2) * t44 + pkin(6) * t30 + t189 * t5 + t192 * t22) - pkin(1) * t44 + qJ(2) * (-t185 * t28 + t186 * t30), t185 * (-pkin(6) * t36 - t189 * t11 + t192 * t9) + t186 * (pkin(6) * t38 + t192 * t11 + t189 * t9 + t340) + t341 + qJ(2) * (-t185 * t36 + t186 * t38), t185 * (-pkin(6) * t3 - t189 * t1 + t192 * t2) + t186 * (-pkin(2) * t6 + pkin(6) * t4 + t192 * t1 + t189 * t2) - pkin(1) * t6 + qJ(2) * (-t185 * t3 + t186 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t235, -t173, -t166, 0, 0, 0, 0, 0, 0, t143, t145, t116, -t139, 0, 0, 0, 0, 0, 0, t296, t57, t45, t15, 0, 0, 0, 0, 0, 0, t296, t44, -t57, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t168 - t167, t169, -t148, -t212, qJDD(3), -t86, -t87, 0, 0, t74, t294, t316, t215, -t327, t204, -pkin(3) * t292 - t263 + t320, pkin(3) * t84 + t270 + t337, pkin(7) * t49 + t16 + t311, -pkin(3) * t61 + pkin(7) * t16, t74, t316, -t294, t204, t327, t215, t191 * t19 + t225 * t292 + t320, pkin(7) * t47 + t191 * t17 + t188 * t20 + t311, -t337 + t188 * t18 - (pkin(3) + t273) * t297, pkin(7) * t7 + (t225 - t273) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, t290, t286, -t254, -t77, t135, -t34, -t35, 0, 0, t254, t286, -t290, t135, t77, -t254, t198, t216, t157 + t282, t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, t286, t291, t25;];
tauJ_reg = t12;
