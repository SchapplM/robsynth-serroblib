% Calculate vector of inverse dynamics joint torques for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:19
% DurationCPUTime: 4.58s
% Computational Cost: add. (9636->424), mult. (5845->485), div. (0->0), fcn. (4304->10), ass. (0->253)
t190 = pkin(9) + qJ(5);
t182 = sin(t190);
t183 = cos(t190);
t143 = rSges(6,1) * t182 + rSges(6,2) * t183;
t191 = qJD(1) + qJD(2);
t192 = qJ(1) + qJ(2);
t185 = sin(t192);
t186 = cos(t192);
t242 = rSges(3,1) * t185 + rSges(3,2) * t186;
t116 = t242 * t191;
t196 = sin(qJ(1));
t310 = pkin(1) * qJD(1);
t267 = t196 * t310;
t104 = t116 + t267;
t184 = pkin(8) + t192;
t176 = sin(t184);
t177 = cos(t184);
t273 = qJD(5) * t191;
t110 = qJDD(5) * t176 + t177 * t273;
t315 = rSges(6,1) * t183;
t144 = -rSges(6,2) * t182 + t315;
t124 = t144 * qJD(5);
t189 = qJDD(1) + qJDD(2);
t188 = t191 ^ 2;
t198 = qJD(1) ^ 2;
t197 = cos(qJ(1));
t326 = pkin(1) * t197;
t327 = pkin(1) * t196;
t246 = -qJDD(1) * t326 + t198 * t327;
t325 = pkin(2) * t185;
t225 = t188 * t325 + t246;
t220 = qJDD(4) * t177 + t225;
t295 = qJ(4) * t176;
t323 = pkin(3) * t177;
t125 = t295 + t323;
t324 = pkin(2) * t186;
t255 = -t125 - t324;
t194 = cos(pkin(9));
t178 = pkin(4) * t194 + pkin(3);
t153 = t177 * t178;
t195 = -pkin(7) - qJ(4);
t281 = qJ(4) + t195;
t76 = t176 * t281 - t153 + t323;
t288 = t177 * t182;
t155 = rSges(6,2) * t288;
t287 = t177 * t183;
t269 = rSges(6,1) * t287;
t223 = rSges(6,3) * t176 + t269;
t90 = -t155 + t223;
t233 = t255 + t76 - t90;
t289 = t176 * t191;
t162 = pkin(3) * t289;
t286 = t177 * t191;
t245 = qJ(4) * t286 - t162;
t170 = qJD(4) * t176;
t250 = -0.2e1 * t170 - t245;
t275 = qJD(5) * t176;
t282 = t191 * t195;
t280 = t177 * t282 + t178 * t289;
t290 = t176 * t183;
t270 = rSges(6,1) * t290;
t349 = t143 * qJD(5);
t263 = -t349 * t177 - t191 * t270;
t291 = t176 * t182;
t344 = rSges(6,2) * t291 + t177 * rSges(6,3);
t55 = t191 * t344 + t263;
t8 = t124 * t275 + t110 * t143 + (t245 + t250 - t55 + t280) * t191 + t233 * t189 + t220;
t361 = -g(2) + t8;
t111 = qJDD(5) * t177 - t176 * t273;
t147 = t176 * t282;
t216 = (-qJDD(1) * t196 - t197 * t198) * pkin(1);
t248 = t344 - t325;
t251 = t281 * t177;
t252 = -t178 - t315;
t274 = qJD(5) * t177;
t284 = t186 * t188;
t171 = qJD(4) * t177;
t239 = -pkin(3) * t176 + qJ(4) * t177;
t346 = -t191 * t125 + t171;
t334 = qJDD(4) * t176 + t189 * t239 + (t171 + t346) * t191;
t262 = t191 * t155 + t349 * t176;
t56 = -t191 * t223 + t262;
t7 = -pkin(2) * t284 - t124 * t274 - t111 * t143 + t216 + (-t251 + (pkin(3) + t252) * t176 + t248) * t189 + (t147 + t56 + (t125 - t153) * t191) * t191 + t334;
t360 = -g(3) + t7;
t317 = rSges(5,1) * t194;
t137 = t289 * t317;
t193 = sin(pkin(9));
t314 = rSges(5,2) * t193;
t163 = t177 * t314;
t312 = rSges(5,3) * t176;
t224 = -t177 * t317 - t312;
t93 = -t163 - t224;
t249 = t255 - t93;
t268 = t176 * t314;
t311 = rSges(5,3) * t177;
t23 = t249 * t189 + (t137 + (-t268 - t311) * t191 + t250) * t191 + t220;
t359 = -g(2) + t23;
t172 = t176 * rSges(4,2);
t318 = rSges(4,1) * t177;
t126 = -t172 + t318;
t254 = -t126 - t324;
t276 = -rSges(4,1) * t289 - rSges(4,2) * t286;
t358 = t189 * t254 - t191 * t276 - g(2) + t225;
t146 = rSges(3,1) * t186 - t185 * rSges(3,2);
t357 = t116 * t191 - t146 * t189 - g(2) + t246;
t138 = t191 * t163;
t202 = t216 + (-t185 * t189 - t284) * pkin(2);
t240 = t314 - t317;
t205 = t176 * t240 + t311;
t22 = t189 * t205 + (t191 * t224 + t138) * t191 + t202 + t334;
t356 = -g(3) + t22;
t159 = rSges(4,2) * t289;
t241 = -rSges(4,1) * t176 - rSges(4,2) * t177;
t355 = -g(3) + t189 * t241 + t191 * (-rSges(4,1) * t286 + t159) + t202;
t283 = t186 * t191;
t285 = t185 * t191;
t117 = -rSges(3,1) * t283 + rSges(3,2) * t285;
t354 = t117 * t191 - t189 * t242 - g(3) + t216;
t353 = -t191 * t93 - t138;
t174 = Icges(6,4) * t183;
t231 = -Icges(6,2) * t182 + t174;
t343 = Icges(6,1) * t182 + t174;
t278 = t343 + t231;
t300 = Icges(6,4) * t182;
t133 = Icges(6,2) * t183 + t300;
t136 = Icges(6,1) * t183 - t300;
t279 = t133 - t136;
t352 = (t182 * t278 + t183 * t279) * t191;
t244 = -t318 - t324;
t351 = -t159 + (-t126 - t244) * t191;
t108 = t143 * t275;
t78 = t191 * t90;
t350 = t191 * t76 + t108 - t147 - t171 - t262 - t78;
t85 = Icges(6,4) * t287 - Icges(6,2) * t288 + Icges(6,6) * t176;
t152 = Icges(6,4) * t288;
t87 = Icges(6,1) * t287 + Icges(6,5) * t176 - t152;
t235 = t182 * t85 - t183 * t87;
t348 = t235 * t177;
t347 = -t191 * t239 - t170;
t345 = -t153 - t324;
t272 = pkin(2) * t285;
t228 = t272 + t347;
t226 = t270 - t344;
t77 = t191 * t226;
t203 = t143 * t274 + t228 - t191 * (-t251 + (pkin(3) - t178) * t176) + t77;
t27 = t203 + t267;
t266 = t197 * t310;
t247 = t171 - t266;
t28 = t191 * t233 + t108 + t247;
t342 = t176 * t28 + t177 * t27;
t84 = Icges(6,6) * t177 - t176 * t231;
t151 = Icges(6,4) * t291;
t86 = -Icges(6,1) * t290 + Icges(6,5) * t177 + t151;
t40 = t182 * t86 + t183 * t84;
t218 = t231 * t191;
t337 = -Icges(6,6) * t191 + qJD(5) * t133;
t52 = t176 * t337 - t177 * t218;
t219 = t136 * t191;
t335 = -Icges(6,5) * t191 + qJD(5) * t343;
t54 = t176 * t335 - t177 * t219;
t132 = Icges(6,5) * t183 - Icges(6,6) * t182;
t82 = Icges(6,3) * t177 - t132 * t176;
t340 = qJD(5) * t40 + t182 * t52 - t183 * t54 - t191 * t82;
t41 = t182 * t87 + t183 * t85;
t51 = -t176 * t218 - t177 * t337;
t53 = -t176 * t219 - t177 * t335;
t83 = Icges(6,5) * t287 - Icges(6,6) * t288 + Icges(6,3) * t176;
t339 = qJD(5) * t41 + t182 * t51 - t183 * t53 - t191 * t83;
t131 = Icges(6,5) * t182 + Icges(6,6) * t183;
t338 = -Icges(6,3) * t191 + qJD(5) * t131;
t122 = t231 * qJD(5);
t123 = t136 * qJD(5);
t230 = t133 * t183 + t182 * t343;
t336 = qJD(5) * t230 + t122 * t182 - t123 * t183 - t131 * t191;
t302 = -t177 * t343 - t85;
t319 = -Icges(6,2) * t287 - t152 + t87;
t332 = t182 * t319 - t183 * t302;
t303 = t176 * t343 - t84;
t320 = Icges(6,2) * t290 + t151 + t86;
t331 = -t182 * t320 + t183 * t303;
t330 = -m(5) - m(6);
t329 = t110 / 0.2e1;
t328 = t111 / 0.2e1;
t322 = t177 * t82 + t84 * t291;
t321 = t176 * t82 + t86 * t287;
t307 = t182 * t84;
t306 = t183 * t86;
t229 = t182 * t133 - t183 * t343;
t96 = t131 * t176;
t48 = -t177 * t229 + t96;
t305 = t48 * t191;
t304 = rSges(5,3) + qJ(4);
t294 = t131 * t177;
t293 = t132 * t191;
t102 = t143 * t176;
t292 = t143 * t177;
t271 = pkin(2) * t283;
t261 = t137 + t162 - t170;
t260 = -pkin(3) - t317;
t259 = -t275 / 0.2e1;
t258 = t275 / 0.2e1;
t257 = -t274 / 0.2e1;
t256 = t274 / 0.2e1;
t253 = -t83 - t306;
t167 = rSges(2,1) * t197 - t196 * rSges(2,2);
t243 = rSges(2,1) * t196 + rSges(2,2) * t197;
t29 = -t290 * t86 + t322;
t30 = t177 * t83 - t290 * t87 + t85 * t291;
t238 = t30 * t176 + t29 * t177;
t31 = -t288 * t84 + t321;
t32 = t176 * t83 - t348;
t237 = t32 * t176 + t31 * t177;
t236 = -t306 + t307;
t227 = -t271 + t346;
t107 = t172 + t244;
t221 = -t170 - t263 + t280;
t105 = -t146 * t191 - t266;
t215 = t267 + t272;
t214 = t266 + t271;
t106 = t241 - t325;
t209 = -t293 * t176 - t177 * t338 + t191 * t235;
t208 = t176 * t338 - t293 * t177 + t191 * t236;
t207 = t132 * qJD(5) + t191 * t229;
t206 = t214 - t346;
t204 = t176 * t226 + t177 * t90;
t60 = -t269 + t155 + (-rSges(6,3) + t195) * t176 + t345;
t59 = t176 * t252 - t177 * t195 + t248;
t66 = -t176 * t304 + t177 * t260 + t163 - t324;
t65 = -t325 + t304 * t177 + (-pkin(3) + t240) * t176;
t47 = t176 * t229 + t294;
t46 = t47 * t191;
t11 = qJD(5) * t238 + t46;
t12 = qJD(5) * t237 + t305;
t16 = -qJD(5) * t236 + t182 * t54 + t183 * t52;
t17 = -qJD(5) * t235 + t182 * t53 + t183 * t51;
t20 = t207 * t176 - t177 * t336;
t21 = t176 * t336 + t207 * t177;
t201 = (t46 + ((t32 + t322 + t348) * t177 + (-t31 + (t253 - t307) * t177 + t30 + t321) * t176) * qJD(5)) * t259 + (-qJD(5) * t229 + t122 * t183 + t123 * t182) * t191 + (t41 + t48) * t329 + (t40 + t47) * t328 + (t12 - t305 + ((t30 + (-t83 + t307) * t177 - t321) * t177 + (t176 * t253 - t29 + t322) * t176) * qJD(5)) * t257 + (t16 + t21) * t256 + (t17 + t20 + t11) * t258 + (Icges(4,3) + Icges(3,3) + t230 + Icges(5,2) * t194 ^ 2 + (Icges(5,1) * t193 + 0.2e1 * Icges(5,4) * t194) * t193) * t189;
t200 = (-t28 * t248 - t27 * (-t223 + t345)) * t191;
t88 = t191 * t205;
t44 = t215 - t88 + t347;
t45 = t191 * t249 + t247;
t199 = (t45 * (-t268 + t325) - t44 * (-t295 - t312 - t324) + (-t260 * t44 - t304 * t45) * t177) * t191;
t113 = t191 * t241;
t80 = t191 * t254 - t266;
t79 = -t113 + t215;
t39 = qJD(5) * t204 + qJD(3);
t13 = qJDD(3) + t111 * t90 + t110 * t226 + (-t176 * t56 + t177 * t55) * qJD(5);
t6 = t176 * t339 + t209 * t177;
t5 = t176 * t340 + t208 * t177;
t4 = t209 * t176 - t177 * t339;
t3 = t208 * t176 - t177 * t340;
t1 = [Icges(2,3) * qJDD(1) + t201 + (t357 * (-t146 - t326) + t354 * (-t242 - t327) + (t105 - t117 + t266) * t104) * m(3) + ((qJDD(1) * t243 + g(3)) * t243 + (qJDD(1) * t167 + g(2)) * t167) * m(2) + (t28 * (t221 + t267) + t200 + t361 * (t60 - t326) + t360 * (t59 - t327) + (-t206 - t28 + t266 + t350) * t27) * m(6) + (t45 * (t261 + t267) + t199 + t359 * (t66 - t326) + t356 * (t65 - t327) + (-t206 - t45 - t247 + t353) * t44) * m(5) + (t80 * (t215 - t276) + t358 * (t107 - t326) + t355 * (t106 - t327) + (-t214 - t80 + t266 + t351) * t79) * m(4); t201 + (t200 + t361 * t60 + t360 * t59 + (t221 - t203) * t28 + (t227 + t350) * t27) * m(6) + (t199 + t359 * t66 + t356 * t65 + (-t228 + t88 + t261) * t45 + (t227 - t171 + t353) * t44) * m(5) + ((-t271 + t351) * t79 + t358 * t107 + t355 * t106 + (-t276 + t113) * t80) * m(4) + (-t104 * t117 + t105 * t116 + (-t104 * t191 - t357) * t146 - (t105 * t191 + t354) * t242) * m(3); m(6) * t13 + (m(4) + m(5)) * qJDD(3) + (-m(4) + t330) * g(1); t330 * (g(2) * t177 + g(3) * t176) + m(5) * (t176 * t22 + t177 * t23) + m(6) * (t176 * t7 + t177 * t8); t189 * (t41 * t176 + t40 * t177) / 0.2e1 + t191 * ((t191 * t41 + t16) * t177 + (-t191 * t40 + t17) * t176) / 0.2e1 - t11 * t289 / 0.2e1 + t177 * (t30 * t110 + t29 * t111 + t47 * t189 + t21 * t191 + (t176 * t6 + t177 * t5) * qJD(5)) / 0.2e1 + t238 * t328 + ((t191 * t30 + t5) * t177 + (-t191 * t29 + t6) * t176) * t256 + t12 * t286 / 0.2e1 + t176 * (t32 * t110 + t31 * t111 + t48 * t189 + t20 * t191 + (t176 * t4 + t177 * t3) * qJD(5)) / 0.2e1 + t237 * t329 + ((t191 * t32 + t3) * t177 + (-t191 * t31 + t4) * t176) * t258 - t191 * ((-t279 * t182 + t278 * t183) * t191 + ((t176 * t319 + t177 * t320) * t183 + (t302 * t176 + t303 * t177) * t182) * qJD(5)) / 0.2e1 + ((t96 * t274 + t293) * t177 + (t352 + (t332 * t176 + (-t331 - t294) * t177) * qJD(5)) * t176) * t257 + ((-t275 * t294 + t293) * t176 + (-t352 + (t331 * t177 + (-t332 + t96) * t176) * qJD(5)) * t177) * t259 + (t13 * t204 + t39 * ((-t56 - t78) * t176 + (t55 + t77) * t177) + t8 * t102 - t7 * t292 + (-t27 * t289 + t28 * t286) * t143 + t342 * t124 - (-t102 * t27 + t28 * t292) * t191 - (t39 * (-t102 * t176 - t177 * t292) + t342 * t144) * qJD(5) - g(1) * t144 - g(2) * t102 + g(3) * t292) * m(6);];
tau = t1;
