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
% m [6x1]
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:51:19
% EndTime: 2022-01-20 09:51:27
% DurationCPUTime: 3.99s
% Computational Cost: add. (9636->411), mult. (5845->490), div. (0->0), fcn. (4304->10), ass. (0->249)
t204 = sin(qJ(1));
t205 = cos(qJ(1));
t206 = qJD(1) ^ 2;
t225 = (-qJDD(1) * t204 - t205 * t206) * pkin(1);
t349 = t225 - g(1);
t200 = qJ(1) + qJ(2);
t192 = pkin(8) + t200;
t183 = sin(t192);
t198 = pkin(9) + qJ(5);
t191 = cos(t198);
t294 = t183 * t191;
t272 = rSges(6,1) * t294;
t193 = sin(t200);
t319 = pkin(2) * t193;
t184 = cos(t192);
t202 = cos(pkin(9));
t185 = pkin(4) * t202 + pkin(3);
t203 = -pkin(7) - qJ(4);
t335 = -t183 * t185 - t184 * t203;
t348 = -t272 - t319 + t335;
t310 = pkin(1) * qJD(1);
t269 = t204 * t310;
t194 = cos(t200);
t142 = rSges(3,1) * t193 + rSges(3,2) * t194;
t199 = qJD(1) + qJD(2);
t296 = t142 * t199;
t105 = -t269 - t296;
t275 = qJD(5) * t199;
t112 = -qJDD(5) * t184 + t183 * t275;
t190 = sin(t198);
t141 = rSges(6,1) * t191 - rSges(6,2) * t190;
t122 = t141 * qJD(5);
t311 = rSges(6,2) * t191;
t140 = rSges(6,1) * t190 + t311;
t292 = t183 * t203;
t144 = t199 * t292;
t149 = t184 * t185;
t197 = qJDD(1) + qJDD(2);
t279 = qJD(4) * t199;
t196 = t199 ^ 2;
t287 = t194 * t196;
t215 = -pkin(2) * t287 + qJDD(4) * t183 + t184 * t279 + t225;
t171 = t184 * qJ(4);
t318 = pkin(3) * t183;
t123 = -t171 + t318;
t258 = -t123 - t319;
t76 = t123 + t335;
t295 = t183 * t190;
t150 = rSges(6,2) * t295;
t334 = -t184 * rSges(6,3) - t150;
t89 = t272 + t334;
t238 = t258 + t76 - t89;
t277 = qJD(5) * t184;
t180 = t184 * pkin(3);
t332 = t183 * qJ(4) + t180;
t290 = t184 * t191;
t230 = rSges(6,1) * t290 + t183 * rSges(6,3);
t267 = qJD(5) * t311;
t291 = t184 * t190;
t270 = rSges(6,2) * t291;
t276 = qJD(5) * t190;
t266 = -t199 * t270 + (-rSges(6,1) * t276 - t267) * t183;
t58 = t199 * t230 + t266;
t169 = qJD(4) * t184;
t92 = t199 * t332 - t169;
t7 = -t122 * t277 + t112 * t140 + t238 * t197 + (-t92 + t144 - t58 + (t332 - t149) * t199) * t199 + t215;
t347 = -g(1) + t7;
t111 = qJDD(5) * t183 + t184 * t275;
t289 = t184 * t199;
t156 = qJ(4) * t289;
t186 = pkin(2) * t194;
t195 = t205 * pkin(1);
t320 = pkin(1) * t204;
t250 = qJDD(1) * t195 - t206 * t320;
t222 = t197 * t186 - t196 * t319 + t250;
t168 = qJD(4) * t183;
t282 = t156 + t168;
t293 = t183 * t199;
t214 = -qJDD(4) * t184 + t197 * t332 + t183 * t279 + t222 + t199 * (-pkin(3) * t293 + t282);
t278 = qJD(5) * t183;
t90 = -t270 + t230;
t339 = t149 - t292 + t90;
t315 = -t332 + t339;
t229 = rSges(6,3) * t289 + t199 * t150 - t184 * t267;
t264 = t184 * t276;
t57 = (-t191 * t293 - t264) * rSges(6,1) + t229;
t8 = -t122 * t278 - t111 * t140 + t315 * t197 + (-t156 + t57 + (t335 + t318) * t199) * t199 + t214;
t346 = -g(2) + t8;
t201 = sin(pkin(9));
t312 = rSges(5,2) * t201;
t271 = t184 * t312;
t137 = t199 * t271;
t313 = rSges(5,1) * t202;
t273 = t183 * t313;
t160 = t183 * t312;
t281 = t184 * rSges(5,3) + t160;
t93 = t273 - t281;
t252 = t258 - t93;
t161 = t184 * t313;
t333 = -t183 * rSges(5,3) - t161;
t22 = (t333 * t199 + t137 - t92) * t199 + t252 * t197 + t215;
t345 = -g(1) + t22;
t284 = rSges(5,3) * t289 + t199 * t160;
t94 = -t271 - t333;
t23 = t197 * t94 + t199 * (-t199 * t273 + t284) + t214;
t344 = -g(2) + t23;
t124 = rSges(4,1) * t183 + rSges(4,2) * t184;
t157 = rSges(4,2) * t293;
t343 = -t197 * t124 - t199 * (rSges(4,1) * t289 - t157) + (-t193 * t197 - t287) * pkin(2) + t349;
t179 = t184 * rSges(4,1);
t126 = -rSges(4,2) * t183 + t179;
t342 = -t124 * t196 + t197 * t126 - g(2) + t222;
t182 = t194 * rSges(3,1);
t288 = t193 * t199;
t116 = -rSges(3,2) * t288 + t182 * t199;
t341 = -t116 * t199 - t142 * t197 + t349;
t143 = -rSges(3,2) * t193 + t182;
t340 = t143 * t197 - t199 * t296 - g(2) + t250;
t257 = t332 + t186;
t338 = t257 + t94;
t337 = t199 * t93 + t284;
t336 = t126 + t186;
t181 = Icges(6,4) * t191;
t235 = -Icges(6,2) * t190 + t181;
t134 = Icges(6,1) * t190 + t181;
t78 = t199 * t89;
t331 = -rSges(6,1) * t264 - t199 * t76 + t168 + t229 + t78;
t131 = Icges(6,5) * t191 - Icges(6,6) * t190;
t130 = Icges(6,5) * t190 + Icges(6,6) * t191;
t217 = Icges(6,3) * t199 - qJD(5) * t130;
t227 = t235 * t184;
t85 = Icges(6,6) * t183 + t227;
t307 = t190 * t85;
t301 = Icges(6,4) * t190;
t135 = Icges(6,1) * t191 - t301;
t228 = t135 * t184;
t87 = Icges(6,5) * t183 + t228;
t240 = -t191 * t87 + t307;
t330 = -t131 * t293 + t184 * t217 + t199 * t240;
t226 = t131 * t184;
t84 = Icges(6,4) * t294 - Icges(6,2) * t295 - Icges(6,6) * t184;
t308 = t190 * t84;
t147 = Icges(6,4) * t295;
t86 = Icges(6,1) * t294 - Icges(6,5) * t184 - t147;
t241 = -t191 * t86 + t308;
t329 = t183 * t217 + (t226 + t241) * t199;
t328 = -t140 * t278 + t199 * (t257 + t315);
t132 = Icges(6,2) * t191 + t301;
t233 = t132 * t190 - t134 * t191;
t327 = t131 * qJD(5) + t199 * t233;
t82 = Icges(6,5) * t294 - Icges(6,6) * t295 - Icges(6,3) * t184;
t30 = -t183 * t241 - t184 * t82;
t304 = -t134 * t183 - t84;
t314 = -Icges(6,2) * t294 - t147 + t86;
t326 = -t190 * t314 + t191 * t304;
t325 = -m(5) - m(6);
t324 = t111 / 0.2e1;
t323 = t112 / 0.2e1;
t322 = t183 / 0.2e1;
t321 = -t184 / 0.2e1;
t317 = -t183 * t82 - t86 * t290;
t83 = Icges(6,3) * t183 + t226;
t316 = t183 * t83 + t87 * t290;
t298 = t130 * t184;
t49 = -t183 * t233 - t298;
t306 = t49 * t199;
t305 = -t132 * t184 + t87;
t303 = -t134 * t184 - t85;
t299 = t130 * t183;
t297 = t131 * t199;
t286 = -t132 + t135;
t285 = t134 + t235;
t113 = t199 * t123;
t280 = t168 - t113;
t268 = t205 * t310;
t248 = -t169 + t268;
t47 = t199 * t338 + t248;
t274 = t47 * t319;
t265 = t140 * t277;
t263 = -pkin(3) - t313;
t262 = -t278 / 0.2e1;
t261 = t278 / 0.2e1;
t260 = -t277 / 0.2e1;
t259 = t277 / 0.2e1;
t107 = -t124 - t319;
t69 = t87 * t294;
t255 = t184 * t83 - t69;
t254 = -t82 + t307;
t249 = t168 - t269;
t247 = t168 - t265;
t166 = rSges(2,1) * t205 - rSges(2,2) * t204;
t165 = rSges(2,1) * t204 + rSges(2,2) * t205;
t28 = t199 * t238 + t247 - t269;
t29 = t248 + t328;
t245 = -t183 * t29 - t184 * t28;
t31 = -t295 * t85 - t255;
t244 = t31 * t183 - t30 * t184;
t32 = -t291 * t84 - t317;
t33 = -t291 * t85 + t316;
t243 = t33 * t183 - t32 * t184;
t242 = t183 * t89 + t184 * t90;
t42 = t190 * t86 + t191 * t84;
t43 = t190 * t87 + t191 * t85;
t239 = t144 + t169 - t266;
t234 = t132 * t191 + t134 * t190;
t224 = -pkin(2) * t288 - t269;
t223 = -t190 * t305 + t191 * t303;
t221 = t224 + t280;
t61 = -t334 + t348;
t220 = (-t190 * t285 + t191 * t286) * t199;
t219 = Icges(6,5) * t199 - qJD(5) * t134;
t218 = Icges(6,6) * t199 - qJD(5) * t132;
t62 = t186 + t339;
t67 = t183 * t263 + t171 + t281 - t319;
t53 = t184 * t218 - t235 * t293;
t55 = -t135 * t293 + t184 * t219;
t213 = -qJD(5) * t43 - t190 * t53 + t191 * t55 + t199 * t83;
t54 = t183 * t218 + t199 * t227;
t56 = t183 * t219 + t199 * t228;
t212 = -qJD(5) * t42 - t190 * t54 + t191 * t56 + t199 * t82;
t120 = t235 * qJD(5);
t121 = t135 * qJD(5);
t211 = -qJD(5) * t234 - t120 * t190 + t121 * t191 + t130 * t199;
t11 = qJD(5) * t244 + t306;
t50 = -t184 * t233 + t299;
t48 = t50 * t199;
t12 = qJD(5) * t243 + t48;
t16 = -qJD(5) * t241 + t190 * t56 + t191 * t54;
t17 = -qJD(5) * t240 + t190 * t55 + t191 * t53;
t20 = t327 * t183 + t211 * t184;
t21 = t211 * t183 - t327 * t184;
t210 = (t48 + ((t31 - t69 + (t83 + t308) * t184 + t317) * t184 + t316 * t183) * qJD(5)) * t259 + (-qJD(5) * t233 + t120 * t191 + t121 * t190) * t199 + (t43 + t50) * t324 + (t42 + t49) * t323 + (-t306 + ((t184 * t254 - t316 + t33) * t184 + (t183 * t254 + t255 + t32) * t183) * qJD(5) + t11) * t262 + (t17 + t20) * t261 + (t12 + t16 + t21) * t260 + (Icges(4,3) + Icges(3,3) + t234 + Icges(5,2) * t202 ^ 2 + (Icges(5,1) * t201 + 0.2e1 * Icges(5,4) * t202) * t201) * t197;
t79 = t107 * t199 - t269;
t80 = t199 * t336 + t268;
t209 = (t79 * (-t179 - t186) + t80 * t107) * t199;
t46 = t199 * t252 + t249;
t208 = (t46 * (-t161 - t180 - t186) - t274 + (t46 * (-rSges(5,3) - qJ(4)) + t47 * t263) * t183) * t199;
t207 = (t28 * (-t149 - t230 - t186) + t29 * t348) * t199;
t114 = t199 * t124;
t106 = t143 * t199 + t268;
t104 = t140 * t184;
t103 = t140 * t183;
t41 = qJD(5) * t242 + qJD(3);
t13 = t111 * t89 - t112 * t90 + qJDD(3) + (t183 * t58 + t184 * t57) * qJD(5);
t6 = t213 * t183 - t330 * t184;
t5 = t212 * t183 - t329 * t184;
t4 = t330 * t183 + t213 * t184;
t3 = t329 * t183 + t212 * t184;
t1 = [Icges(2,3) * qJDD(1) + t210 + (t340 * (t143 + t195) + t341 * (-t142 - t320) + (-t116 - t268 + t106) * t105) * m(3) + ((t165 ^ 2 + t166 ^ 2) * qJDD(1) + g(1) * t165 - g(2) * t166) * m(2) + (t28 * (t239 - t268) + t207 + t346 * (t195 + t62) + t347 * (t61 - t320) + (-t221 + t28 + t265 - t269 + t331) * t29) * m(6) + (t46 * (t137 - t248) + t208 + t344 * (t195 + t338) + t345 * (t67 - t320) + (-t221 + t46 + t156 + t249 + t337) * t47) * m(5) + (t79 * (t157 - t268) + t209 + t342 * (t336 + t195) + t343 * (t107 - t320) + (t114 - t224 + t79 - t269) * t80) * m(4); t210 + (t207 + t346 * t62 + t347 * t61 + (t199 * t319 + t113 - t247 + t331) * t29 + (-t169 + t239 + t328) * t28) * m(6) + (t46 * t137 + t208 - (-t338 * t46 - t274) * t199 + t344 * t338 + t345 * t67 + (t282 - t280 + t337) * t47) * m(5) + (t79 * t157 + t209 + t80 * t114 - (-t319 * t80 - t336 * t79) * t199 + t342 * t336 + t343 * t107) * m(4) + (-t105 * t116 - t106 * t296 + (t105 * t199 + t340) * t143 + (t106 * t199 - t341) * t142) * m(3); m(6) * t13 + (m(4) + m(5)) * qJDD(3) + (-m(4) + t325) * g(3); t325 * (g(1) * t183 - g(2) * t184) + 0.2e1 * (t321 * t8 + t322 * t7) * m(6) + 0.2e1 * (t22 * t322 + t23 * t321) * m(5); t12 * t289 / 0.2e1 + (t111 * t33 + t32 * t112 + t50 * t197 + t20 * t199 + (t183 * t4 - t184 * t3) * qJD(5)) * t322 + t243 * t324 + ((t199 * t33 - t3) * t184 + (t199 * t32 + t4) * t183) * t261 + t11 * t293 / 0.2e1 + (t31 * t111 + t30 * t112 + t49 * t197 + t21 * t199 + (t183 * t6 - t184 * t5) * qJD(5)) * t321 + t244 * t323 + ((t199 * t31 - t5) * t184 + (t199 * t30 + t6) * t183) * t260 + t197 * (t183 * t43 - t184 * t42) / 0.2e1 + t199 * ((t199 * t43 - t16) * t184 + (t199 * t42 + t17) * t183) / 0.2e1 + ((-t278 * t298 + t297) * t183 + (t220 + (-t326 * t184 + (t299 + t223) * t183) * qJD(5)) * t184) * t262 + ((-t277 * t299 - t297) * t184 + (t220 + (t223 * t183 + (-t326 + t298) * t184) * qJD(5)) * t183) * t259 - t199 * ((t286 * t190 + t285 * t191) * t199 + ((t183 * t305 - t184 * t314) * t191 + (t183 * t303 - t304 * t184) * t190) * qJD(5)) / 0.2e1 + (t13 * t242 + t41 * ((t57 + t78) * t184 + (-t199 * t90 + t58) * t183) + t245 * t122 + ((-t199 * t29 - t7) * t184 + (t199 * t28 - t8) * t183) * t140 - (t103 * t28 - t104 * t29) * t199 - (t41 * (-t103 * t183 - t104 * t184) + t245 * t141) * qJD(5) + g(1) * t104 + g(2) * t103 - g(3) * t141) * m(6);];
tau = t1;
