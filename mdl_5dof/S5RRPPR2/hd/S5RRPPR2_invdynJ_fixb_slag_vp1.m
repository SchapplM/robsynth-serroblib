% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:27
% EndTime: 2022-01-20 10:05:40
% DurationCPUTime: 6.83s
% Computational Cost: add. (12170->460), mult. (10102->604), div. (0->0), fcn. (9456->10), ass. (0->236)
t208 = qJ(1) + qJ(2);
t201 = pkin(8) + t208;
t196 = cos(t201);
t193 = t196 * pkin(3);
t195 = sin(t201);
t155 = t195 * qJ(4) + t193;
t187 = qJD(4) * t196;
t207 = qJD(1) + qJD(2);
t111 = t155 * t207 - t187;
t209 = sin(pkin(9));
t278 = qJD(5) * t207;
t118 = (qJDD(5) * t195 + t196 * t278) * t209;
t210 = cos(pkin(9));
t211 = sin(qJ(5));
t213 = cos(qJ(5));
t146 = -rSges(6,3) * t210 + (rSges(6,1) * t213 - rSges(6,2) * t211) * t209;
t161 = (-rSges(6,1) * t211 - rSges(6,2) * t213) * t209;
t150 = qJD(5) * t161;
t206 = qJDD(1) + qJDD(2);
t183 = -qJDD(5) * t210 + t206;
t184 = -qJD(5) * t210 + t207;
t212 = sin(qJ(1));
t214 = cos(qJ(1));
t215 = qJD(1) ^ 2;
t231 = (-qJDD(1) * t212 - t214 * t215) * pkin(1);
t279 = qJD(4) * t207;
t203 = cos(t208);
t205 = t207 ^ 2;
t290 = t203 * t205;
t220 = -pkin(2) * t290 + qJDD(4) * t195 + t196 * t279 + t231;
t312 = pkin(4) * t210;
t245 = pkin(7) * t209 + t312;
t136 = t245 * t195;
t189 = t196 * qJ(4);
t153 = pkin(3) * t195 - t189;
t202 = sin(t208);
t313 = pkin(2) * t202;
t262 = -t153 - t313;
t251 = -t136 + t262;
t277 = qJD(5) * t209;
t267 = t195 * t277;
t294 = t196 * t207;
t287 = t210 * t213;
t133 = t195 * t287 - t196 * t211;
t288 = t210 * t211;
t134 = t195 * t213 - t196 * t288;
t94 = -qJD(5) * t133 + t134 * t207;
t132 = t195 * t288 + t196 * t213;
t135 = t195 * t211 + t196 * t287;
t95 = -qJD(5) * t132 + t135 * t207;
t246 = rSges(6,1) * t95 + rSges(6,2) * t94;
t289 = t207 * t209;
t269 = t196 * t289;
t51 = rSges(6,3) * t269 + t246;
t286 = -t133 * rSges(6,1) + t132 * rSges(6,2);
t295 = t195 * t209;
t83 = rSges(6,3) * t295 - t286;
t16 = t150 * t267 + t118 * t146 - t183 * t83 - t184 * t51 + (-t245 * t294 - t111) * t207 + t251 * t206 + t220;
t349 = t16 - g(1);
t348 = -t146 * t267 + t184 * t83;
t74 = Icges(6,5) * t133 - Icges(6,6) * t132 + Icges(6,3) * t295;
t122 = Icges(6,4) * t133;
t77 = -Icges(6,2) * t132 + Icges(6,6) * t295 + t122;
t121 = Icges(6,4) * t132;
t81 = -Icges(6,1) * t133 - Icges(6,5) * t295 + t121;
t30 = -(t211 * t77 + t213 * t81) * t209 - t210 * t74;
t259 = t189 - t313;
t346 = t259 + t286;
t345 = -t132 * t77 - t133 * t81;
t344 = t134 * t77 - t135 * t81;
t342 = t231 - g(1);
t293 = t196 * t209;
t340 = t74 * t293;
t292 = t196 * t210;
t175 = rSges(5,1) * t292;
t327 = -t195 * rSges(5,3) - t175;
t115 = -rSges(5,2) * t293 - t327;
t197 = pkin(2) * t203;
t326 = t197 + t155;
t339 = t115 + t326;
t137 = pkin(4) * t292 + pkin(7) * t293;
t338 = t137 + t326;
t308 = pkin(1) * qJD(1);
t273 = t212 * t308;
t164 = rSges(3,1) * t202 + rSges(3,2) * t203;
t297 = t164 * t207;
t128 = -t273 - t297;
t119 = (qJDD(5) * t196 - t195 * t278) * t209;
t204 = t214 * pkin(1);
t314 = pkin(1) * t212;
t254 = qJDD(1) * t204 - t215 * t314;
t224 = t206 * t197 - t205 * t313 + t254;
t170 = qJ(4) * t294;
t186 = qJD(4) * t195;
t282 = t170 + t186;
t296 = t195 * t207;
t221 = t207 * (-pkin(3) * t296 + t282) + t206 * t155 + t195 * t279 + t224;
t270 = t195 * t289;
t92 = -qJD(5) * t135 + t132 * t207;
t93 = qJD(5) * t134 - t133 * t207;
t310 = t93 * rSges(6,1) + t92 * rSges(6,2);
t50 = -rSges(6,3) * t270 + t310;
t85 = t135 * rSges(6,1) + t134 * rSges(6,2) + rSges(6,3) * t293;
t17 = -t119 * t146 + t206 * t137 + t183 * t85 + t184 * t50 + (-t150 * t277 - qJDD(4)) * t196 - t205 * t136 + t221;
t337 = -g(2) + t17;
t163 = rSges(5,2) * t269;
t309 = rSges(5,1) * t210;
t275 = t195 * t309;
t281 = rSges(5,2) * t295 + t196 * rSges(5,3);
t114 = t275 - t281;
t253 = -t114 + t262;
t32 = (t207 * t327 - t111 + t163) * t207 + t253 * t206 + t220;
t336 = t32 - g(1);
t283 = rSges(5,2) * t270 + rSges(5,3) * t294;
t33 = -qJDD(4) * t196 + t206 * t115 + t207 * (-t207 * t275 + t283) + t221;
t335 = t33 - g(2);
t154 = rSges(4,1) * t195 + rSges(4,2) * t196;
t171 = rSges(4,2) * t296;
t334 = -t206 * t154 - t207 * (rSges(4,1) * t294 - t171) + (-t202 * t206 - t290) * pkin(2) + t342;
t192 = t196 * rSges(4,1);
t156 = -rSges(4,2) * t195 + t192;
t333 = -t154 * t205 + t206 * t156 - g(2) + t224;
t194 = t203 * rSges(3,1);
t291 = t202 * t207;
t145 = -rSges(3,2) * t291 + t194 * t207;
t332 = -t145 * t207 - t164 * t206 + t342;
t165 = -rSges(3,2) * t202 + t194;
t331 = t165 * t206 - t207 * t297 - g(2) + t254;
t26 = t340 + t344;
t330 = t26 - t340;
t329 = t207 * t114 + t283;
t328 = t156 + t197;
t229 = -pkin(2) * t291 - t273;
t249 = t186 - t273;
t280 = -t207 * t153 + t186;
t325 = t170 + t249 - t229 - t280;
t324 = t282 - t280;
t76 = Icges(6,5) * t135 + Icges(6,6) * t134 + Icges(6,3) * t293;
t300 = Icges(6,4) * t135;
t79 = Icges(6,2) * t134 + Icges(6,6) * t293 + t300;
t123 = Icges(6,4) * t134;
t82 = Icges(6,1) * t135 + Icges(6,5) * t293 + t123;
t25 = -t132 * t79 + t133 * t82 + t76 * t295;
t230 = -t312 - pkin(3) + (-rSges(6,3) - pkin(7)) * t209;
t272 = t214 * t308;
t322 = -t196 * t146 * t277 + t184 * t85 + t207 * t338 - t187;
t37 = t272 + t322;
t303 = t37 * t202;
t36 = t207 * t251 + t249 - t348;
t323 = ((-t36 * t203 - t303) * pkin(2) + t36 * t230 * t196 + (-t36 * qJ(4) + t37 * (-rSges(6,3) * t209 - pkin(3) - t245)) * t195) * t207 + t349 * t195 * t230;
t321 = t207 * t136 + t310 + t348;
t320 = t195 * (-Icges(6,2) * t133 - t121 - t81) + t196 * (-Icges(6,2) * t135 + t123 + t82);
t319 = -m(5) - m(6);
t318 = t118 / 0.2e1;
t317 = t119 / 0.2e1;
t316 = t195 / 0.2e1;
t315 = -t196 / 0.2e1;
t298 = Icges(6,4) * t213;
t140 = -Icges(6,6) * t210 + (-Icges(6,2) * t211 + t298) * t209;
t299 = Icges(6,4) * t211;
t141 = -Icges(6,5) * t210 + (Icges(6,1) * t213 - t299) * t209;
t158 = (-Icges(6,5) * t211 - Icges(6,6) * t213) * t209;
t147 = qJD(5) * t158;
t159 = (-Icges(6,2) * t213 - t299) * t209;
t148 = qJD(5) * t159;
t160 = (-Icges(6,1) * t211 - t298) * t209;
t149 = qJD(5) * t160;
t42 = -t147 * t210 + (-t148 * t211 + t149 * t213 + (-t140 * t213 - t141 * t211) * qJD(5)) * t209;
t139 = -Icges(6,3) * t210 + (Icges(6,5) * t213 - Icges(6,6) * t211) * t209;
t69 = -t139 * t210 + (-t140 * t211 + t141 * t213) * t209;
t311 = t69 * t183 + t42 * t184;
t52 = -t132 * t140 + t133 * t141 + t139 * t295;
t307 = t184 * t52;
t24 = t295 * t74 + t345;
t306 = t195 * t24;
t305 = t30 * t118;
t31 = -t210 * t76 + (-t211 * t79 + t213 * t82) * t209;
t304 = t31 * t119;
t285 = -t140 + t160;
t284 = t141 + t159;
t248 = -t187 + t272;
t68 = t207 * t339 + t248;
t276 = t68 * t313;
t27 = t134 * t79 + t135 * t82 + t76 * t293;
t265 = -pkin(3) - t309;
t264 = -t277 / 0.2e1;
t263 = t277 / 0.2e1;
t130 = -t154 - t313;
t258 = t195 * t264;
t257 = t195 * t263;
t256 = t196 * t264;
t255 = t196 * t263;
t181 = rSges(2,1) * t214 - rSges(2,2) * t212;
t180 = rSges(2,1) * t212 + rSges(2,2) * t214;
t243 = t36 * t195 - t196 * t37;
t242 = -t195 * t50 + t196 * t51;
t241 = -t195 * t85 + t196 * t83;
t240 = t195 * (-Icges(6,5) * t132 - Icges(6,6) * t133) + t196 * (Icges(6,5) * t134 - Icges(6,6) * t135);
t237 = t187 - t246;
t234 = t209 * t240;
t233 = (t196 * t25 + t306) * t209;
t232 = (t195 * t26 + t196 * t27) * t209;
t227 = (Icges(6,1) * t134 - t300 - t79) * t196 + (-Icges(6,1) * t132 - t122 - t77) * t195;
t55 = t85 + t338;
t96 = t195 * t265 + t259 + t281;
t53 = t134 * t140 + t135 * t141 + t139 * t293;
t43 = t53 * t184;
t10 = qJD(5) * t232 + t43;
t45 = Icges(6,5) * t95 + Icges(6,6) * t94 + Icges(6,3) * t269;
t47 = Icges(6,4) * t95 + Icges(6,2) * t94 + Icges(6,6) * t269;
t49 = Icges(6,1) * t95 + Icges(6,4) * t94 + Icges(6,5) * t269;
t13 = -t210 * t45 + (-t211 * t47 + t213 * t49 + (t211 * t81 - t213 * t77) * qJD(5)) * t209;
t44 = Icges(6,5) * t93 + Icges(6,6) * t92 - Icges(6,3) * t270;
t46 = Icges(6,4) * t93 + Icges(6,2) * t92 - Icges(6,6) * t270;
t48 = Icges(6,1) * t93 + Icges(6,4) * t92 - Icges(6,5) * t270;
t14 = -t210 * t44 + (-t211 * t46 + t213 * t48 + (-t211 * t82 - t213 * t79) * qJD(5)) * t209;
t21 = t134 * t148 + t135 * t149 + t140 * t92 + t141 * t93 + (-t139 * t296 + t147 * t196) * t209;
t22 = -t132 * t148 + t133 * t149 + t140 * t94 + t141 * t95 + (t139 * t294 + t147 * t195) * t209;
t9 = qJD(5) * t233 + t307;
t219 = (t43 + ((t24 + t27 - t345) * t196 + t330 * t195) * t277) * t258 + t305 / 0.2e1 + t304 / 0.2e1 + t52 * t318 + t53 * t317 + t311 + (t14 + t21) * t255 + (t13 + t22 + t10) * t257 + (-t307 + ((-t25 + t330 - t344) * t196 - t306) * t277 + t9) * t256 + (Icges(5,2) * t210 ^ 2 + (Icges(5,1) * t209 + 0.2e1 * Icges(5,4) * t210) * t209 + Icges(4,3) + Icges(3,3)) * t206;
t106 = t130 * t207 - t273;
t107 = t207 * t328 + t272;
t218 = (t106 * (-t192 - t197) + t107 * t130) * t207;
t67 = t207 * t253 + t249;
t217 = (t67 * (-t175 - t193 - t197) - t276 + (t67 * (-rSges(5,3) - qJ(4)) + t68 * t265) * t195) * t207;
t143 = t207 * t154;
t129 = t165 * t207 + t272;
t105 = rSges(6,1) * t134 - rSges(6,2) * t135;
t104 = -rSges(6,1) * t132 - rSges(6,2) * t133;
t40 = t241 * t277 + qJD(3);
t15 = -t118 * t85 + t119 * t83 + t242 * t277 + qJDD(3);
t6 = -t132 * t46 + t133 * t48 + t79 * t94 + t82 * t95 + (t195 * t44 + t294 * t76) * t209;
t5 = -t132 * t47 + t133 * t49 + t77 * t94 - t81 * t95 + (t195 * t45 + t294 * t74) * t209;
t4 = t134 * t46 + t135 * t48 + t79 * t92 + t82 * t93 + (t196 * t44 - t296 * t76) * t209;
t3 = t134 * t47 + t135 * t49 + t77 * t92 - t81 * t93 + (t196 * t45 - t296 * t74) * t209;
t1 = [Icges(2,3) * qJDD(1) + t219 + (t331 * (t165 + t204) + t332 * (-t164 - t314) + (-t145 - t272 + t129) * t128) * m(3) + ((t180 ^ 2 + t181 ^ 2) * qJDD(1) + g(1) * t180 - g(2) * t181) * m(2) + (t36 * (t237 - t272) + t337 * (t204 + t55) + (t36 + t321 + t325) * t37 + t323 + t349 * (t346 - t314)) * m(6) + (t67 * (t163 - t248) + t217 + t335 * (t204 + t339) + t336 * (t96 - t314) + (t67 + t325 + t329) * t68) * m(5) + (t106 * (t171 - t272) + t218 + t333 * (t328 + t204) + t334 * (t130 - t314) + (-t273 + t106 + t143 - t229) * t107) * m(4); t219 + (pkin(2) * t303 * t207 + t337 * t55 + (t321 + t324) * t37 + (t237 + t322) * t36 + t323 + t349 * t346) * m(6) + (t67 * t163 + t217 - (-t339 * t67 - t276) * t207 + t335 * t339 + t336 * t96 + (t324 + t329) * t68) * m(5) + (t106 * t171 + t218 + t107 * t143 - (-t106 * t328 - t107 * t313) * t207 + t333 * t328 + t334 * t130) * m(4) + (-t128 * t145 - t129 * t297 + (t128 * t207 + t331) * t165 + (t129 * t207 - t332) * t164) * m(3); m(6) * t15 + (m(4) + m(5)) * qJDD(3) + (-m(4) + t319) * g(3); t319 * (g(1) * t195 - g(2) * t196) + 0.2e1 * (t16 * t316 + t17 * t315) * m(6) + 0.2e1 * (t315 * t33 + t316 * t32) * m(5); -t10 * t270 / 0.2e1 + (-t210 * t53 + t232) * t317 + (-t21 * t210 + ((t207 * t26 + t4) * t196 + (-t207 * t27 + t3) * t195) * t209) * t255 + (t118 * t24 + t119 * t25 + t183 * t52 + t184 * t22 + (t195 * t5 + t196 * t6) * t277) * t295 / 0.2e1 + (-t210 * t52 + t233) * t318 + (-t210 * t22 + ((t207 * t24 + t6) * t196 + (-t207 * t25 + t5) * t195) * t209) * t257 - t210 * (t305 + t304 + (t13 * t195 + t14 * t196) * t277 + t311) / 0.2e1 + t183 * (-t210 * t69 + (t195 * t30 + t196 * t31) * t209) / 0.2e1 + t184 * (-t210 * t42 + ((t207 * t30 + t14) * t196 + (-t207 * t31 + t13) * t195) * t209) / 0.2e1 + ((t134 * t284 + t135 * t285 + t158 * t293) * t184 + (t134 * t320 + t227 * t135 + t196 * t234) * t277) * t256 + ((-t132 * t284 + t133 * t285 + t158 * t295) * t184 + (-t132 * t320 + t133 * t227 + t195 * t234) * t277) * t258 - t184 * (-t210 * t158 * t184 + ((-t211 * t284 + t213 * t285) * t184 + ((-t211 * t320 + t213 * t227) * t209 - t240 * t210) * qJD(5)) * t209) / 0.2e1 + (t118 * t26 + t119 * t27 + t183 * t53 + t184 * t21 + (t195 * t3 + t196 * t4) * t277 + t207 * t9) * t293 / 0.2e1 + ((t16 * t83 - t17 * t85 + t36 * t51 - t37 * t50) * t210 + (t15 * t241 + t40 * (-t294 * t85 - t296 * t83 + t242) + t243 * t150 + ((t207 * t36 - t17) * t196 + (t207 * t37 + t16) * t195) * t146) * t209 - (-t104 * t36 + t105 * t37) * t184 - (t40 * (t104 * t196 - t105 * t195) + t243 * t161) * t277 - g(1) * t105 - g(2) * t104 - g(3) * t161) * m(6);];
tau = t1;
