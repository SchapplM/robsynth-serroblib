% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:59
% EndTime: 2019-12-05 15:03:15
% DurationCPUTime: 11.95s
% Computational Cost: add. (7705->460), mult. (11388->723), div. (0->0), fcn. (10691->6), ass. (0->241)
t194 = pkin(8) + qJ(3);
t190 = sin(t194);
t191 = cos(t194);
t362 = (Icges(5,5) - Icges(4,6)) * t191 + (Icges(5,4) - Icges(4,5)) * t190;
t195 = sin(pkin(7));
t192 = t195 ^ 2;
t196 = cos(pkin(7));
t193 = t196 ^ 2;
t299 = t192 + t193;
t268 = qJD(3) * t299;
t287 = qJD(5) * t191;
t295 = qJD(3) * t195;
t168 = t196 * t287 + t295;
t197 = sin(qJ(5));
t198 = cos(qJ(5));
t263 = rSges(6,1) * t197 + rSges(6,2) * t198;
t217 = -rSges(6,3) * t190 + t191 * t263;
t365 = t168 * t217;
t293 = qJD(3) * t196;
t169 = t195 * t287 - t293;
t364 = t169 * t217;
t363 = t362 * qJD(3);
t296 = qJD(3) * t191;
t297 = qJD(3) * t190;
t316 = Icges(5,6) * t191;
t322 = Icges(4,4) * t191;
t323 = Icges(4,4) * t190;
t361 = (-Icges(5,6) * t190 - t323 + (-Icges(4,2) - Icges(5,3)) * t191) * t297 + (t316 + t322 + (Icges(4,1) + Icges(5,2)) * t190) * t296;
t174 = rSges(4,1) * t190 + rSges(4,2) * t191;
t360 = t174 * t268;
t359 = t363 * t195;
t358 = t363 * t196;
t357 = t362 * t195;
t356 = t362 * t196;
t240 = Icges(5,3) * t190 - t316;
t108 = Icges(5,5) * t195 + t196 * t240;
t109 = -Icges(5,5) * t196 + t195 * t240;
t312 = t190 * t196;
t186 = Icges(5,6) * t312;
t310 = t191 * t196;
t110 = Icges(5,4) * t195 - Icges(5,2) * t310 + t186;
t313 = t190 * t195;
t185 = Icges(5,6) * t313;
t311 = t191 * t195;
t111 = -Icges(5,4) * t196 - Icges(5,2) * t311 + t185;
t247 = -Icges(4,2) * t190 + t322;
t250 = Icges(4,1) * t191 - t323;
t349 = -(Icges(4,6) * t195 + t196 * t247) * t191 - (Icges(4,5) * t195 + t196 * t250) * t190;
t350 = (-Icges(4,6) * t196 + t195 * t247) * t191 + (-Icges(4,5) * t196 + t195 * t250) * t190;
t355 = (t108 * t195 - t109 * t196) * t191 + ((Icges(5,3) * t310 + t110 + t186) * t195 - (Icges(5,3) * t311 + t111 + t185) * t196) * t190 + t195 * t349 + t196 * t350;
t354 = t361 * t196 + (-t108 * t191 - t110 * t190 - t349) * qJD(3);
t353 = t361 * t195 + (-t109 * t191 - t111 * t190 + t350) * qJD(3);
t262 = rSges(5,2) * t190 + rSges(5,3) * t191;
t348 = t262 * t268;
t340 = t168 / 0.2e1;
t338 = t169 / 0.2e1;
t242 = Icges(6,5) * t197 + Icges(6,6) * t198;
t214 = -Icges(6,3) * t190 + t191 * t242;
t319 = Icges(6,4) * t197;
t244 = Icges(6,2) * t198 + t319;
t215 = -Icges(6,6) * t190 + t191 * t244;
t318 = Icges(6,4) * t198;
t248 = Icges(6,1) * t197 + t318;
t216 = -Icges(6,5) * t190 + t191 * t248;
t306 = t196 * t198;
t309 = t195 * t197;
t162 = t190 * t306 - t309;
t154 = Icges(6,4) * t162;
t307 = t196 * t197;
t308 = t195 * t198;
t164 = t190 * t308 + t307;
t155 = Icges(6,4) * t164;
t159 = (Icges(6,2) * t197 - t318) * t191;
t163 = t190 * t307 + t308;
t165 = -t190 * t309 + t306;
t288 = qJD(5) * t190;
t66 = Icges(6,1) * t163 + Icges(6,5) * t310 + t154;
t67 = -Icges(6,1) * t165 + Icges(6,5) * t311 + t155;
t207 = t168 * (-Icges(6,2) * t163 + t154 + t66) + t169 * (Icges(6,2) * t165 + t155 + t67) + t288 * (-t216 + t159);
t160 = (-Icges(6,1) * t198 + t319) * t191;
t320 = Icges(6,4) * t165;
t321 = Icges(6,4) * t163;
t64 = Icges(6,2) * t162 + Icges(6,6) * t310 + t321;
t65 = Icges(6,2) * t164 + Icges(6,6) * t311 - t320;
t208 = t168 * (-Icges(6,1) * t162 + t321 + t64) + t169 * (-Icges(6,1) * t164 - t320 + t65) - t288 * (t160 + t215);
t158 = (-Icges(6,5) * t198 + Icges(6,6) * t197) * t191;
t96 = Icges(6,3) * t191 + t190 * t242;
t58 = qJD(3) * t96 + qJD(5) * t158;
t231 = t191 * t58 + t214 * t297;
t62 = Icges(6,5) * t163 + Icges(6,6) * t162 + Icges(6,3) * t310;
t19 = t162 * t64 + t163 * t66 + t310 * t62;
t63 = -Icges(6,5) * t165 + Icges(6,6) * t164 + Icges(6,3) * t311;
t20 = t162 * t65 + t163 * t67 + t310 * t63;
t328 = t195 * t20;
t259 = t19 * t196 + t328;
t38 = -t162 * t215 - t163 * t216 - t214 * t310;
t326 = t38 * t191;
t98 = Icges(6,6) * t191 + t190 * t244;
t59 = qJD(3) * t98 + qJD(5) * t159;
t100 = Icges(6,5) * t191 + t190 * t248;
t60 = qJD(3) * t100 + qJD(5) * t160;
t279 = t197 * t296;
t90 = qJD(5) * t162 + t196 * t279;
t278 = t198 * t296;
t91 = -qJD(5) * t163 + t196 * t278;
t206 = (t162 * t59 + t163 * t60 + t196 * t231 - t215 * t91 - t216 * t90) * t190 + (-t190 * t259 + t326) * qJD(3);
t280 = t190 * t293;
t42 = Icges(6,5) * t90 + Icges(6,6) * t91 - Icges(6,3) * t280;
t233 = t191 * t42 - t297 * t62;
t44 = Icges(6,4) * t90 + Icges(6,2) * t91 - Icges(6,6) * t280;
t46 = Icges(6,1) * t90 + Icges(6,4) * t91 - Icges(6,5) * t280;
t7 = t162 * t44 + t163 * t46 + t196 * t233 + t64 * t91 + t66 * t90;
t281 = t190 * t295;
t92 = qJD(5) * t164 + t195 * t279;
t93 = qJD(5) * t165 + t195 * t278;
t43 = Icges(6,5) * t92 + Icges(6,6) * t93 - Icges(6,3) * t281;
t232 = t191 * t43 - t297 * t63;
t45 = Icges(6,4) * t92 + Icges(6,2) * t93 - Icges(6,6) * t281;
t47 = Icges(6,1) * t92 + Icges(6,4) * t93 - Icges(6,5) * t281;
t8 = t162 * t45 + t163 * t47 + t196 * t232 + t65 * t91 + t67 * t90;
t343 = qJD(5) * t206 / 0.2e1 + t7 * t340 + t8 * t338;
t342 = -t190 * t268 / 0.2e1;
t341 = -t168 / 0.2e1;
t339 = -t169 / 0.2e1;
t337 = -t191 / 0.2e1;
t336 = pkin(6) * qJD(3) ^ 2;
t290 = qJD(4) * t195;
t181 = t190 * t290;
t172 = pkin(3) * t190 - qJ(4) * t191;
t226 = qJD(3) * t172;
t94 = -t195 * t226 + t181;
t289 = qJD(4) * t196;
t183 = t190 * t289;
t95 = -t196 * t226 + t183;
t331 = t195 * t94 + t196 * t95;
t329 = t190 * t214;
t21 = t164 * t64 - t165 * t66 + t311 * t62;
t327 = t196 * t21;
t39 = -t164 * t215 + t165 * t216 - t214 * t311;
t325 = t39 * t191;
t175 = pkin(3) * t191 + qJ(4) * t190;
t291 = qJD(4) * t191;
t122 = qJD(3) * t175 - t291;
t176 = -rSges(5,2) * t191 + rSges(5,3) * t190;
t304 = -t176 * qJD(3) - t122;
t156 = t175 * t195;
t157 = t175 * t196;
t303 = t195 * t156 + t196 * t157;
t302 = -t172 + t262;
t301 = -t175 - t176;
t189 = qJD(2) * t195;
t300 = t183 + t189;
t298 = qJD(2) * t196;
t286 = qJD(3) * qJD(4);
t285 = t191 * t336;
t284 = t190 * t286 + t95 * t293 + t94 * t295;
t282 = qJD(4) * t190 - t226 * t299;
t277 = t191 * t286;
t276 = t296 / 0.2e1;
t275 = -t295 / 0.2e1;
t273 = -t288 / 0.2e1;
t272 = t288 / 0.2e1;
t271 = -pkin(6) * t190 - t172;
t270 = qJD(3) * t304;
t269 = qJD(3) * t302;
t267 = t181 - t298;
t266 = t217 + t271;
t265 = qJD(3) * t273;
t264 = qJD(5) * t276;
t178 = t195 * t277;
t48 = rSges(6,1) * t90 + rSges(6,2) * t91 - rSges(6,3) * t280;
t102 = rSges(6,3) * t191 + t190 * t263;
t161 = (-rSges(6,1) * t198 + rSges(6,2) * t197) * t191;
t61 = qJD(3) * t102 + qJD(5) * t161;
t68 = rSges(6,1) * t163 + rSges(6,2) * t162 + rSges(6,3) * t310;
t16 = -t195 * t285 + t48 * t288 - t168 * t61 + t178 + (-t122 * t195 + (t191 * t68 - t217 * t312) * qJD(5)) * qJD(3);
t179 = t196 * t277;
t49 = rSges(6,1) * t92 + rSges(6,2) * t93 - rSges(6,3) * t281;
t69 = -rSges(6,1) * t165 + rSges(6,2) * t164 + rSges(6,3) * t311;
t17 = -t196 * t285 - t49 * t288 + t169 * t61 + t179 + (-t122 * t196 + (-t191 * t69 + t217 * t313) * qJD(5)) * qJD(3);
t261 = -t16 * t196 + t17 * t195;
t260 = -t62 * t168 - t63 * t169;
t22 = t164 * t65 - t165 * t67 + t311 * t63;
t258 = t195 * t22 + t327;
t254 = t197 * t66 + t198 * t64;
t23 = t190 * t62 - t191 * t254;
t253 = t197 * t67 + t198 * t65;
t24 = t190 * t63 - t191 * t253;
t257 = t24 * t195 + t23 * t196;
t252 = qJD(3) * t271;
t35 = t196 * t252 - t288 * t69 + t300 - t364;
t36 = t195 * t252 + t288 * t68 + t267 + t365;
t256 = -t195 * t36 - t196 * t35;
t255 = t195 * t68 - t196 * t69;
t251 = -t197 * t216 - t198 * t215;
t237 = -pkin(6) * t296 - t122 - t61;
t236 = t195 * t265;
t235 = t196 * t265;
t234 = t156 * t295 + t157 * t293 + qJD(1) - t291;
t170 = pkin(4) * t195 + pkin(6) * t310;
t171 = -pkin(4) * t196 + pkin(6) * t311;
t18 = t168 * t69 - t169 * t68 + (t170 * t196 + t171 * t195) * qJD(3) + t234;
t230 = t18 * t255;
t120 = rSges(5,1) * t195 + t176 * t196;
t121 = -rSges(5,1) * t196 + t176 * t195;
t40 = (t120 * t196 + t121 * t195) * qJD(3) + t234;
t229 = t40 * t262;
t219 = (t251 + t96) * t190;
t218 = t158 * t288 + t168 * (Icges(6,5) * t162 - Icges(6,6) * t163) + t169 * (Icges(6,5) * t164 + Icges(6,6) * t165);
t213 = t191 * t218;
t205 = (t164 * t59 - t165 * t60 + t195 * t231 - t215 * t93 - t216 * t92) * t190 + (-t190 * t258 + t325) * qJD(3);
t41 = -t191 * t251 - t329;
t204 = ((qJD(3) * t251 + t58) * t190 + (-qJD(3) * t214 - t197 * t60 - t198 * t59 + (-t197 * t215 + t198 * t216) * qJD(5)) * t191) * t190 + (-t190 * t257 + t41 * t191) * qJD(3);
t203 = (t214 * t196 + t254) * t168 + (t214 * t195 + t253) * t169;
t200 = (qJD(5) * t219 + t203) * t191;
t184 = t191 * t289;
t182 = t191 * t290;
t117 = -t174 * t295 - t298;
t116 = -t174 * t293 + t189;
t87 = t216 * t196;
t86 = t216 * t195;
t85 = t215 * t196;
t84 = t215 * t195;
t77 = rSges(6,1) * t164 + rSges(6,2) * t165;
t76 = rSges(6,1) * t162 - rSges(6,2) * t163;
t57 = t195 * t269 + t267;
t56 = t196 * t269 + t300;
t55 = t360 * qJD(3);
t53 = t196 * t270 + t179;
t52 = t195 * t270 + t178;
t37 = qJD(3) * t348 + t284;
t12 = t168 * t49 - t169 * t48 + (qJD(3) * qJD(5) * t255 - t299 * t336) * t190 + t284;
t11 = t168 * t23 + t169 * t24 + t288 * t41;
t10 = t164 * t45 - t165 * t47 + t195 * t232 + t65 * t93 + t67 * t92;
t9 = t164 * t44 - t165 * t46 + t195 * t233 + t64 * t93 + t66 * t92;
t6 = t168 * t21 + t169 * t22 + t288 * t39;
t5 = t168 * t19 + t169 * t20 + t288 * t38;
t4 = (qJD(3) * t253 + t43) * t190 + (qJD(3) * t63 - t197 * t47 - t198 * t45 + (t197 * t65 - t198 * t67) * qJD(5)) * t191;
t3 = (qJD(3) * t254 + t42) * t190 + (qJD(3) * t62 - t197 * t46 - t198 * t44 + (t197 * t64 - t198 * t66) * qJD(5)) * t191;
t2 = qJD(5) * t205 + t169 * t10 + t168 * t9;
t1 = [-m(4) * t55 + m(5) * t37 + m(6) * t12; m(5) * (t195 * t53 - t196 * t52) + m(6) * t261; -t11 * t287 / 0.2e1 + (((-t197 * t87 - t198 * t85 + t62) * t168 + (-t197 * t86 - t198 * t84 + t63) * t169 + t41 * qJD(5)) * t191 + ((t219 + (-t100 * t197 - t198 * t98 - t214) * t191 - t257) * qJD(5) + t203) * t190) * t273 - t196 * t2 / 0.2e1 + (t195 * t23 - t196 * t24) * t264 + (t19 * t195 - t196 * t20) * t235 + (t195 * t21 - t196 * t22) * t236 + (-t10 * t196 + t195 * t9) * t338 + ((t164 * t85 - t165 * t87) * t168 + (t164 * t84 - t165 * t86) * t169 + (t325 + (-t100 * t165 + t164 * t98 - t327) * t190) * qJD(5) + (((-t22 + t329) * qJD(5) + t260) * t190 + t200) * t195) * t339 + (t195 * t7 - t196 * t8) * t340 + ((t162 * t85 + t163 * t87) * t168 + (t162 * t84 + t163 * t86) * t169 + (t326 + (t100 * t163 + t162 * t98 - t328) * t190) * qJD(5) + (((-t19 + t329) * qJD(5) + t260) * t190 + t200) * t196) * t341 + t195 * t343 + (t353 * t193 + (t358 * t195 + (-t354 - t359) * t196) * t195) * t295 + (-t359 * t193 + (t354 * t195 + (-t353 + t358) * t196) * t195) * t293 + (t356 * qJD(3) * t192 + (-t195 * t357 + t355) * t293) * t275 + ((-t196 * t356 + t355) * t295 + t357 * qJD(3) * t193) * t293 / 0.2e1 + (-t35 * (t102 * t169 + t184) - t36 * (-t102 * t168 + t182) - t18 * (t195 * t365 - t196 * t364 + t282) - (t256 * t175 + (-t18 * t190 * t299 + t191 * t256) * pkin(6)) * qJD(3) - ((-t35 * t69 + t36 * t68) * t191 + t230 * t190) * qJD(5) + t12 * t303 + t18 * (-pkin(6) * t297 * t299 + t331) + (t17 * t266 + t35 * t237 + t12 * (t170 + t68) + t18 * t48) * t196 + (t16 * t266 + t36 * t237 + t12 * (t171 + t69) + t18 * t49) * t195) * m(6) + (-t56 * t184 - t57 * t182 - ((t196 * t229 + t301 * t56) * t196 + (t195 * t229 + t301 * t57) * t195) * qJD(3) + t37 * t303 + (t37 * t120 + t302 * t53 + t304 * t56) * t196 + (t37 * t121 + t302 * t52 + t304 * t57) * t195 + (-t282 + t331 + t348) * t40) * m(5) + (-t55 * t299 + t116 * t293 + t117 * t295 + (-t116 * t196 - t117 * t195 + t360) * qJD(3)) * (rSges(4,1) * t191 - rSges(4,2) * t190) * m(4) + ((-t4 + t5) * t196 + (t3 + t6) * t195) * t272; 0.2e1 * (t12 * t337 + t18 * t342) * m(6) + 0.2e1 * (t37 * t337 + t40 * t342) * m(5) + 0.2e1 * (m(5) * (qJD(3) * t40 + t195 * t52 + t196 * t53) / 0.2e1 + m(6) * (qJD(3) * t18 + t16 * t195 + t17 * t196) / 0.2e1) * t190; -t5 * t280 / 0.2e1 + t310 * t343 + (t190 * t38 + t191 * t259) * t235 + ((t195 * t8 + t196 * t7) * t191 + t206) * t340 + t190 * t6 * t275 + t2 * t311 / 0.2e1 + (t190 * t39 + t191 * t258) * t236 + ((t10 * t195 + t196 * t9) * t191 + t205) * t338 + t11 * t276 + t190 * (qJD(5) * t204 + t168 * t3 + t169 * t4) / 0.2e1 + (t190 * t41 + t191 * t257) * t264 + ((t195 * t4 + t196 * t3) * t191 + t204) * t272 + (t207 * t162 - t163 * t208 + t196 * t213) * t341 + (t164 * t207 + t165 * t208 + t195 * t213) * t339 + (t218 * t190 + (t208 * t197 - t198 * t207) * t191) * t273 + ((t16 * t68 - t17 * t69 - t35 * t49 + t36 * t48 + (t230 - (-t195 * t35 + t196 * t36) * t217) * qJD(3)) * t190 + (t35 * (-qJD(3) * t69 + t195 * t61) + t36 * (qJD(3) * t68 - t196 * t61) - t12 * t255 + t18 * (-t195 * t48 + t196 * t49) - t261 * t217) * t191 - t35 * (t161 * t169 - t288 * t77) - t36 * (-t161 * t168 + t288 * t76) - t18 * (t168 * t77 - t169 * t76)) * m(6);];
tauc = t1(:);
