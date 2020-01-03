% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR10_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:26:03
% DurationCPUTime: 5.90s
% Computational Cost: add. (9159->383), mult. (10935->472), div. (0->0), fcn. (10970->8), ass. (0->214)
t195 = cos(qJ(3));
t193 = sin(qJ(3));
t343 = cos(qJ(1));
t284 = t343 * t193;
t342 = sin(qJ(1));
t374 = t342 * t195 - t284;
t194 = cos(qJ(5));
t192 = sin(qJ(5));
t291 = qJ(3) + pkin(8);
t268 = sin(t291);
t269 = cos(t291);
t128 = -t342 * t268 - t343 * t269;
t129 = t343 * t268 - t342 * t269;
t323 = Icges(6,4) * t194;
t244 = -Icges(6,2) * t192 + t323;
t62 = -Icges(6,6) * t129 + t128 * t244;
t327 = t192 * t62;
t324 = Icges(6,4) * t192;
t246 = Icges(6,1) * t194 - t324;
t66 = -Icges(6,5) * t129 + t128 * t246;
t247 = t194 * t66 - t327;
t63 = Icges(6,6) * t128 + t129 * t244;
t328 = t192 * t63;
t67 = Icges(6,5) * t128 + t129 * t246;
t373 = -t194 * t67 + t328;
t309 = t129 * t194;
t242 = Icges(6,5) * t194 - Icges(6,6) * t192;
t59 = Icges(6,3) * t128 + t129 * t242;
t338 = -t128 * t59 - t309 * t67;
t58 = -Icges(6,3) * t129 + t128 * t242;
t372 = -(t58 - t328) * t129 + t338;
t295 = qJD(5) * t129;
t296 = qJD(5) * t128;
t191 = qJD(1) - qJD(3);
t344 = -t191 / 0.2e1;
t36 = -t192 * t67 - t194 * t63;
t37 = -t192 * t66 - t194 * t62;
t371 = ((-t37 * t128 + t36 * t129) * qJD(5) - t36 * t295 + t37 * t296) * t344;
t369 = t129 * t58;
t364 = t59 * t129;
t362 = t59 + t327;
t101 = t191 * t128;
t102 = t191 * t129;
t330 = t102 * rSges(5,1) - t101 * rSges(5,2);
t355 = t129 * rSges(5,1) - t128 * rSges(5,2);
t361 = t191 * t355 - t330;
t255 = rSges(6,1) * t194 - rSges(6,2) * t192;
t138 = t255 * qJD(5);
t180 = pkin(3) * t195 + pkin(2);
t275 = qJD(1) * t342;
t276 = qJD(1) * t343;
t357 = t374 * qJD(3) + t193 * t276;
t261 = t357 * pkin(3) - t180 * t275;
t264 = pkin(2) * t275;
t210 = t264 + t261;
t196 = qJD(1) ^ 2;
t286 = t342 * pkin(2);
t182 = qJD(2) * t342;
t299 = qJ(2) * t276 + t182;
t224 = -pkin(1) * t275 + t299;
t305 = (t182 + t224) * qJD(1);
t228 = -t196 * t286 + t305;
t218 = t191 * t210 + t228;
t254 = rSges(6,1) * t192 + rSges(6,2) * t194;
t278 = -t102 * pkin(4) - t101 * pkin(7);
t293 = qJD(5) * t192;
t231 = t102 * t194 + t128 * t293;
t230 = t231 * rSges(6,1) + t101 * rSges(6,3);
t292 = qJD(5) * t194;
t232 = -t102 * t192 + t128 * t292;
t30 = rSges(6,2) * t232 + t230;
t11 = (-t278 + t30) * t191 + (t101 * t254 + t129 * t138) * qJD(5) + t218;
t189 = t343 * pkin(2);
t183 = qJD(2) * t343;
t297 = t343 * pkin(1) + t342 * qJ(2);
t351 = qJD(1) * t297 - t183;
t267 = (t183 - t351) * qJD(1);
t211 = -t189 * t196 + t267;
t277 = -t101 * pkin(4) + t102 * pkin(7);
t233 = -t101 * t194 + t129 * t293;
t229 = t233 * rSges(6,1) + t102 * rSges(6,3);
t234 = t101 * t192 + t129 * t292;
t29 = rSges(6,2) * t234 + t229;
t282 = t342 * t193;
t173 = pkin(3) * t282;
t283 = t343 * t195;
t358 = pkin(3) * t283 + t173;
t262 = -qJD(1) * t173 + t358 * qJD(3) - t180 * t276;
t265 = pkin(2) * t276;
t90 = -t265 - t262;
t12 = (-t102 * t254 + t128 * t138) * qJD(5) + (-t277 - t29 - t90) * t191 + t211;
t239 = pkin(3) * t284 - t342 * t180;
t116 = t286 + t239;
t185 = t343 * qJ(2);
t287 = t342 * pkin(1);
t154 = t287 - t185;
t214 = t182 + (-t286 - t154) * qJD(1);
t294 = qJD(5) * t254;
t104 = -t129 * pkin(4) - t128 * pkin(7);
t123 = t128 * rSges(6,3);
t306 = rSges(6,1) * t309 + t123;
t310 = t129 * t192;
t72 = rSges(6,2) * t310 - t306;
t325 = t104 + t72;
t21 = t128 * t294 + (t116 - t325) * t191 + t214;
t300 = t343 * t180 + t173;
t117 = -t189 + t300;
t354 = t189 + t297;
t212 = qJD(1) * t354 - t183;
t270 = -t128 * pkin(4) + pkin(7) * t129;
t120 = t129 * rSges(6,3);
t312 = t128 * t194;
t304 = -rSges(6,1) * t312 + t120;
t313 = t128 * t192;
t73 = rSges(6,2) * t313 + t304;
t22 = t129 * t294 + (t117 + t270 + t73) * t191 + t212;
t360 = (t192 * (-t101 * t21 - t102 * t22 + t11 * t128 - t12 * t129) + (t128 * t22 - t129 * t21) * t292) * rSges(6,2);
t216 = t283 + t282;
t257 = -rSges(4,1) * t374 + rSges(4,2) * t216;
t359 = t191 * t257;
t241 = Icges(6,5) * t192 + Icges(6,6) * t194;
t81 = t241 * t128;
t80 = t241 * t129;
t356 = 0.2e1 * qJD(5);
t243 = Icges(6,2) * t194 + t324;
t245 = Icges(6,1) * t192 + t323;
t240 = -t192 * t243 + t194 * t245;
t46 = t129 * t240 + t81;
t42 = t46 * t191;
t47 = t128 * t240 - t80;
t43 = t47 * t191;
t303 = -t128 * rSges(5,1) - t129 * rSges(5,2);
t298 = -qJD(1) * t154 + t182;
t225 = -t264 + t298;
t353 = -t191 * t116 + t224 - t225;
t221 = t343 * rSges(3,1) + t342 * rSges(3,3);
t352 = t297 + t221;
t332 = t243 * t129 - t67;
t334 = -t245 * t129 - t63;
t350 = t192 * t332 + t194 * t334;
t337 = t128 * t58 + t309 * t66;
t336 = t312 * t67 - t364;
t335 = t312 * t66 - t369;
t333 = -t245 * t128 - t62;
t331 = t243 * t128 - t66;
t329 = t101 * rSges(5,1) + t102 * rSges(5,2);
t326 = t22 * t254;
t308 = t242 * t191;
t307 = t355 + t116;
t302 = -t243 + t246;
t301 = -t244 - t245;
t290 = -t90 + t329;
t289 = -t343 / 0.2e1;
t288 = t342 / 0.2e1;
t285 = t342 * rSges(3,1);
t272 = t296 / 0.2e1;
t271 = -t295 / 0.2e1;
t263 = -t300 - t303;
t107 = t191 * t216;
t108 = -t195 * t275 + t357;
t57 = t108 * rSges(4,1) + t107 * rSges(4,2);
t56 = t107 * rSges(4,1) - t108 * rSges(4,2);
t256 = -rSges(4,1) * t216 - rSges(4,2) * t374;
t252 = -t128 * t21 - t129 * t22;
t251 = t128 * t73 + t129 * t72;
t15 = -t310 * t63 - t338;
t16 = -t310 * t62 + t337;
t227 = (-t128 * t15 + t129 * t16) * qJD(5);
t17 = -t313 * t63 + t336;
t18 = -t313 * t62 + t335;
t226 = (-t128 * t17 + t129 * t18) * qJD(5);
t223 = -t287 - t286;
t222 = -t285 - t287;
t217 = -t270 - t300 - t304;
t213 = t192 * t331 + t194 * t333;
t207 = t101 * t72 - t102 * t73 + t128 * t30 + t129 * t29;
t206 = (t192 * t301 + t194 * t302) * t191;
t205 = t104 - t239 - t306;
t134 = t244 * qJD(5);
t135 = t246 * qJD(5);
t200 = -t134 * t192 + t135 * t194 + (-t192 * t245 - t194 * t243) * qJD(5);
t199 = -t230 - t261 + t278;
t198 = t229 - t262 + t277;
t10 = qJD(5) * t247 - t192 * (Icges(6,1) * t231 + Icges(6,4) * t232 + Icges(6,5) * t101) - t194 * (Icges(6,4) * t231 + Icges(6,2) * t232 + Icges(6,6) * t101);
t133 = t242 * qJD(5);
t13 = t101 * t240 - t102 * t241 + t128 * t133 + t129 * t200;
t14 = -t101 * t241 - t102 * t240 + t128 * t200 - t129 * t133;
t5 = t42 + t227;
t6 = t43 + t226;
t9 = -t373 * qJD(5) - t192 * (Icges(6,1) * t233 + Icges(6,4) * t234 + Icges(6,5) * t102) - t194 * (Icges(6,4) * t233 + Icges(6,2) * t234 + Icges(6,6) * t102);
t197 = t191 * (-t135 * t192 + t243 * t293 + (-qJD(5) * t245 - t134) * t194) + t5 * t295 / 0.2e1 + (t14 + t10) * t271 + (t13 + t6 + t9) * t272 - ((t47 - t37) * t101 + (t46 - t36) * t102) * qJD(5) / 0.2e1;
t187 = t343 * rSges(3,3);
t181 = rSges(3,3) * t276;
t127 = t191 * t374 * pkin(3);
t100 = -t196 * t221 + t267;
t99 = qJD(1) * (-rSges(3,1) * t275 + t181) + t305;
t87 = t254 * t128;
t86 = t254 * t129;
t78 = -t191 * t256 + t212;
t77 = t214 + t359;
t71 = t129 * t255 + t123;
t70 = t128 * t255 - t120;
t45 = (t117 + t303) * t191 + t212;
t44 = t191 * t307 + t214;
t39 = -t191 * t56 + t211;
t38 = t191 * t57 + t228;
t33 = t191 * t290 + t211;
t32 = t191 * t330 + t218;
t31 = qJD(5) * t251 - qJD(4);
t24 = Icges(6,5) * t231 + Icges(6,6) * t232 + Icges(6,3) * t101;
t23 = Icges(6,5) * t233 + Icges(6,6) * t234 + Icges(6,3) * t102;
t8 = (-t128 * t63 - t129 * t62) * t192 + t336 + t337;
t7 = t207 * qJD(5);
t1 = [(-t42 + ((t17 + (-t247 + t59) * t129) * t129 + (t362 * t128 + t18 - t335 - t369) * t128) * qJD(5)) * t271 + (t43 + ((t373 * t129 + t15 + t335) * t129 + (-t362 * t129 + t16 - t8) * t128) * qJD(5)) * t272 + t371 - t197 + (t12 * (-t205 - t154) + t21 * (-pkin(1) * t276 - qJ(2) * t275 + t183 - t198) + t11 * (-t217 + t297) + t360 - (t326 + t31 * (t71 + t72)) * t296 + (-t199 + t21 - (t71 - t104) * t191 + t353) * t22) * m(6) + (t33 * (t239 - t154 + t355) + t44 * (t262 + t329 - t351) + t32 * (-t263 + t297) + (t261 + t353 - t361 + t44) * t45) * m(5) + (t39 * (t185 + t223 + t257) + t38 * (-t256 + t354) + (-t212 - t56) * t77 + (t223 * qJD(1) - t225 + t299 - t359 + t57 + t77) * t78) * m(4) + (t100 * (t185 + t187 + t222) + t99 * t352 + (t181 - t298 + t299 + (t285 - t187 + t222) * qJD(1)) * (qJD(1) * t352 - t183)) * m(3); 0.2e1 * (t11 * t289 + t12 * t288) * m(6) + 0.2e1 * (t288 * t33 + t289 * t32) * m(5) + 0.2e1 * (t288 * t39 + t289 * t38) * m(4) + 0.2e1 * (t100 * t288 + t289 * t99) * m(3); (t42 + ((-t17 + t8) * t129 + (t247 * t128 - t18 + t372) * t128) * qJD(5)) * t271 + (-t43 + ((-t15 - t372) * t129 + (-t16 + (-t373 + t58) * t128 - t364) * t128) * qJD(5)) * t272 + t371 + t197 + (t12 * (-t286 + t205) + t21 * (-t265 + t198) + t11 * (t189 + t217) - t360 - (-t128 * t326 + (t21 * t254 + t31 * (t70 + t73)) * t129) * qJD(5) - t21 * (t358 - t70 + t270) * t191 + (-t191 * t325 - t127 + t199 - t264) * t22) * m(6) + (-t33 * t307 + t32 * (t189 + t263) + (-t127 - t210 + t361) * t45 + (-t290 - (t358 + t303) * t191) * t44) * m(5) + (-t257 * t39 + t38 * t256 + t56 * t77 - t57 * t78 - (-t256 * t77 - t257 * t78) * t191) * m(4); -m(6) * t7; t191 * (t10 * t129 - t101 * t37 - t102 * t36 - t128 * t9) / 0.2e1 + ((t81 * t295 - t308) * t129 + (t206 + (-t350 * t128 + (-t80 + t213) * t129) * qJD(5)) * t128) * t271 + ((t80 * t296 + t308) * t128 + (t206 + (t213 * t129 + (-t350 - t81) * t128) * qJD(5)) * t129) * t272 + ((t192 * t302 - t194 * t301) * t191 + ((t128 * t332 - t129 * t331) * t194 + (-t128 * t334 + t129 * t333) * t192) * qJD(5)) * t344 - (t13 * t191 + (-(-t101 * t373 - t102 * t59 - t128 * t23) * t128 + t101 * t16 + t102 * t15 + t129 * (t101 * t247 - t102 * t58 - t128 * t24)) * t356) * t128 / 0.2e1 + (t14 * t191 + (t101 * t18 + t102 * t17 - t128 * (-t101 * t59 + t102 * t373 + t129 * t23) + t129 * (-t101 * t58 - t102 * t247 + t129 * t24)) * t356) * t129 / 0.2e1 + (t6 + t226) * t101 / 0.2e1 + (t5 + t227) * t102 / 0.2e1 + (t7 * t251 + t31 * t207 - t252 * t138 - (-t101 * t22 + t102 * t21 - t11 * t129 - t12 * t128) * t254 - (-t21 * t86 + t22 * t87) * t191 - (t31 * (t128 * t87 + t129 * t86) - t252 * t255) * qJD(5)) * m(6);];
tauc = t1(:);
