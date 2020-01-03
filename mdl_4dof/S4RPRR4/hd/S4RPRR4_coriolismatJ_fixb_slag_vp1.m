% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:15
% EndTime: 2019-12-31 16:50:25
% DurationCPUTime: 8.52s
% Computational Cost: add. (28153->442), mult. (34394->660), div. (0->0), fcn. (38001->8), ass. (0->258)
t256 = qJ(1) + pkin(7);
t253 = sin(t256);
t258 = sin(qJ(3));
t355 = t253 * t258;
t254 = cos(t256);
t304 = -sin(qJ(1)) * pkin(1) + t254 * pkin(5);
t261 = cos(qJ(3));
t385 = pkin(3) * t261;
t388 = rSges(5,3) + pkin(6);
t409 = t388 * t258 + pkin(2) + t385;
t257 = sin(qJ(4));
t348 = t257 * t261;
t260 = cos(qJ(4));
t352 = t254 * t260;
t201 = t253 * t348 + t352;
t344 = t260 * t261;
t202 = t253 * t344 - t254 * t257;
t412 = -t202 * rSges(5,1) + t201 * rSges(5,2);
t124 = -t253 * t409 + t304 + t412;
t203 = t253 * t260 - t254 * t348;
t204 = t253 * t257 + t254 * t344;
t299 = t204 * rSges(5,1) + t203 * rSges(5,2);
t386 = cos(qJ(1)) * pkin(1);
t125 = t253 * pkin(5) + t254 * t409 + t299 + t386;
t167 = -rSges(5,1) * t201 - rSges(5,2) * t202;
t168 = rSges(5,1) * t203 - rSges(5,2) * t204;
t429 = m(5) * (-t124 * t167 + t125 * t168);
t353 = t254 * t258;
t367 = Icges(5,4) * t260;
t292 = -Icges(5,2) * t257 + t367;
t207 = -Icges(5,6) * t261 + t292 * t258;
t368 = Icges(5,4) * t257;
t293 = Icges(5,1) * t260 - t368;
t209 = -Icges(5,5) * t261 + t293 * t258;
t290 = Icges(5,5) * t260 - Icges(5,6) * t257;
t205 = -Icges(5,3) * t261 + t290 * t258;
t347 = t258 * t205;
t113 = -t201 * t207 + t202 * t209 + t253 * t347;
t145 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t355;
t197 = Icges(5,4) * t202;
t148 = -Icges(5,2) * t201 + Icges(5,6) * t355 + t197;
t196 = Icges(5,4) * t201;
t152 = -Icges(5,1) * t202 - Icges(5,5) * t355 + t196;
t80 = t145 * t355 - t148 * t201 - t152 * t202;
t147 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t353;
t369 = Icges(5,4) * t204;
t150 = Icges(5,2) * t203 + Icges(5,6) * t353 + t369;
t198 = Icges(5,4) * t203;
t153 = Icges(5,1) * t204 + Icges(5,5) * t353 + t198;
t81 = t147 * t355 - t201 * t150 + t202 * t153;
t297 = t253 * t80 + t254 * t81;
t12 = t261 * t113 - t297 * t258;
t428 = m(5) * qJD(4);
t376 = rSges(5,2) * t257;
t377 = rSges(5,1) * t260;
t298 = -t376 + t377;
t219 = -rSges(5,3) * t261 + t298 * t258;
t154 = rSges(5,3) * t355 - t412;
t360 = t154 * t261;
t126 = t219 * t355 + t360;
t115 = t203 * t207 + t204 * t209 + t254 * t347;
t82 = t145 * t353 + t203 * t148 - t152 * t204;
t83 = t147 * t353 + t203 * t150 + t204 * t153;
t296 = t253 * t82 + t254 * t83;
t426 = -t261 * t115 + t296 * t258;
t42 = t253 * t81 - t254 * t80;
t362 = t145 * t261;
t423 = t148 * t257 + t152 * t260;
t94 = t423 * t258 + t362;
t424 = t253 * t83 - t254 * t82;
t354 = t253 * t261;
t191 = Icges(4,4) * t354 - Icges(4,2) * t355 - Icges(4,6) * t254;
t255 = Icges(4,4) * t261;
t365 = Icges(4,2) * t258;
t192 = Icges(4,6) * t253 + (t255 - t365) * t254;
t370 = Icges(4,4) * t258;
t233 = Icges(4,1) * t261 - t370;
t194 = Icges(4,5) * t253 + t233 * t254;
t175 = t194 * t354;
t229 = Icges(4,5) * t261 - Icges(4,6) * t258;
t190 = Icges(4,3) * t253 + t229 * t254;
t303 = t190 * t254 - t175;
t189 = Icges(4,5) * t354 - Icges(4,6) * t355 - Icges(4,3) * t254;
t240 = Icges(4,4) * t355;
t193 = Icges(4,1) * t354 - Icges(4,5) * t254 - t240;
t351 = t254 * t261;
t332 = -t253 * t189 - t193 * t351;
t416 = -t191 * t353 - t192 * t355 - t303 - t332;
t394 = t253 / 0.2e1;
t392 = -t254 / 0.2e1;
t390 = t254 / 0.2e1;
t158 = (-t388 * t261 + (pkin(3) + t298) * t258) * t253;
t246 = pkin(6) * t351;
t325 = rSges(5,3) * t351 + t353 * t376;
t159 = t246 + (-pkin(3) - t377) * t353 + t325;
t378 = rSges(4,1) * t261;
t305 = pkin(2) + t378;
t324 = rSges(4,2) * t355 + t254 * rSges(4,3);
t169 = -t305 * t253 + t304 + t324;
t242 = rSges(4,2) * t353;
t170 = t386 - t242 + t305 * t254 + (rSges(4,3) + pkin(5)) * t253;
t234 = rSges(4,1) * t258 + rSges(4,2) * t261;
t217 = t234 * t253;
t218 = t234 * t254;
t411 = -m(4) * (t169 * t217 - t170 * t218) - m(5) * (t124 * t158 + t125 * t159);
t225 = (-Icges(5,2) * t260 - t368) * t258;
t226 = (-Icges(5,1) * t257 - t367) * t258;
t408 = -(t209 / 0.2e1 + t225 / 0.2e1) * t257 + (t226 / 0.2e1 - t207 / 0.2e1) * t260;
t407 = t253 ^ 2;
t406 = t254 ^ 2;
t405 = m(5) / 0.2e1;
t404 = -t426 / 0.2e1;
t403 = t42 / 0.2e1;
t402 = t424 / 0.2e1;
t161 = -Icges(5,5) * t201 - Icges(5,6) * t202;
t335 = -Icges(5,2) * t202 - t152 - t196;
t337 = -Icges(5,1) * t201 - t148 - t197;
t72 = -t161 * t261 + (-t335 * t257 + t337 * t260) * t258;
t401 = t72 / 0.2e1;
t184 = t219 * t253;
t220 = rSges(5,3) * t258 + t298 * t261;
t285 = t219 * t261 + t220 * t258;
t103 = -t154 * t258 - t184 * t261 + t285 * t253;
t156 = rSges(5,3) * t353 + t299;
t185 = -rSges(5,1) * t258 * t352 + t325;
t104 = t156 * t258 - t185 * t261 - t285 * t254;
t109 = (t154 * t254 - t156 * t253) * t258;
t359 = t156 * t261;
t128 = t219 * t353 + t359;
t84 = (-t184 * t258 + t360) * t254 + (-t185 * t258 - t359) * t253;
t399 = m(5) * (t103 * t126 - t104 * t128 + t109 * t84);
t398 = m(5) * (t103 * t124 + t104 * t125 + t126 * t158 - t128 * t159);
t244 = pkin(3) * t258 - pkin(6) * t261;
t327 = t219 + t244;
t171 = t327 * t253;
t173 = t327 * t254;
t227 = (-rSges(5,1) * t257 - rSges(5,2) * t260) * t258;
t396 = m(5) * (t167 * t173 - t168 * t171 + (-t124 * t254 - t125 * t253) * t227);
t393 = t253 / 0.4e1;
t391 = -t254 / 0.4e1;
t389 = -t261 / 0.2e1;
t374 = t253 * t12;
t373 = t254 * t426;
t371 = Icges(4,1) * t258;
t361 = t147 * t261;
t357 = t191 * t258;
t356 = t205 * t261;
t350 = t257 * t207;
t208 = Icges(5,6) * t258 + t292 * t261;
t349 = t257 * t208;
t346 = t260 * t209;
t210 = Icges(5,5) * t258 + t293 * t261;
t345 = t260 * t210;
t224 = (-Icges(5,5) * t257 - Icges(5,6) * t260) * t258;
t342 = t261 * t224;
t122 = t167 * t253 + t168 * t254;
t63 = 0.2e1 * (t84 / 0.4e1 - t122 / 0.4e1) * m(5);
t341 = t63 * qJD(2);
t336 = Icges(5,1) * t203 - t150 - t369;
t334 = -Icges(5,2) * t204 + t153 + t198;
t331 = t253 * t190 + t194 * t351;
t329 = -t207 + t226;
t328 = t209 + t225;
t245 = pkin(6) * t258 + t385;
t326 = -t220 - t245;
t323 = qJD(1) * t258;
t322 = qJD(1) * t261;
t321 = qJD(4) * t258;
t320 = qJD(4) * t261;
t59 = t161 * t355 - t335 * t201 + t337 * t202;
t162 = Icges(5,5) * t203 - Icges(5,6) * t204;
t60 = t162 * t355 - t334 * t201 + t336 * t202;
t28 = t253 * t60 - t254 * t59;
t180 = t207 * t253;
t182 = t209 * t253;
t279 = -t205 * t253 + t423;
t268 = t279 * t258 + t362;
t55 = t180 * t201 - t182 * t202 + t253 * t268;
t181 = t207 * t254;
t183 = t209 * t254;
t288 = -t150 * t257 + t153 * t260;
t278 = -t205 * t254 - t288;
t267 = t278 * t258 + t361;
t56 = t181 * t201 - t183 * t202 + t253 * t267;
t206 = Icges(5,3) * t258 + t290 * t261;
t286 = t346 - t350;
t277 = t206 - t286;
t266 = t277 * t258 + t356;
t87 = -t201 * t208 + t202 * t210 + t253 * t266;
t8 = (t297 - t87) * t261 + (t253 * t55 + t254 * t56 + t113) * t258;
t318 = t28 / 0.2e1 - t8 / 0.2e1;
t61 = t161 * t353 + t335 * t203 + t337 * t204;
t62 = t162 * t353 + t334 * t203 + t336 * t204;
t29 = t253 * t62 - t254 * t61;
t57 = -t180 * t203 - t182 * t204 + t254 * t268;
t58 = -t181 * t203 - t183 * t204 + t254 * t267;
t88 = t203 * t208 + t204 * t210 + t254 * t266;
t9 = (t296 - t88) * t261 + (t253 * t57 + t254 * t58 + t115) * t258;
t317 = t29 / 0.2e1 - t9 / 0.2e1;
t313 = t404 + t426 / 0.2e1;
t311 = t355 / 0.4e1;
t302 = t192 * t258 - t189;
t73 = -t162 * t261 + (-t334 * t257 + t336 * t260) * t258;
t96 = -t328 * t201 + t329 * t202 + t224 * t355;
t97 = t328 * t203 + t329 * t204 + t224 * t353;
t301 = t396 / 0.2e1 + (t73 + t97) * t393 + (t72 + t96) * t391;
t95 = t288 * t258 - t361;
t295 = -t94 * t253 + t95 * t254;
t294 = -t255 - t371;
t230 = Icges(4,2) * t261 + t370;
t291 = -Icges(4,5) * t258 - Icges(4,6) * t261;
t121 = -t192 * t353 + t331;
t284 = (t302 * t254 + t121 - t331) * t390 + (-(-t193 * t261 + t357) * t253 - t189 * t254) * t392 - t42 / 0.2e1 + t403 + (t302 * t253 + t303 + t416) * t394;
t283 = t121 * t394 - t331 * t253 / 0.2e1 + t402 - t424 / 0.2e1 + (-t175 + (t190 + t357) * t254 + t332 + t416) * t392;
t213 = -Icges(4,2) * t354 - t240;
t215 = t294 * t253;
t275 = (t191 - t215) * t261 + (t193 + t213) * t258;
t214 = t230 * t254;
t216 = t294 * t254;
t274 = (-t192 + t216) * t261 + (-t194 + t214) * t258;
t273 = t12 * t393 + t426 * t391 - t374 / 0.4e1 + t373 / 0.4e1 + (t311 - t355 / 0.4e1) * t424;
t102 = -t277 * t261 + (t205 + t345 - t349) * t258;
t130 = t286 * t258 - t356;
t67 = -t279 * t261 + (t180 * t257 - t182 * t260 + t145) * t258;
t68 = -t278 * t261 + (t181 * t257 - t183 * t260 + t147) * t258;
t269 = t102 * t389 + t130 * t258 / 0.2e1 + t398 / 0.2e1 + (t67 + t87) * t311 + (t113 - t94) * t354 / 0.4e1 + (t68 + t88) * t353 / 0.4e1 + (t115 + t95) * t351 / 0.4e1;
t265 = -t346 / 0.2e1 + t350 / 0.2e1 - t255 - t371 / 0.2e1 + t365 / 0.2e1 + t206 / 0.2e1;
t264 = -t345 / 0.2e1 + t349 / 0.2e1 - t233 / 0.2e1 + t230 / 0.2e1 - t205 / 0.2e1;
t236 = -rSges(4,2) * t258 + t378;
t212 = t291 * t254;
t211 = t291 * t253;
t174 = t326 * t254;
t172 = t326 * t253;
t160 = -t217 * t253 - t218 * t254;
t138 = -t168 * t261 - t227 * t353;
t137 = t167 * t261 + t227 * t355;
t116 = (t167 * t254 - t168 * t253) * t258;
t112 = -t342 + (-t328 * t257 + t329 * t260) * t258;
t106 = (-pkin(3) * t353 + t185 + t246) * t254 + (-t244 * t253 - t184) * t253;
t99 = (t245 * t253 + t154) * t253 + (t245 * t254 + t156) * t254;
t64 = (t122 + t84) * t405;
t51 = -t342 / 0.2e1 + t429 + t408 * t258;
t49 = t122 * t99 + (t171 * t253 + t173 * t254) * t227;
t39 = -t130 * t261 + t295 * t258;
t32 = -t258 * t264 - t261 * t265 - t411;
t27 = t253 * t58 - t254 * t57;
t26 = t253 * t56 - t254 * t55;
t22 = -t261 * t97 + (t253 * t61 + t254 * t62) * t258;
t21 = -t261 * t96 + (t253 * t59 + t254 * t60) * t258;
t18 = (-t102 + t295) * t261 + (t67 * t253 + t68 * t254 + t130) * t258;
t7 = m(5) * t49 + t28 * t392 + t29 * t394;
t6 = t313 * t355;
t5 = t284 * t253 + t283 * t254;
t4 = t399 + (t373 / 0.2e1 - t374 / 0.2e1 - t18 / 0.2e1) * t261 + (t9 * t390 + t8 * t394 + t39 / 0.2e1) * t258;
t3 = t269 - t396 / 0.2e1 + (t96 / 0.4e1 + t72 / 0.4e1) * t254 + (-t97 / 0.4e1 - t73 / 0.4e1) * t253 + t273;
t2 = t269 + t301;
t1 = (t102 / 0.2e1 + (-t115 / 0.4e1 - t95 / 0.4e1) * t254 + (-t113 / 0.4e1 + t94 / 0.4e1) * t253) * t261 - t398 / 0.2e1 + (-t130 / 0.2e1 + (-t88 / 0.4e1 - t68 / 0.4e1) * t254 + (-t87 / 0.4e1 - t67 / 0.4e1) * t253) * t258 + t273 + t301;
t10 = [t32 * qJD(3) + t51 * qJD(4), 0, t32 * qJD(1) + t2 * qJD(4) + (m(5) * (t124 * t174 + t125 * t172 - t158 * t173 - t159 * t171) + (m(4) * (-t169 * t236 - t217 * t234) - t67 / 0.2e1 - t87 / 0.2e1 + t229 * t390 + (-t193 / 0.2e1 - t213 / 0.2e1) * t261 + (t191 / 0.2e1 - t215 / 0.2e1) * t258 - t283) * t254 + (m(4) * (-t170 * t236 + t218 * t234) + t68 / 0.2e1 + t88 / 0.2e1 + t229 * t394 + (t194 / 0.2e1 - t214 / 0.2e1) * t261 + (-t192 / 0.2e1 + t216 / 0.2e1) * t258 - t284) * t253) * qJD(3), t51 * qJD(1) + t2 * qJD(3) - t112 * t320 + (t124 * t137 + t125 * t138 - t126 * t167 - t128 * t168) * t428 + ((t73 / 0.2e1 + t97 / 0.2e1) * t254 + (t401 + t96 / 0.2e1 - t313) * t253) * t321; 0, 0, 0.2e1 * (m(4) * t160 / 0.2e1 + t106 * t405) * qJD(3) + t64 * qJD(4), t64 * qJD(3) + t116 * t428; t411 * qJD(1) + t5 * qJD(3) + t1 * qJD(4) + t264 * t323 + t265 * t322, -qJD(4) * t63, t5 * qJD(1) + (m(4) * ((t253 * (rSges(4,1) * t354 - t324) + t254 * (rSges(4,1) * t351 + t253 * rSges(4,3) - t242)) * t160 + (t406 + t407) * t236 * t234) + m(5) * (t106 * t99 - t171 * t172 - t173 * t174) + (t407 * t212 + (t275 * t254 + (-t211 + t274) * t253) * t254 + t27) * t394 + (t406 * t211 + (t274 * t253 + (-t212 + t275) * t254) * t253 + t26) * t392) * qJD(3) + t7 * qJD(4), t1 * qJD(1) - t341 + t7 * qJD(3) + (m(5) * (t109 * t122 + t116 * t99 - t137 * t173 - t138 * t171 + (-t126 * t254 + t128 * t253) * t227) + t22 * t394 + t21 * t392 - t399) * qJD(4) + (t18 / 0.2e1 + (t401 + t404) * t254 + (-t73 / 0.2e1 + t12 / 0.2e1) * t253) * t320 + (-t39 / 0.2e1 + t317 * t254 + t318 * t253) * t321; t224 * t322 / 0.2e1 + t3 * qJD(3) + t6 * qJD(4) - qJD(1) * t429 - t408 * t323, qJD(3) * t63, t3 * qJD(1) + t341 + (((t402 + t67 / 0.2e1) * t261 + (t27 / 0.2e1 + t94 / 0.2e1) * t258 + t318) * t254 + ((t403 - t68 / 0.2e1) * t261 + (t26 / 0.2e1 + t95 / 0.2e1) * t258 - t317) * t253 + (-t103 * t173 - t104 * t171 + t106 * t109 + t126 * t174 - t128 * t172 + t84 * t99 - t49) * m(5)) * qJD(3) + t4 * qJD(4), t6 * qJD(1) + t4 * qJD(3) + (m(5) * (t109 * t116 + t126 * t137 - t128 * t138) + t261 ^ 2 * t112 / 0.2e1 + (t22 * t390 + t21 * t394 + (t72 * t253 + t73 * t254) * t389) * t258) * qJD(4);];
Cq = t10;
