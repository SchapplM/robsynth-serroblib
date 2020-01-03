% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:24
% EndTime: 2019-12-31 16:57:33
% DurationCPUTime: 6.19s
% Computational Cost: add. (6112->307), mult. (8329->412), div. (0->0), fcn. (7575->6), ass. (0->195)
t401 = Icges(3,3) + Icges(4,3);
t235 = sin(qJ(1));
t234 = sin(qJ(2));
t335 = pkin(2) * t234;
t231 = qJ(2) + pkin(6);
t216 = cos(t231);
t215 = sin(t231);
t336 = rSges(5,1) + pkin(3);
t277 = t336 * t215;
t322 = rSges(5,3) + qJ(4);
t384 = t322 * t216 - t277;
t258 = t335 - t384;
t398 = t258 * t235;
t400 = t235 * t398;
t237 = cos(qJ(1));
t369 = t322 * t215 + t336 * t216;
t399 = -t237 * rSges(5,2) + t235 * t369;
t171 = Icges(4,5) * t216 - Icges(4,6) * t215;
t236 = cos(qJ(2));
t187 = Icges(3,5) * t236 - Icges(3,6) * t234;
t397 = (t171 + t187) * t237 + t401 * t235;
t210 = Icges(5,5) * t215;
t318 = Icges(5,1) * t216;
t252 = t210 + t318;
t131 = Icges(5,4) * t235 + t237 * t252;
t316 = Icges(4,4) * t215;
t178 = Icges(4,1) * t216 - t316;
t133 = Icges(4,5) * t235 + t178 * t237;
t396 = t131 + t133;
t299 = t235 * t236;
t301 = t234 * t235;
t303 = t216 * t235;
t305 = t215 * t235;
t392 = -Icges(3,5) * t299 - Icges(4,5) * t303 + Icges(3,6) * t301 + Icges(4,6) * t305 + t401 * t237;
t317 = Icges(3,4) * t234;
t191 = Icges(3,1) * t236 - t317;
t145 = Icges(3,5) * t235 + t191 * t237;
t395 = -t133 * t303 - t145 * t299;
t356 = m(5) / 0.2e1;
t181 = rSges(4,1) * t215 + rSges(4,2) * t216;
t240 = t181 + t335;
t374 = t240 * t237;
t375 = t240 * t235;
t379 = t235 * t375 + t237 * t374;
t86 = t258 * t237;
t330 = (-t237 * t86 - t400) * t356 - m(4) * t379 / 0.2e1;
t357 = m(4) / 0.2e1;
t302 = t216 * t237;
t282 = t322 * t302;
t81 = (-t277 - t335) * t237 + t282;
t332 = (-t237 * t81 + t400) * t356 + t379 * t357;
t11 = t332 - t330;
t394 = t11 * qJD(1);
t391 = t237 * t397 + t395;
t173 = Icges(4,2) * t216 + t316;
t312 = Icges(5,3) * t216;
t248 = t312 - t210;
t390 = (-t173 - t248) * t237 + t396;
t200 = Icges(4,4) * t305;
t132 = Icges(4,1) * t303 - Icges(4,5) * t237 - t200;
t207 = Icges(3,4) * t301;
t144 = Icges(3,1) * t299 - Icges(3,5) * t237 - t207;
t298 = t236 * t237;
t389 = -t132 * t302 - t144 * t298 + t235 * t392;
t128 = Icges(4,4) * t303 - Icges(4,2) * t305 - Icges(4,6) * t237;
t142 = Icges(3,4) * t299 - Icges(3,2) * t301 - Icges(3,6) * t237;
t388 = t128 * t215 + t142 * t234;
t199 = Icges(5,5) * t302;
t304 = t215 * t237;
t123 = Icges(5,6) * t235 + Icges(5,3) * t304 + t199;
t172 = Icges(5,4) * t216 + Icges(5,6) * t215;
t127 = Icges(5,2) * t235 + t172 * t237;
t387 = t123 * t304 + t145 * t298 + t396 * t302 + (t127 + t397) * t235;
t386 = -Icges(3,5) * t234 - Icges(3,6) * t236 + (-Icges(4,6) + Icges(5,6)) * t216 + (-Icges(5,4) - Icges(4,5)) * t215;
t232 = t235 ^ 2;
t233 = t237 ^ 2;
t279 = t232 + t233;
t211 = Icges(4,4) * t216;
t313 = Icges(4,2) * t215;
t129 = Icges(4,6) * t235 + (t211 - t313) * t237;
t226 = Icges(3,4) * t236;
t189 = -Icges(3,2) * t234 + t226;
t143 = Icges(3,6) * t235 + t189 * t237;
t385 = t129 * t215 + t143 * t234 + t392;
t300 = t234 * t237;
t383 = -t129 * t304 - t143 * t300 + t387;
t310 = (-Icges(5,2) * t237 + t235 * t172) * t237;
t382 = t310 + t387;
t381 = -t128 * t304 - t129 * t305 - t142 * t300 - t143 * t301 - t389 - t391;
t315 = Icges(5,5) * t216;
t170 = Icges(5,3) * t215 + t315;
t175 = Icges(5,1) * t215 - t315;
t188 = Icges(3,2) * t236 + t317;
t319 = Icges(4,1) * t215;
t380 = -(t191 / 0.2e1 - t188 / 0.2e1) * t234 - (t211 + t319 / 0.2e1 - t313 / 0.2e1 + t175 / 0.2e1 - t170 / 0.2e1) * t216 - (t178 / 0.2e1 - t173 / 0.2e1 + t210 + t318 / 0.2e1 - t312 / 0.2e1) * t215;
t378 = -t235 / 0.2e1;
t339 = t235 / 0.2e1;
t377 = -t237 / 0.2e1;
t122 = -Icges(5,6) * t237 + t170 * t235;
t130 = -Icges(5,4) * t237 + t252 * t235;
t371 = (t122 * t215 + t130 * t216) * t235;
t190 = Icges(3,1) * t234 + t226;
t368 = t390 * t235;
t253 = -t211 - t319;
t269 = (t253 * t237 - t129) * t235;
t271 = (-Icges(5,1) * t304 + t123 + t199) * t235;
t367 = t269 + t271;
t366 = qJD(2) * t237;
t365 = t386 * t235;
t364 = t386 * t237;
t164 = t188 * t237;
t166 = t190 * t237;
t361 = (-t145 + t164) * t301 + (-t143 - t166) * t299;
t163 = -Icges(3,2) * t299 - t207;
t165 = t190 * t235;
t360 = (t144 + t163) * t234 + (t142 + t165) * t236;
t359 = 0.4e1 * qJD(1);
t358 = 0.2e1 * qJD(2);
t230 = t237 * pkin(5);
t325 = rSges(3,1) * t236;
t276 = pkin(1) + t325;
t280 = rSges(3,2) * t301 + t237 * rSges(3,3);
t106 = -t276 * t235 + t230 + t280;
t209 = rSges(3,2) * t300;
t107 = -t209 + t276 * t237 + (rSges(3,3) + pkin(5)) * t235;
t192 = rSges(3,1) * t234 + rSges(3,2) * t236;
t167 = t192 * t235;
t168 = t192 * t237;
t355 = m(3) * (t106 * t167 - t107 * t168);
t243 = rSges(4,1) * t303 - rSges(4,2) * t305 - t237 * rSges(4,3);
t334 = pkin(2) * t236;
t214 = pkin(1) + t334;
t333 = -qJ(3) - pkin(5);
t281 = -t235 * t214 - t237 * t333;
t89 = -t243 + t281;
t212 = t235 * t333;
t274 = -rSges(4,2) * t304 + t235 * rSges(4,3);
t324 = rSges(4,1) * t216;
t90 = -t212 + (t214 + t324) * t237 + t274;
t353 = m(4) * (-t374 * t90 + t375 * t89);
t352 = m(4) * (t90 * t235 + t237 * t89);
t69 = t281 - t399;
t323 = t235 * rSges(5,2);
t70 = t323 - t212 + (t214 + t369) * t237;
t329 = t69 * t302 + t70 * t303;
t346 = m(5) * ((t235 * t81 + t237 * t398) * t215 + t329);
t345 = m(5) * (-t304 * t398 + t86 * t305 + t329);
t344 = m(5) * (t398 * t69 + t70 * t81);
t343 = m(5) * (t70 * t235 + t237 * t69);
t328 = -t86 * t302 - t303 * t398;
t327 = m(5) * qJD(2);
t326 = m(5) * qJD(4);
t306 = t215 * t216;
t160 = t279 * t215;
t82 = m(5) * t160;
t296 = t82 * qJD(1);
t294 = -t235 * (pkin(1) * t235 - t230 + t281) + t237 * (-t235 * pkin(5) - t212 + (-pkin(1) + t214) * t237);
t291 = -t175 * t235 + t122;
t290 = -t253 * t235 + t128;
t289 = -t248 * t235 + t130;
t287 = -Icges(4,2) * t303 + t132 - t200;
t283 = t279 * t306;
t34 = t70 * t304 - t69 * t305;
t278 = m(5) * t34 * qJD(1);
t275 = rSges(4,2) * t215 - t324 - t334;
t272 = t291 * t237;
t270 = t290 * t237;
t268 = t289 * t237;
t266 = t287 * t237;
t259 = t279 * t335;
t257 = -t334 - t369;
t256 = -t123 * t305 + t127 * t237 - t131 * t303;
t255 = t172 / 0.2e1 + t171 / 0.2e1 + t187 / 0.2e1;
t43 = -t310 + t371;
t239 = t256 * t378 + (t43 + t392 * t237 + (t132 * t216 + t144 * t236 - t388) * t235) * t377 + (t385 * t237 - t382 + t383) * t237 / 0.2e1 + (t122 * t304 + t130 * t302 + t385 * t235 + t256 + t381 + t391) * t339;
t238 = (t43 - t371 + t382) * t378 + t383 * t339 + ((t388 + t397) * t237 + t381 + t389 + t395) * t377;
t194 = -rSges(3,2) * t234 + t325;
t120 = t275 * t237;
t118 = t275 * t235;
t91 = t283 - t306;
t87 = t257 * t237;
t85 = t257 * t235;
t38 = -t259 + (-t237 * t277 + t282) * t237 + t384 * t232;
t32 = (t237 * t369 + t323) * t237 + t294 + t399 * t235;
t24 = t343 + t352;
t22 = t345 / 0.2e1;
t20 = t346 / 0.2e1;
t17 = t160 * t32 + t328;
t13 = t330 + t332;
t5 = (t190 / 0.2e1 + t189 / 0.2e1) * t236 + t355 + t353 + t344 - t380;
t4 = t22 - t346 / 0.2e1;
t3 = t22 + t20;
t2 = t20 - t345 / 0.2e1;
t1 = t239 * t235 + t237 * t238;
t6 = [t5 * qJD(2) + t24 * qJD(3) + t34 * t326, t5 * qJD(1) + t13 * qJD(3) + t3 * qJD(4) + ((t69 * t87 + t70 * t85 + (-t81 - t86) * t398) * t356 + (t118 * t90 + t120 * t89) * t357) * t358 + (m(3) * (-t106 * t194 - t167 * t192) + t255 * t237 + (-t144 / 0.2e1 - t163 / 0.2e1) * t236 + (t142 / 0.2e1 + t165 / 0.2e1) * t234 - t238) * t366 + ((m(3) * (-t107 * t194 + t168 * t192) + t255 * t235 + (t145 / 0.2e1 - t164 / 0.2e1) * t236 + (-t143 / 0.2e1 - t166 / 0.2e1) * t234 - t239) * t235 + (t271 / 0.2e1 + t269 / 0.2e1 - t272 / 0.2e1 + t270 / 0.2e1) * t215 + (-t268 / 0.2e1 - t266 / 0.2e1 + t390 * t339) * t216) * qJD(2), qJD(1) * t24 + qJD(2) * t13, t3 * qJD(2) + t278; t1 * qJD(2) - t11 * qJD(3) + t4 * qJD(4) + (-t355 / 0.4e1 - t353 / 0.4e1 - t344 / 0.4e1) * t359 + (-(t190 + t189) * t236 / 0.2e1 + t380) * qJD(1), t1 * qJD(1) + t17 * t326 - (((t270 - t272 + t367) * t216 - t368 * t215 + ((t287 + t289) * t215 + t360 - t364) * t237 + t361) * t235 + t365 * t233) * t366 / 0.2e1 + (m(3) * ((t235 * (rSges(3,1) * t299 - t280) + t237 * (rSges(3,1) * t298 + t235 * rSges(3,3) - t209)) * (-t235 * t167 - t168 * t237) + t279 * t194 * t192) + m(5) * (t32 * t38 - t398 * t85 - t86 * t87) + m(4) * (-t374 * t120 - t375 * t118 + (t235 * t243 + t237 * (rSges(4,1) * t302 + t274) + t294) * (-t181 * t279 - t259)) + ((t360 * t237 + (t266 + t268 - t368) * t215 - t365 * t235 + ((t290 - t291) * t237 + t367) * t216 + t361) * t237 + t364 * t232) * t339) * qJD(2), -t394, t4 * qJD(1) + t17 * t327 + (-t160 * t216 + t283 - t91) * t326; t11 * qJD(2) - t82 * qJD(4) + (-t352 / 0.4e1 - t343 / 0.4e1) * t359, t394 + ((t235 * t87 - t237 * t85) * t356 + (-t118 * t237 + t235 * t120) * t357) * t358, 0, -t296; t2 * qJD(2) + t82 * qJD(3) - t278, t2 * qJD(1) + (-t216 * t38 + (t235 * t85 + t237 * t87 + t32) * t215 - t17 + t328) * t327 + t91 * t326, t296, t91 * t327;];
Cq = t6;
