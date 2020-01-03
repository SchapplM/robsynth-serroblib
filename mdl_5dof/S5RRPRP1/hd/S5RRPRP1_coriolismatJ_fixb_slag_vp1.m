% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:55
% EndTime: 2020-01-03 11:59:02
% DurationCPUTime: 4.71s
% Computational Cost: add. (19343->277), mult. (12755->350), div. (0->0), fcn. (11248->8), ass. (0->185)
t269 = qJ(1) + qJ(2);
t263 = pkin(8) + t269;
t259 = sin(t263);
t260 = cos(t263);
t271 = sin(qJ(4));
t273 = cos(qJ(4));
t294 = rSges(6,2) * t273 + (rSges(6,1) + pkin(4)) * t271;
t432 = t294 * t260;
t433 = t294 * t259;
t436 = m(6) * (t259 * t433 + t260 * t432);
t121 = -t436 / 0.2e1;
t110 = t436 / 0.2e1;
t435 = Icges(5,2) + Icges(6,2);
t347 = Icges(6,4) * t271;
t348 = Icges(5,4) * t271;
t427 = t435 * t273 + t347 + t348;
t266 = Icges(6,4) * t273;
t267 = Icges(5,4) * t273;
t431 = t266 + t267 + (Icges(5,1) + Icges(6,1)) * t271;
t327 = t260 * t273;
t328 = t260 * t271;
t170 = Icges(6,4) * t327 - Icges(6,2) * t328 + Icges(6,6) * t259;
t172 = Icges(5,4) * t327 - Icges(5,2) * t328 + Icges(5,6) * t259;
t430 = t170 + t172;
t228 = Icges(6,4) * t328;
t174 = Icges(6,1) * t327 + Icges(6,5) * t259 - t228;
t229 = Icges(5,4) * t328;
t176 = Icges(5,1) * t327 + Icges(5,5) * t259 - t229;
t429 = t174 + t176;
t241 = Icges(6,1) * t273 - t347;
t243 = Icges(5,1) * t273 - t348;
t426 = t241 + t243;
t425 = (Icges(5,6) + Icges(6,6)) * t273 + (Icges(5,5) + Icges(6,5)) * t271;
t424 = t435 * t327 + t228 + t229 - t429;
t173 = -Icges(6,5) * t260 + t241 * t259;
t175 = -Icges(5,5) * t260 + t243 * t259;
t423 = -t427 * t259 + t173 + t175;
t422 = t431 * t260 + t430;
t237 = -Icges(6,2) * t271 + t266;
t169 = -Icges(6,6) * t260 + t237 * t259;
t239 = -Icges(5,2) * t271 + t267;
t171 = -Icges(5,6) * t260 + t239 * t259;
t421 = -t431 * t259 - t169 - t171;
t420 = t237 + t239 + t431;
t404 = rSges(6,3) + qJ(5) + pkin(7);
t419 = t404 * t260;
t332 = t259 * t271;
t290 = -rSges(6,2) * t332 - t419;
t359 = pkin(4) * t273;
t262 = pkin(3) + t359;
t352 = rSges(6,1) * t273;
t292 = t262 + t352;
t264 = sin(t269);
t361 = pkin(2) * t264;
t139 = t292 * t259 + t290 + t361;
t362 = sin(qJ(1)) * pkin(1);
t135 = t139 + t362;
t418 = -t135 + t139;
t265 = cos(t269);
t261 = pkin(2) * t265;
t281 = rSges(6,1) * t327 - rSges(6,2) * t328 + t260 * t262;
t140 = t404 * t259 + t261 + t281;
t65 = -t139 * t433 - t140 * t432;
t254 = t260 * pkin(7);
t305 = -rSges(5,2) * t332 - t260 * rSges(5,3);
t353 = rSges(5,1) * t273;
t147 = t361 - t254 + (pkin(3) + t353) * t259 + t305;
t279 = rSges(5,1) * t327 - rSges(5,2) * t328 + rSges(5,3) * t259;
t300 = t260 * pkin(3) + t259 * pkin(7);
t148 = t261 + t279 + t300;
t245 = rSges(5,1) * t271 + rSges(5,2) * t273;
t201 = t245 * t259;
t202 = t245 * t260;
t69 = -t147 * t201 - t148 * t202;
t415 = -m(5) * t69 - m(6) * t65;
t389 = m(5) / 0.2e1;
t388 = m(6) / 0.2e1;
t268 = cos(qJ(1)) * pkin(1);
t368 = m(3) * (-t268 * (rSges(3,1) * t264 + rSges(3,2) * t265) + t362 * (t265 * rSges(3,1) - rSges(3,2) * t264));
t366 = m(4) * (-t268 * (rSges(4,1) * t259 + rSges(4,2) * t260 + t361) + t362 * (t260 * rSges(4,1) - rSges(4,2) * t259 + t261));
t248 = -rSges(5,2) * t271 + t353;
t411 = t248 * t389;
t410 = t425 * t259;
t409 = t425 * t260;
t407 = -t423 * t271 + t421 * t273;
t406 = t424 * t271 - t422 * t273;
t405 = (-t426 + t427) * t273 + t420 * t271;
t298 = qJD(1) + qJD(2);
t256 = t259 ^ 2;
t257 = t260 ^ 2;
t299 = t256 + t257;
t136 = t268 + t140;
t144 = t147 + t362;
t145 = t148 + t268;
t358 = ((-t136 + t140) * t432 + t418 * t433) * t388 + ((-t145 + t148) * t260 + (-t144 + t147) * t259) * t245 * t389;
t64 = -t135 * t433 - t136 * t432;
t68 = -t144 * t201 - t145 * t202;
t401 = (t65 + t64) * t388 + (t69 + t68) * t389;
t278 = t420 * t273 / 0.2e1 + (-t427 / 0.2e1 + t426 / 0.2e1) * t271;
t394 = 4 * qJD(1);
t392 = 2 * qJD(4);
t62 = t144 * t148 - t147 * t145;
t385 = m(5) * t62;
t383 = m(5) * t68;
t128 = t136 * t259;
t129 = t140 * t259;
t378 = m(6) * (t418 * t260 + t128 - t129);
t377 = m(6) * (t128 + t129 + (-t135 - t139) * t260);
t53 = t135 * t140 - t139 * t136;
t376 = m(6) * t53;
t374 = m(6) * t64;
t372 = -t259 / 0.2e1;
t370 = -t260 / 0.2e1;
t369 = t260 / 0.2e1;
t354 = m(6) * qJD(5);
t342 = t169 * t271;
t341 = t170 * t271;
t340 = t171 * t271;
t339 = t172 * t271;
t338 = t173 * t273;
t337 = t175 * t273;
t234 = Icges(6,5) * t273 - Icges(6,6) * t271;
t334 = t234 * t259;
t235 = Icges(5,5) * t273 - Icges(5,6) * t271;
t333 = t235 * t259;
t331 = t259 * t273;
t153 = t173 * t331;
t154 = t174 * t331;
t157 = t169 * t328;
t165 = -Icges(6,3) * t260 + t334;
t284 = t174 * t273 - t341;
t86 = -t165 * t259 - t173 * t327 + t157;
t166 = Icges(6,5) * t327 - Icges(6,6) * t328 + Icges(6,3) * t259;
t87 = t166 * t259 + t284 * t260;
t10 = (t86 + t154 - t157 + (t165 - t341) * t259) * t259 + (-t153 - t87 + (t165 + t284) * t260 + (t338 + t342) * t259) * t260;
t155 = t175 * t331;
t156 = t176 * t331;
t158 = t171 * t328;
t167 = -Icges(5,3) * t260 + t333;
t282 = t176 * t273 - t339;
t88 = -t167 * t259 - t175 * t327 + t158;
t168 = Icges(5,5) * t327 - Icges(5,6) * t328 + Icges(5,3) * t259;
t89 = t168 * t259 + t282 * t260;
t11 = (t88 + t156 - t158 + (t167 - t339) * t259) * t259 + (-t155 - t89 + (t167 + t282) * t260 + (t337 + t340) * t259) * t260;
t82 = -t165 * t260 - t169 * t332 + t153;
t83 = t166 * t260 + t170 * t332 - t154;
t12 = (t157 - t83 + (t166 - t338) * t260) * t260 + (-t153 + t82 + (t166 + t342) * t259) * t259;
t84 = -t167 * t260 - t171 * t332 + t155;
t85 = t168 * t260 + t172 * t332 - t156;
t13 = (t158 - t85 + (t168 - t337) * t260) * t260 + (-t155 + t84 + (t168 + t340) * t259) * t259;
t48 = -t259 * t83 - t260 * t82;
t49 = -t259 * t85 - t260 * t84;
t50 = -t259 * t87 - t260 * t86;
t51 = -t259 * t89 - t260 * t88;
t2 = (-t13 / 0.2e1 - t51 / 0.2e1 - t12 / 0.2e1 - t50 / 0.2e1) * t260 + (t49 / 0.2e1 - t11 / 0.2e1 + t48 / 0.2e1 - t10 / 0.2e1) * t259;
t61 = 0.2e1 * t121;
t297 = t2 * qJD(4) + t61 * qJD(5);
t71 = -t135 * t260 + t128;
t296 = m(6) * t71 * qJD(1);
t75 = -t139 * t260 + t129;
t295 = t75 * m(6) * qJD(2);
t293 = -rSges(6,2) * t271 + t352 + t359;
t277 = t278 + t401;
t276 = -t278 + (t429 * t271 + t430 * t273) * (t369 + t370);
t60 = t121 + t110;
t275 = t60 * qJD(5) + ((t10 + t11) * t259 / 0.2e1 + (t423 * t273 + t421 * t271 + (-t234 - t235) * t260 - t405 * t259) * t370 + (t12 + t13 + t50 + t51) * t369 + (t405 * t260 + t422 * t271 + t424 * t273 - t333 - t334 + t48 + t49) * t372) * qJD(4);
t182 = t293 * t260;
t180 = t293 * t259;
t141 = -t201 * t259 - t202 * t260;
t111 = t294 * t299;
t59 = 0.2e1 * t110;
t57 = t60 * qJD(4);
t55 = t59 * qJD(4);
t35 = t377 / 0.2e1;
t34 = t378 / 0.2e1;
t21 = t278 - t415;
t20 = t278 + t374 + t383;
t17 = t366 + t368 + t376 + t385;
t16 = t35 - t378 / 0.2e1;
t15 = t35 + t34;
t14 = t34 - t377 / 0.2e1;
t5 = t277 + t358;
t4 = t277 - t358;
t3 = t276 + t358 - t401;
t1 = [t17 * qJD(2) + t20 * qJD(4) + t71 * t354, t17 * qJD(1) + t5 * qJD(4) + t15 * qJD(5) + 0.2e1 * (t366 / 0.2e1 + t368 / 0.2e1 + t53 * t388 + t62 * t389) * qJD(2), 0, t20 * qJD(1) + t5 * qJD(2) + ((t144 * t260 - t145 * t259) * t411 + (t135 * t182 - t136 * t180) * t388) * t392 + t275, t15 * qJD(2) + t296 + t57; t4 * qJD(4) + t16 * qJD(5) + (-t368 / 0.4e1 - t366 / 0.4e1 - t385 / 0.4e1 - t376 / 0.4e1) * t394, t21 * qJD(4) + t75 * t354, 0, t4 * qJD(1) + t21 * qJD(2) + ((t139 * t182 - t140 * t180) * t388 + (t147 * t260 - t148 * t259) * t411) * t392 + t275, t16 * qJD(1) + t295 + t57; 0, 0, 0, (-t111 * t388 + t141 * t389) * t392, 0; t276 * qJD(1) + t3 * qJD(2) + (-t374 / 0.4e1 - t383 / 0.4e1) * t394 + t297, t3 * qJD(1) + t297 + (t276 + t415) * qJD(2), 0, (m(5) * (t245 * t248 * t299 + (t260 * t279 + t259 * (rSges(5,1) * t331 + t305)) * t141) + m(6) * (-t111 * ((t281 - t300) * t260 + (t419 + t254 + t290 + (-pkin(3) + t292) * t259) * t259) + t433 * t180 + t432 * t182) + (t409 * t256 + (t407 * t260 + (-t406 - t410) * t259) * t260) * t372 + (-t410 * t257 + (t406 * t259 + (-t407 + t409) * t260) * t259) * t370) * qJD(4) + t298 * t2, t298 * t61; t14 * qJD(2) - t296 + t55, t14 * qJD(1) - t295 + t55, 0, m(6) * (t180 * t260 - t182 * t259) * qJD(4) + t298 * t59, 0;];
Cq = t1;
