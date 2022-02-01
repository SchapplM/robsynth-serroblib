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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:19:34
% EndTime: 2022-01-20 10:19:44
% DurationCPUTime: 4.54s
% Computational Cost: add. (19343->270), mult. (12755->349), div. (0->0), fcn. (11248->8), ass. (0->177)
t280 = qJ(1) + qJ(2);
t275 = pkin(8) + t280;
t272 = sin(t275);
t273 = cos(t275);
t281 = sin(qJ(4));
t283 = cos(qJ(4));
t438 = rSges(6,2) * t283 + (rSges(6,1) + pkin(4)) * t281;
t445 = t438 * t273;
t446 = t438 * t272;
t449 = m(6) * (t272 * t446 + t273 * t445);
t121 = -t449 / 0.2e1;
t110 = t449 / 0.2e1;
t448 = Icges(5,2) + Icges(6,2);
t362 = Icges(6,4) * t281;
t363 = Icges(5,4) * t281;
t441 = t448 * t283 + t362 + t363;
t278 = Icges(6,4) * t283;
t279 = Icges(5,4) * t283;
t444 = t278 + t279 + (Icges(5,1) + Icges(6,1)) * t281;
t345 = t272 * t283;
t346 = t272 * t281;
t176 = Icges(6,4) * t345 - Icges(6,2) * t346 - Icges(6,6) * t273;
t178 = Icges(5,4) * t345 - Icges(5,2) * t346 - Icges(5,6) * t273;
t443 = t176 + t178;
t235 = Icges(6,4) * t346;
t180 = Icges(6,1) * t345 - Icges(6,5) * t273 - t235;
t236 = Icges(5,4) * t346;
t182 = Icges(5,1) * t345 - Icges(5,5) * t273 - t236;
t442 = t180 + t182;
t248 = Icges(6,1) * t283 - t362;
t250 = Icges(5,1) * t283 - t363;
t440 = t248 + t250;
t244 = -Icges(6,2) * t281 + t278;
t246 = -Icges(5,2) * t281 + t279;
t439 = -t246 - t244 - t444;
t435 = (-Icges(5,6) - Icges(6,6)) * t283 + (-Icges(5,5) - Icges(6,5)) * t281;
t181 = Icges(6,5) * t272 + t248 * t273;
t183 = Icges(5,5) * t272 + t250 * t273;
t434 = -t273 * t441 + t181 + t183;
t433 = -t448 * t345 - t235 - t236 + t442;
t177 = Icges(6,6) * t272 + t244 * t273;
t179 = Icges(5,6) * t272 + t246 * t273;
t432 = -t273 * t444 - t177 - t179;
t431 = t444 * t272 + t443;
t376 = pkin(4) * t283;
t274 = pkin(3) + t376;
t375 = -qJ(5) - pkin(7);
t295 = rSges(6,1) * t345 - rSges(6,2) * t346 + t272 * t274 + (-rSges(6,3) + t375) * t273;
t276 = sin(t280);
t379 = pkin(2) * t276;
t142 = -t295 - t379;
t367 = rSges(6,1) * t283;
t309 = t274 + t367;
t344 = t273 * t281;
t316 = -rSges(6,2) * t344 - t272 * t375;
t277 = cos(t280);
t378 = pkin(2) * t277;
t143 = t272 * rSges(6,3) + t309 * t273 + t316 + t378;
t65 = t142 * t446 - t143 * t445;
t268 = t273 * pkin(7);
t368 = rSges(5,1) * t283;
t312 = pkin(3) + t368;
t317 = rSges(5,2) * t346 + t273 * rSges(5,3);
t150 = -t312 * t272 + t268 + t317 - t379;
t240 = rSges(5,2) * t344;
t151 = t378 - t240 + t312 * t273 + (rSges(5,3) + pkin(7)) * t272;
t254 = rSges(5,1) * t281 + rSges(5,2) * t283;
t208 = t254 * t272;
t209 = t254 * t273;
t69 = t150 * t208 - t151 * t209;
t430 = -m(5) * t69 - m(6) * t65;
t411 = m(5) / 0.2e1;
t410 = m(6) / 0.2e1;
t380 = cos(qJ(1)) * pkin(1);
t381 = sin(qJ(1)) * pkin(1);
t388 = m(3) * (t380 * (-rSges(3,1) * t276 - rSges(3,2) * t277) + (t277 * rSges(3,1) - t276 * rSges(3,2)) * t381);
t386 = m(4) * (t380 * (-rSges(4,1) * t272 - rSges(4,2) * t273 - t379) + (t273 * rSges(4,1) - t272 * rSges(4,2) + t378) * t381);
t75 = t142 * t273 + t143 * t272;
t426 = t75 * m(6) * qJD(2);
t425 = (t440 - t441) * t283 + t439 * t281;
t138 = t142 - t381;
t139 = t143 + t380;
t71 = t138 * t273 + t139 * t272;
t423 = t71 * m(6) * qJD(1);
t422 = t435 * t272;
t421 = t435 * t273;
t269 = t272 ^ 2;
t270 = t273 ^ 2;
t315 = t269 + t270;
t420 = -t281 * t434 + t283 * t432;
t419 = t281 * t433 + t283 * t431;
t314 = qJD(1) + qJD(2);
t147 = t150 - t381;
t148 = t151 + t380;
t374 = ((-t139 + t143) * t445 + (t138 - t142) * t446) * t410 + ((-t148 + t151) * t273 + (t147 - t150) * t272) * t254 * t411;
t64 = t138 * t446 - t139 * t445;
t68 = t147 * t208 - t148 * t209;
t418 = (t65 + t64) * t410 + (t69 + t68) * t411;
t294 = -t439 * t283 / 0.2e1 + (-t441 / 0.2e1 + t440 / 0.2e1) * t281;
t416 = 0.4e1 * qJD(1);
t414 = 2 * qJD(4);
t62 = -t151 * t147 + t148 * t150;
t407 = m(5) * t62;
t405 = m(5) * t68;
t400 = m(6) * (t71 - t75);
t399 = m(6) * (t75 + t71);
t53 = -t143 * t138 + t139 * t142;
t398 = m(6) * t53;
t396 = m(6) * t64;
t393 = -t272 / 0.2e1;
t392 = t272 / 0.2e1;
t391 = -t273 / 0.2e1;
t369 = m(6) * qJD(5);
t354 = t176 * t281;
t353 = t178 * t281;
t241 = Icges(6,5) * t283 - Icges(6,6) * t281;
t350 = t241 * t273;
t242 = Icges(5,5) * t283 - Icges(5,6) * t281;
t349 = t242 * t273;
t343 = t273 * t283;
t172 = Icges(6,5) * t345 - Icges(6,6) * t346 - Icges(6,3) * t273;
t329 = -t272 * t172 - t180 * t343;
t173 = Icges(6,3) * t272 + t350;
t328 = t272 * t173 + t181 * t343;
t174 = Icges(5,5) * t345 - Icges(5,6) * t346 - Icges(5,3) * t273;
t327 = -t272 * t174 - t182 * t343;
t175 = Icges(5,3) * t272 + t349;
t326 = t272 * t175 + t183 * t343;
t306 = t177 * t281 - t172;
t156 = t181 * t345;
t308 = t173 * t273 - t156;
t86 = -t176 * t344 - t329;
t87 = -t177 * t344 + t328;
t10 = (t306 * t273 - t328 + t87) * t273 + (t306 * t272 + t308 + t86) * t272;
t305 = t179 * t281 - t174;
t157 = t183 * t345;
t307 = t175 * t273 - t157;
t88 = -t178 * t344 - t327;
t89 = -t179 * t344 + t326;
t11 = (t305 * t273 - t326 + t89) * t273 + (t305 * t272 + t307 + t88) * t272;
t83 = -t177 * t346 - t308;
t12 = (t83 - t156 + (t173 + t354) * t273 + t329) * t273 + t328 * t272;
t85 = -t179 * t346 - t307;
t13 = (t85 - t157 + (t175 + t353) * t273 + t327) * t273 + t326 * t272;
t48 = t272 * t83 - t273 * (-(-t180 * t283 + t354) * t272 - t172 * t273);
t49 = t272 * t85 - t273 * (-(-t182 * t283 + t353) * t272 - t174 * t273);
t50 = t272 * t87 - t273 * t86;
t51 = t272 * t89 - t273 * t88;
t2 = (t51 / 0.2e1 - t13 / 0.2e1 + t50 / 0.2e1 - t12 / 0.2e1) * t273 + (t11 / 0.2e1 + t49 / 0.2e1 + t10 / 0.2e1 + t48 / 0.2e1) * t272;
t61 = 0.2e1 * t121;
t313 = t2 * qJD(4) + t61 * qJD(5);
t310 = rSges(6,2) * t281 - t367 - t376;
t296 = (-t208 * t273 + t209 * t272) * t254;
t287 = t294 + t418;
t286 = -t294 + (t442 * t281 + t443 * t283) * (t392 + t393);
t60 = t121 + t110;
t285 = t60 * qJD(5) + ((t12 + t13) * t273 / 0.2e1 + (t10 + t11 + t48 + t49) * t393 + (t434 * t283 + t432 * t281 + t425 * t273 + (t241 + t242) * t272) * t392 + (t425 * t272 - t281 * t431 + t283 * t433 - t349 - t350 + t50 + t51) * t391) * qJD(4);
t257 = -rSges(5,2) * t281 + t368;
t189 = t310 * t273;
t187 = t310 * t272;
t144 = -t208 * t272 - t209 * t273;
t111 = t438 * t315;
t59 = 0.2e1 * t110;
t57 = t60 * qJD(4);
t55 = t59 * qJD(4);
t35 = t399 / 0.2e1;
t34 = t400 / 0.2e1;
t21 = t294 - t430;
t20 = t294 + t396 + t405;
t17 = t386 + t388 + t398 + t407;
t16 = t35 - t400 / 0.2e1;
t15 = t35 + t34;
t14 = t34 - t399 / 0.2e1;
t5 = t287 - t374;
t4 = t287 + t374;
t3 = t286 + t374 - t418;
t1 = [t17 * qJD(2) + t20 * qJD(4) + t71 * t369, t17 * qJD(1) + t4 * qJD(4) + t15 * qJD(5) + 0.2e1 * (t386 / 0.2e1 + t388 / 0.2e1 + t53 * t410 + t62 * t411) * qJD(2), 0, t20 * qJD(1) + t4 * qJD(2) + ((t138 * t189 + t139 * t187) * t410 + ((-t147 * t273 - t148 * t272) * t257 + t296) * t411) * t414 + t285, t15 * qJD(2) + t423 + t57; t5 * qJD(4) + t16 * qJD(5) + (-t398 / 0.4e1 - t407 / 0.4e1 - t386 / 0.4e1 - t388 / 0.4e1) * t416, t21 * qJD(4) + t75 * t369, 0, t5 * qJD(1) + t21 * qJD(2) + ((t142 * t189 + t143 * t187) * t410 + ((-t150 * t273 - t151 * t272) * t257 + t296) * t411) * t414 + t285, t16 * qJD(1) + t426 + t57; 0, 0, 0, (-t111 * t410 + t144 * t411) * t414, 0; t286 * qJD(1) + t3 * qJD(2) + (-t396 / 0.4e1 - t405 / 0.4e1) * t416 + t313, t3 * qJD(1) + t313 + (t286 + t430) * qJD(2), 0, (m(5) * (t254 * t257 * t315 + (t272 * (rSges(5,1) * t345 - t317) + t273 * (rSges(5,1) * t343 + t272 * rSges(5,3) - t240)) * t144) + m(6) * (-t111 * ((-pkin(3) * t272 + t268 + t295) * t272 + ((-pkin(3) + t309) * t273 + (rSges(6,3) - pkin(7)) * t272 + t316) * t273) - t446 * t187 - t445 * t189) + (t421 * t269 + (t419 * t273 + (t420 - t422) * t272) * t273) * t392 + (t422 * t270 + (t420 * t272 + (t419 - t421) * t273) * t272) * t391) * qJD(4) + t314 * t2, t314 * t61; t14 * qJD(2) - t423 + t55, t14 * qJD(1) - t426 + t55, 0, m(6) * (-t187 * t273 + t189 * t272) * qJD(4) + t314 * t59, 0;];
Cq = t1;
