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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:22:00
% EndTime: 2019-12-05 18:22:11
% DurationCPUTime: 4.52s
% Computational Cost: add. (19343->291), mult. (12755->360), div. (0->0), fcn. (11248->8), ass. (0->190)
t454 = Icges(5,1) + Icges(6,1);
t453 = Icges(5,2) + Icges(6,2);
t284 = cos(qJ(4));
t278 = Icges(6,4) * t284;
t279 = Icges(5,4) * t284;
t282 = sin(qJ(4));
t441 = t454 * t282 + t278 + t279;
t280 = qJ(1) + qJ(2);
t275 = pkin(8) + t280;
t272 = sin(t275);
t273 = cos(t275);
t341 = t273 * t284;
t342 = t273 * t282;
t180 = Icges(6,4) * t341 - Icges(6,2) * t342 + Icges(6,6) * t272;
t182 = Icges(5,4) * t341 - Icges(5,2) * t342 + Icges(5,6) * t272;
t452 = t180 + t182;
t239 = Icges(6,4) * t342;
t184 = Icges(6,1) * t341 + Icges(6,5) * t272 - t239;
t241 = Icges(5,4) * t342;
t186 = Icges(5,1) * t341 + Icges(5,5) * t272 - t241;
t451 = t184 + t186;
t450 = (-Icges(5,4) - Icges(6,4)) * t282;
t449 = t453 * t284 - t450;
t448 = t454 * t284 + t450;
t447 = (Icges(5,6) + Icges(6,6)) * t284 + (Icges(5,5) + Icges(6,5)) * t282;
t446 = -t453 * t341 - t239 - t241 + t451;
t347 = t272 * t282;
t238 = Icges(6,4) * t347;
t346 = t272 * t284;
t183 = -Icges(6,1) * t346 + Icges(6,5) * t273 + t238;
t240 = Icges(5,4) * t347;
t185 = -Icges(5,1) * t346 + Icges(5,5) * t273 + t240;
t445 = t453 * t346 + t183 + t185 + t238 + t240;
t444 = t441 * t273 + t452;
t420 = Icges(6,2) * t282 - t278;
t179 = Icges(6,6) * t273 + t420 * t272;
t419 = Icges(5,2) * t282 - t279;
t181 = Icges(5,6) * t273 + t419 * t272;
t443 = t441 * t272 - t179 - t181;
t442 = t420 + t419;
t372 = pkin(4) * t284;
t274 = pkin(3) + t372;
t290 = rSges(6,1) * t341 - rSges(6,2) * t342 + t273 * t274;
t277 = cos(t280);
t374 = pkin(2) * t277;
t426 = -rSges(6,3) - qJ(5) - pkin(7);
t143 = t426 * t272 - t290 - t374;
t376 = cos(qJ(1)) * pkin(1);
t139 = t143 - t376;
t440 = t139 - t143;
t300 = -rSges(6,2) * t347 + t426 * t273;
t365 = rSges(6,1) * t284;
t304 = t274 + t365;
t276 = sin(t280);
t375 = pkin(2) * t276;
t142 = t304 * t272 + t300 + t375;
t248 = pkin(4) * t347;
t314 = rSges(6,1) * t347 + rSges(6,2) * t346;
t171 = t248 + t314;
t364 = rSges(6,2) * t284;
t172 = (t364 + (rSges(6,1) + pkin(4)) * t282) * t273;
t65 = -t142 * t171 + t143 * t172;
t247 = rSges(5,2) * t342;
t366 = rSges(5,1) * t284;
t306 = pkin(3) + t366;
t151 = -t374 + t247 - t306 * t273 + (-rSges(5,3) - pkin(7)) * t272;
t260 = rSges(5,1) * t282 + rSges(5,2) * t284;
t212 = t260 * t273;
t268 = t273 * pkin(7);
t315 = -rSges(5,2) * t347 - t273 * rSges(5,3);
t150 = t306 * t272 - t268 + t315 + t375;
t211 = t260 * t272;
t355 = t150 * t211;
t69 = t151 * t212 - t355;
t439 = -m(5) * t69 - m(6) * t65;
t405 = m(5) / 0.2e1;
t404 = m(6) / 0.2e1;
t377 = sin(qJ(1)) * pkin(1);
t384 = m(3) * (-t376 * (rSges(3,1) * t276 + rSges(3,2) * t277) - (-t277 * rSges(3,1) + t276 * rSges(3,2)) * t377);
t382 = m(4) * (-t376 * (rSges(4,1) * t272 + rSges(4,2) * t273 + t375) - (-t273 * rSges(4,1) + t272 * rSges(4,2) - t374) * t377);
t345 = t273 * t142;
t75 = -t143 * t272 - t345;
t434 = t75 * m(6) * qJD(2);
t433 = t447 * t272;
t432 = t447 * t273;
t430 = -t445 * t282 + t443 * t284;
t429 = t446 * t282 + t444 * t284;
t428 = (-t448 + t449) * t284 + (t441 - t442) * t282;
t308 = qJD(1) + qJD(2);
t138 = t142 + t377;
t129 = t273 * t138;
t71 = -t139 * t272 - t129;
t427 = t71 * m(6) * qJD(1);
t263 = -rSges(5,2) * t282 + t366;
t425 = t263 * t405;
t424 = (t182 * t282 - t186 * t284) * t273;
t423 = (t180 * t282 - t184 * t284) * t273;
t148 = t151 - t376;
t259 = rSges(6,1) * t282 + t364;
t189 = t259 * t272 + t248;
t373 = pkin(4) * t282;
t191 = (t259 + t373) * t273;
t147 = t150 + t377;
t356 = t147 * t211;
t371 = (t440 * t191 + (-t138 + t142) * t189) * t404 + (-t356 + t355 + (t148 - t151) * t212) * t405;
t64 = -t138 * t171 + t139 * t172;
t68 = t148 * t212 - t356;
t418 = (t65 + t64) * t404 + (t69 + t68) * t405;
t289 = (-t442 / 0.2e1 + t441 / 0.2e1) * t284 + (-t449 / 0.2e1 + t448 / 0.2e1) * t282;
t269 = t272 ^ 2;
t270 = t273 ^ 2;
t411 = 0.4e1 * qJD(1);
t408 = 2 * qJD(4);
t62 = -t151 * t147 + t148 * t150;
t401 = m(5) * t62;
t399 = m(5) * t68;
t395 = m(6) * (-t440 * t272 - t129 + t345);
t394 = m(6) * (-t129 - t345 + (-t139 - t143) * t272);
t53 = -t143 * t138 + t139 * t142;
t393 = m(6) * t53;
t391 = m(6) * t64;
t387 = t272 / 0.2e1;
t386 = -t273 / 0.2e1;
t385 = t273 / 0.2e1;
t380 = m(6) * (t171 * t272 + t172 * t273);
t378 = m(6) * (-t189 * t272 - t273 * t191);
t367 = m(6) * qJD(5);
t353 = t179 * t282;
t352 = t181 * t282;
t249 = Icges(6,5) * t284 - Icges(6,6) * t282;
t349 = t249 * t272;
t250 = Icges(5,5) * t284 - Icges(5,6) * t282;
t348 = t250 * t272;
t328 = -t191 * t171 + t189 * t172;
t175 = Icges(6,3) * t273 - t349;
t327 = t273 * t175 + t179 * t347;
t177 = Icges(5,3) * t273 - t348;
t326 = t273 * t177 + t181 * t347;
t325 = t272 * t175 + t183 * t341;
t324 = t272 * t177 + t185 * t341;
t309 = t269 + t270;
t176 = Icges(6,5) * t341 - Icges(6,6) * t342 + Icges(6,3) * t272;
t303 = -t183 * t284 - t176;
t83 = t273 * t176 + t180 * t347 - t184 * t346;
t86 = -t179 * t342 + t325;
t87 = t176 * t272 - t423;
t10 = (t327 + t87 + t423) * t273 + (-t86 + (t303 - t353) * t273 + t83 + t325) * t272;
t178 = Icges(5,5) * t341 - Icges(5,6) * t342 + Icges(5,3) * t272;
t302 = -t185 * t284 - t178;
t85 = t273 * t178 + t182 * t347 - t186 * t346;
t88 = -t181 * t342 + t324;
t89 = t178 * t272 - t424;
t11 = (t326 + t89 + t424) * t273 + (-t88 + (t302 - t352) * t273 + t85 + t324) * t272;
t82 = -t183 * t346 + t327;
t12 = (t83 + (-t176 + t353) * t273 - t325) * t273 + (t303 * t272 + t327 - t82) * t272;
t84 = -t185 * t346 + t326;
t13 = (t85 + (-t178 + t352) * t273 - t324) * t273 + (t302 * t272 + t326 - t84) * t272;
t48 = t272 * t83 + t273 * t82;
t49 = t272 * t85 + t273 * t84;
t50 = t272 * t87 + t273 * t86;
t51 = t272 * t89 + t273 * t88;
t2 = (t13 / 0.2e1 + t51 / 0.2e1 + t12 / 0.2e1 + t50 / 0.2e1) * t273 + (-t49 / 0.2e1 + t11 / 0.2e1 - t48 / 0.2e1 + t10 / 0.2e1) * t272;
t121 = t378 / 0.2e1;
t61 = t121 - t380 / 0.2e1;
t307 = t2 * qJD(4) + t61 * qJD(5);
t305 = -rSges(6,2) * t282 + t365 + t372;
t288 = t289 + t418;
t287 = -t289 + (t451 * t282 + t452 * t284) * (t385 + t386);
t110 = t380 / 0.2e1;
t60 = t121 + t110;
t286 = t60 * qJD(5) + (-(t10 + t11) * t272 / 0.2e1 + (t12 + t13 + t50 + t51) * t386 + (t445 * t284 + t443 * t282 + (t249 + t250) * t273 + t428 * t272) * t385 + (-t428 * t273 - t444 * t282 + t446 * t284 + t348 + t349 + t48 + t49) * t387) * qJD(4);
t192 = t305 * t273;
t190 = t305 * t272;
t144 = -t211 * t272 - t212 * t273;
t111 = -t259 * t270 - t272 * t314 - t309 * t373;
t59 = t110 - t378 / 0.2e1;
t57 = t60 * qJD(4);
t55 = t59 * qJD(4);
t35 = t394 / 0.2e1;
t34 = t395 / 0.2e1;
t21 = t289 - t439;
t20 = t289 + t391 + t399;
t17 = t382 + t384 + t393 + t401;
t16 = t35 - t395 / 0.2e1;
t15 = t35 + t34;
t14 = t34 - t394 / 0.2e1;
t5 = t288 + t371;
t4 = t288 - t371;
t3 = t287 + t371 - t418;
t1 = [t17 * qJD(2) + t20 * qJD(4) + t71 * t367, t17 * qJD(1) + t5 * qJD(4) + t15 * qJD(5) + 0.2e1 * (t382 / 0.2e1 + t384 / 0.2e1 + t53 * t404 + t62 * t405) * qJD(2), 0, t20 * qJD(1) + t5 * qJD(2) + ((t147 * t273 + t148 * t272) * t425 + (t138 * t192 + t139 * t190 + t328) * t404) * t408 + t286, t15 * qJD(2) + t427 + t57; t4 * qJD(4) + t16 * qJD(5) + (-t384 / 0.4e1 - t382 / 0.4e1 - t401 / 0.4e1 - t393 / 0.4e1) * t411, t21 * qJD(4) + t75 * t367, 0, t4 * qJD(1) + t21 * qJD(2) + ((t142 * t192 + t143 * t190 + t328) * t404 + (t150 * t273 + t151 * t272) * t425) * t408 + t286, t16 * qJD(1) + t434 + t57; 0, 0, 0, (t111 * t404 + t144 * t405) * t408, 0; t287 * qJD(1) + t3 * qJD(2) + (-t399 / 0.4e1 - t391 / 0.4e1) * t411 + t307, t3 * qJD(1) + t307 + (t287 + t439) * qJD(2), 0, (m(5) * (t309 * t263 * t260 + (t273 * (rSges(5,1) * t341 + t272 * rSges(5,3) - t247) - t272 * (-rSges(5,1) * t346 - t315)) * t144) + m(6) * (t111 * ((-t273 * pkin(3) + t290) * t273 + ((-pkin(7) - t426) * t273 + t268 + t300 + (-pkin(3) + t304) * t272) * t272) + t189 * t190 + t191 * t192) + (-t432 * t269 + (t430 * t273 + (-t429 + t433) * t272) * t273) * t387 + (t433 * t270 + (t429 * t272 + (-t430 - t432) * t273) * t272) * t385) * qJD(4) + t308 * t2, t308 * t61; t14 * qJD(2) - t427 + t55, t14 * qJD(1) - t434 + t55, 0, m(6) * (t190 * t273 - t192 * t272) * qJD(4) + t308 * t59, 0;];
Cq = t1;
