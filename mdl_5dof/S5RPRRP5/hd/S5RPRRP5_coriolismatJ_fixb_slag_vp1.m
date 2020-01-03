% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:49
% DurationCPUTime: 4.29s
% Computational Cost: add. (21502->293), mult. (14724->377), div. (0->0), fcn. (13143->8), ass. (0->197)
t290 = qJ(1) + pkin(8);
t287 = qJ(3) + t290;
t283 = cos(t287);
t291 = sin(qJ(4));
t361 = t283 * t291;
t282 = sin(t287);
t293 = cos(qJ(4));
t392 = rSges(6,1) + pkin(4);
t323 = t392 * t291;
t381 = rSges(6,3) + qJ(5);
t453 = -t381 * t293 + t323;
t458 = t453 * t282;
t459 = t361 * t458;
t288 = Icges(6,5) * t291;
t311 = Icges(6,3) * t293 - t288;
t378 = Icges(5,4) * t291;
t456 = Icges(5,2) * t293 + t311 + t378;
t289 = Icges(5,4) * t293;
t260 = Icges(5,1) * t291 + t289;
t377 = Icges(6,5) * t293;
t457 = Icges(6,1) * t291 + t260 - t377;
t257 = -Icges(5,2) * t291 + t289;
t455 = t257 + t457;
t261 = Icges(5,1) * t293 - t378;
t430 = Icges(6,1) * t293 + t288;
t454 = t261 + t430;
t280 = t282 ^ 2;
t281 = t283 ^ 2;
t327 = t280 + t281;
t431 = t381 * t291 + t392 * t293;
t279 = t283 * pkin(7);
t426 = pkin(3) + t431;
t145 = t283 * rSges(6,2) - t426 * t282 + t279;
t319 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t290);
t128 = t145 + t319;
t451 = t128 - t145;
t449 = (-Icges(5,6) + Icges(6,6)) * t293 + (-Icges(6,4) - Icges(5,5)) * t291;
t196 = Icges(6,4) * t282 + t283 * t430;
t198 = Icges(5,5) * t282 + t261 * t283;
t448 = -t456 * t283 + t196 + t198;
t195 = -Icges(6,4) * t283 + t282 * t430;
t363 = t282 * t291;
t248 = Icges(5,4) * t363;
t362 = t282 * t293;
t197 = Icges(5,1) * t362 - Icges(5,5) * t283 - t248;
t447 = -Icges(5,2) * t362 - t311 * t282 + t195 + t197 - t248;
t360 = t283 * t293;
t247 = Icges(6,5) * t360;
t188 = Icges(6,6) * t282 + Icges(6,3) * t361 + t247;
t194 = Icges(5,6) * t282 + t257 * t283;
t446 = -Icges(6,1) * t361 - t260 * t283 + t188 - t194 + t247;
t253 = Icges(6,3) * t291 + t377;
t187 = -Icges(6,6) * t283 + t253 * t282;
t193 = Icges(5,4) * t362 - Icges(5,2) * t363 - Icges(5,6) * t283;
t445 = t457 * t282 - t187 + t193;
t146 = (rSges(6,2) + pkin(7)) * t282 + t426 * t283;
t163 = -t283 * t323 + t381 * t360;
t67 = t145 * t458 + t146 * t163;
t382 = rSges(5,1) * t293;
t322 = pkin(3) + t382;
t330 = rSges(5,2) * t363 + t283 * rSges(5,3);
t158 = -t322 * t282 + t279 + t330;
t250 = rSges(5,2) * t361;
t159 = -t250 + t322 * t283 + (rSges(5,3) + pkin(7)) * t282;
t264 = rSges(5,1) * t291 + rSges(5,2) * t293;
t218 = t264 * t282;
t220 = t264 * t283;
t84 = t158 * t218 - t159 * t220;
t444 = -m(5) * t84 - m(6) * t67;
t421 = m(5) / 0.2e1;
t420 = m(6) / 0.2e1;
t443 = -t291 / 0.2e1;
t318 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t290);
t390 = m(4) * (t318 * (-rSges(4,1) * t282 - rSges(4,2) * t283) - (t283 * rSges(4,1) - t282 * rSges(4,2)) * t319);
t324 = t146 * t361;
t81 = -t145 * t363 + t324;
t440 = t81 * m(6) * qJD(3);
t255 = Icges(6,4) * t293 + Icges(6,6) * t291;
t367 = t255 * t282;
t191 = -Icges(6,2) * t283 + t367;
t178 = t282 * t191;
t102 = t187 * t361 + t195 * t360 + t178;
t439 = t102 * t283;
t438 = (t454 - t456) * t293 + (t253 - t455) * t291;
t435 = (t187 * t291 + t195 * t293) * t282;
t154 = t158 + t319;
t155 = t159 + t318;
t75 = -t159 * t154 + t155 * t158;
t433 = t449 * t282;
t432 = t449 * t283;
t429 = -t448 * t291 + t446 * t293;
t428 = t447 * t291 + t445 * t293;
t129 = t146 + t318;
t184 = t453 * t283;
t388 = ((-t129 + t146) * t184 + t451 * t458) * t420 + ((-t155 + t159) * t283 + (t154 - t158) * t282) * t264 * t421;
t61 = t128 * t458 + t129 * t163;
t78 = t154 * t218 - t155 * t220;
t427 = (t67 + t61) * t420 + (t84 + t78) * t421;
t304 = t456 * t443 + t454 * t291 / 0.2e1 + (-t253 / 0.2e1 + t455 / 0.2e1) * t293;
t425 = 4 * qJD(1);
t423 = 2 * qJD(4);
t417 = m(5) * t75;
t415 = m(5) * t78;
t124 = t129 * t361;
t409 = m(6) * (-t451 * t363 + t124 - t324);
t408 = m(6) * (t124 + t324 + (-t128 - t145) * t363);
t42 = -t146 * t128 + t129 * t145;
t407 = m(6) * t42;
t343 = t163 * t363 + t459;
t346 = t128 * t360 + t129 * t362;
t406 = m(6) * (t343 + t346);
t317 = t184 * t363 - t459;
t404 = m(6) * (t317 + t346);
t344 = t145 * t360 + t146 * t362;
t403 = m(6) * (t343 + t344);
t402 = m(6) * t61;
t401 = m(6) * (t317 + t344);
t221 = t327 * t291;
t397 = t221 / 0.2e1;
t396 = -t282 / 0.2e1;
t395 = t282 / 0.2e1;
t394 = -t283 / 0.2e1;
t366 = t255 * t283;
t192 = Icges(6,2) * t282 + t366;
t103 = t188 * t361 + t282 * t192 + t196 * t360;
t371 = t191 * t283;
t306 = t103 + t371;
t316 = -t188 * t363 + t192 * t283 - t196 * t362;
t16 = (t103 - t306) * t283 + (t102 + t316 - t178) * t282;
t189 = Icges(5,5) * t362 - Icges(5,6) * t363 - Icges(5,3) * t283;
t341 = -t282 * t189 - t197 * t360;
t104 = -t193 * t361 - t341;
t254 = Icges(5,5) * t293 - Icges(5,6) * t291;
t368 = t254 * t283;
t190 = Icges(5,3) * t282 + t368;
t340 = t282 * t190 + t198 * t360;
t105 = -t194 * t361 + t340;
t320 = t194 * t291 - t189;
t169 = t198 * t362;
t321 = t190 * t283 - t169;
t17 = (t320 * t283 + t105 - t340) * t283 + (t320 * t282 + t104 + t321) * t282;
t98 = -t371 + t435;
t18 = -t439 + (t306 + t98 - t435) * t282;
t101 = -t194 * t363 - t321;
t370 = t193 * t291;
t19 = (t101 - t169 + (t190 + t370) * t283 + t341) * t283 + t340 * t282;
t56 = -t282 * t316 - t283 * t98;
t57 = -(-(-t197 * t293 + t370) * t282 - t189 * t283) * t283 + t101 * t282;
t58 = t103 * t282 - t439;
t59 = -t104 * t283 + t105 * t282;
t2 = (t59 / 0.2e1 - t19 / 0.2e1 + t58 / 0.2e1 - t18 / 0.2e1) * t283 + (t17 / 0.2e1 + t57 / 0.2e1 + t16 / 0.2e1 + t56 / 0.2e1) * t282;
t391 = t2 * qJD(4);
t384 = m(6) * qJD(4);
t383 = m(6) * qJD(5);
t355 = t291 * t293;
t345 = (-t163 - t184) * t458;
t342 = -t184 * t360 - t362 * t458;
t331 = t327 * t355;
t202 = (t397 + t443) * m(6);
t326 = t202 * qJD(2);
t77 = -t128 * t363 + t124;
t325 = m(6) * t77 * qJD(1);
t305 = (-t218 * t283 + t220 * t282) * t264;
t297 = t304 + t427;
t296 = -t304 + ((t187 + t193) * t293 + (-t195 + t197) * t291) * (t395 + t396);
t295 = ((t18 + t19) * t283 / 0.2e1 + (t16 + t17 + t56 + t57) * t396 + (t254 * t282 + t438 * t283 + t446 * t291 + t448 * t293 + t367) * t395 + (t438 * t282 - t445 * t291 + t447 * t293 - t366 - t368 + t58 + t59) * t394) * qJD(4);
t268 = -rSges(5,2) * t291 + t382;
t201 = m(6) * t397 + t291 * t420;
t186 = t331 - t355;
t185 = t431 * t283;
t183 = t431 * t282;
t151 = -t218 * t282 - t220 * t283;
t106 = t163 * t283 - t280 * t453;
t88 = t327 * t431;
t63 = t221 * t88 + t342;
t60 = t401 / 0.2e1;
t49 = t403 / 0.2e1;
t48 = t404 / 0.2e1;
t43 = t406 / 0.2e1;
t37 = t408 / 0.2e1;
t36 = t409 / 0.2e1;
t27 = t304 - t444;
t26 = t304 + t402 + t415;
t23 = t390 + t407 + t417;
t22 = t60 - t403 / 0.2e1;
t21 = t60 + t49;
t20 = t49 - t401 / 0.2e1;
t11 = t48 - t406 / 0.2e1;
t10 = t48 + t43;
t9 = t43 - t404 / 0.2e1;
t8 = t37 - t409 / 0.2e1;
t7 = t37 + t36;
t6 = t36 - t408 / 0.2e1;
t5 = t297 + t388;
t4 = t297 - t388;
t3 = t296 + t388 - t427;
t1 = [t23 * qJD(3) + t26 * qJD(4) + t77 * t383, 0, t23 * qJD(1) + t5 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t390 / 0.2e1 + t42 * t420 + t75 * t421) * qJD(3), t26 * qJD(1) + t5 * qJD(3) + t10 * qJD(5) + ((-t128 * t185 - t129 * t183 + t345) * t420 + ((-t154 * t283 - t155 * t282) * t268 + t305) * t421) * t423 + t295, t7 * qJD(3) + t10 * qJD(4) + t325; 0, 0, 0, (t106 * t420 + t151 * t421) * t423 + t201 * qJD(5), t201 * qJD(4); t4 * qJD(4) + t8 * qJD(5) + (-t407 / 0.4e1 - t417 / 0.4e1 - t390 / 0.4e1) * t425, 0, t27 * qJD(4) + t81 * t383, t4 * qJD(1) + t27 * qJD(3) + t21 * qJD(5) + ((-t145 * t185 - t146 * t183 + t345) * t420 + ((-t158 * t283 - t159 * t282) * t268 + t305) * t421) * t423 + t295, t8 * qJD(1) + t21 * qJD(4) + t440; t296 * qJD(1) + t3 * qJD(3) + t11 * qJD(5) + (-t402 / 0.4e1 - t415 / 0.4e1) * t425 + t391, qJD(5) * t202, t3 * qJD(1) + t22 * qJD(5) + t391 + (t296 + t444) * qJD(3), (m(5) * (t327 * t268 * t264 + (t282 * (rSges(5,1) * t362 - t330) + t283 * (rSges(5,1) * t360 + t282 * rSges(5,3) - t250)) * t151) + m(6) * (t106 * t88 + t183 * t458 + t184 * t185) + (t432 * t280 + (t428 * t283 + (t429 - t433) * t282) * t283) * t395 + (t433 * t281 + (t429 * t282 + (t428 - t432) * t283) * t282) * t394) * qJD(4) + t63 * t383 + (qJD(1) + qJD(3)) * t2, t11 * qJD(1) + t326 + t22 * qJD(3) + t63 * t384 + (-t221 * t293 - t186 + t331) * t383; t6 * qJD(3) + t9 * qJD(4) - t325, -t202 * qJD(4), t6 * qJD(1) + t20 * qJD(4) - t440, t9 * qJD(1) - t326 + t20 * qJD(3) + (-t106 * t293 + (-t183 * t282 - t185 * t283 + t88) * t291 - t63 + t342) * t384 + t186 * t383, t186 * t384;];
Cq = t1;
