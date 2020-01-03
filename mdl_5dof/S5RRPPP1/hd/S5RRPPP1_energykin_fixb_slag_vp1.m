% Calculate kinetic energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:19
% EndTime: 2019-12-31 19:23:20
% DurationCPUTime: 1.56s
% Computational Cost: add. (1078->188), mult. (2919->275), div. (0->0), fcn. (3339->8), ass. (0->106)
t462 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t461 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t460 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t459 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t458 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t457 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t456 = rSges(6,1) + pkin(4);
t455 = rSges(6,3) + qJ(5);
t399 = sin(pkin(5));
t401 = sin(qJ(1));
t400 = sin(qJ(2));
t402 = cos(qJ(2));
t441 = cos(pkin(8));
t442 = cos(pkin(5));
t419 = t442 * t441;
t440 = sin(pkin(8));
t405 = t400 * t419 + t440 * t402;
t403 = cos(qJ(1));
t422 = t403 * t441;
t360 = t399 * t422 + t405 * t401;
t418 = t442 * t440;
t407 = t400 * t418;
t423 = t401 * t441;
t424 = t399 * t440;
t361 = -t401 * t407 + t402 * t423 - t403 * t424;
t437 = t399 * t400;
t383 = t401 * t437 - t403 * t442;
t454 = t460 * t360 - t462 * t361 + t459 * t383;
t362 = -t399 * t423 + t405 * t403;
t363 = t401 * t424 + t402 * t422 - t403 * t407;
t384 = t401 * t442 + t403 * t437;
t453 = -t460 * t362 + t462 * t363 - t459 * t384;
t452 = -t457 * t360 + t460 * t361 - t458 * t383;
t451 = t457 * t362 - t460 * t363 + t458 * t384;
t450 = t458 * t360 - t459 * t361 + t461 * t383;
t449 = t458 * t362 - t459 * t363 + t461 * t384;
t381 = t400 * t440 - t402 * t419;
t382 = t400 * t441 + t402 * t418;
t436 = t399 * t402;
t448 = -t460 * t381 + t462 * t382 + t459 * t436;
t447 = t457 * t381 - t460 * t382 - t458 * t436;
t446 = t458 * t381 - t459 * t382 - t461 * t436;
t444 = pkin(2) * t402;
t439 = Icges(3,4) * t400;
t438 = Icges(3,4) * t402;
t435 = rSges(6,2) * t360 + t455 * t361 + t456 * t383;
t434 = rSges(6,2) * t362 + t455 * t363 + t456 * t384;
t385 = pkin(2) * t400 - qJ(3) * t436;
t433 = -pkin(3) * t382 - qJ(4) * t381 - t385;
t380 = qJD(3) * t384;
t432 = qJD(4) * t362 + t380;
t368 = t383 * qJ(3) + t401 * t444;
t398 = pkin(1) * t401 - pkin(7) * t403;
t431 = -t368 - t398;
t430 = qJD(2) * t401;
t429 = qJD(2) * t403;
t338 = pkin(3) * t361 + qJ(4) * t360;
t428 = -t338 + t431;
t369 = t384 * qJ(3) + t403 * t444;
t388 = qJD(1) * (pkin(1) * t403 + pkin(7) * t401);
t427 = qJD(1) * t369 + qJD(3) * t383 + t388;
t421 = qJD(2) * (-rSges(4,1) * t382 + rSges(4,2) * t381 + rSges(4,3) * t436 - t385);
t420 = qJD(2) * (rSges(5,1) * t436 + rSges(5,2) * t382 - rSges(5,3) * t381 + t433);
t417 = rSges(3,1) * t402 - rSges(3,2) * t400;
t339 = pkin(3) * t363 + qJ(4) * t362;
t416 = qJD(1) * t339 + qJD(4) * t360 + t427;
t415 = Icges(3,1) * t402 - t439;
t414 = -Icges(3,2) * t400 + t438;
t413 = Icges(3,5) * t402 - Icges(3,6) * t400;
t372 = -Icges(3,6) * t403 + t414 * t401;
t374 = -Icges(3,5) * t403 + t415 * t401;
t412 = t372 * t400 - t374 * t402;
t373 = Icges(3,6) * t401 + t414 * t403;
t375 = Icges(3,5) * t401 + t415 * t403;
t411 = -t373 * t400 + t375 * t402;
t390 = Icges(3,2) * t402 + t439;
t391 = Icges(3,1) * t400 + t438;
t410 = -t390 * t400 + t391 * t402;
t409 = qJD(2) * (-rSges(6,2) * t381 - t455 * t382 + t456 * t436 + t433);
t408 = -qJD(3) * t436 + t368 * t430 + t369 * t429;
t406 = qJD(4) * t381 + t338 * t430 + t339 * t429 + t408;
t394 = rSges(2,1) * t403 - rSges(2,2) * t401;
t393 = rSges(2,1) * t401 + rSges(2,2) * t403;
t392 = rSges(3,1) * t400 + rSges(3,2) * t402;
t389 = Icges(3,5) * t400 + Icges(3,6) * t402;
t377 = rSges(3,3) * t401 + t417 * t403;
t376 = -rSges(3,3) * t403 + t417 * t401;
t371 = Icges(3,3) * t401 + t413 * t403;
t370 = -Icges(3,3) * t403 + t413 * t401;
t344 = qJD(1) * t377 - t392 * t430 + t388;
t343 = -t392 * t429 + (-t376 - t398) * qJD(1);
t342 = (t376 * t401 + t377 * t403) * qJD(2);
t334 = rSges(5,1) * t384 - rSges(5,2) * t363 + rSges(5,3) * t362;
t332 = rSges(5,1) * t383 - rSges(5,2) * t361 + rSges(5,3) * t360;
t330 = rSges(4,1) * t363 - rSges(4,2) * t362 + rSges(4,3) * t384;
t329 = rSges(4,1) * t361 - rSges(4,2) * t360 + rSges(4,3) * t383;
t310 = qJD(1) * t330 + t401 * t421 + t427;
t309 = t380 + t403 * t421 + (-t329 + t431) * qJD(1);
t308 = (t329 * t401 + t330 * t403) * qJD(2) + t408;
t307 = qJD(1) * t334 + t401 * t420 + t416;
t306 = t403 * t420 + (-t332 + t428) * qJD(1) + t432;
t305 = (t332 * t401 + t334 * t403) * qJD(2) + t406;
t304 = t434 * qJD(1) + qJD(5) * t361 + t401 * t409 + t416;
t303 = qJD(5) * t363 + t403 * t409 + (t428 - t435) * qJD(1) + t432;
t302 = qJD(5) * t382 + (t435 * t401 + t434 * t403) * qJD(2) + t406;
t1 = m(3) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(4) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(5) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + m(6) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + (m(2) * (t393 ^ 2 + t394 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t402 * t372 - t400 * t374 + t452 * t381 + t382 * t454 + t450 * t436) * t403 + (t402 * t373 + t400 * t375 + t381 * t451 + t382 * t453 - t436 * t449) * t401) * qJD(2) + (t381 * t447 + t382 * t448 + t402 * t390 + t400 * t391 - t436 * t446) * qJD(1)) * qJD(1) / 0.2e1 + (((t452 * t362 + t454 * t363 - t450 * t384 + t412 * t403) * t403 + ((-t370 + t411) * t403 + t401 * t371 + t449 * t384 + t453 * t363 + t451 * t362) * t401) * qJD(2) + (t362 * t447 + t363 * t448 + t384 * t446 + t401 * t389 + t410 * t403) * qJD(1)) * t430 / 0.2e1 - (((t452 * t360 + t454 * t361 + t403 * t370 - t450 * t383) * t403 + (t411 * t401 + (-t371 + t412) * t403 + t449 * t383 + t453 * t361 + t451 * t360) * t401) * qJD(2) + (t360 * t447 + t361 * t448 + t383 * t446 - t403 * t389 + t410 * t401) * qJD(1)) * t429 / 0.2e1;
T = t1;
