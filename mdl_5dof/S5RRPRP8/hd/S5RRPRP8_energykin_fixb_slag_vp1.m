% Calculate kinetic energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:19
% EndTime: 2019-12-31 20:03:20
% DurationCPUTime: 1.92s
% Computational Cost: add. (627->174), mult. (1474->262), div. (0->0), fcn. (1471->6), ass. (0->108)
t455 = Icges(3,4) - Icges(4,5);
t454 = Icges(3,1) + Icges(4,1);
t453 = Icges(3,2) + Icges(4,3);
t361 = sin(qJ(2));
t452 = t455 * t361;
t364 = cos(qJ(2));
t451 = t455 * t364;
t450 = Icges(4,4) + Icges(3,5);
t449 = Icges(3,6) - Icges(4,6);
t448 = t453 * t361 - t451;
t447 = t454 * t364 - t452;
t446 = Icges(5,1) + Icges(6,1);
t445 = Icges(5,4) + Icges(6,4);
t444 = Icges(5,5) + Icges(6,5);
t443 = Icges(4,2) + Icges(3,3);
t442 = Icges(5,2) + Icges(6,2);
t441 = Icges(5,6) + Icges(6,6);
t440 = Icges(5,3) + Icges(6,3);
t362 = sin(qJ(1));
t365 = cos(qJ(1));
t439 = t448 * t362 + t449 * t365;
t438 = -t449 * t362 + t448 * t365;
t437 = -t447 * t362 + t450 * t365;
t436 = t450 * t362 + t447 * t365;
t435 = -t453 * t364 - t452;
t434 = t454 * t361 + t451;
t433 = -t449 * t361 + t450 * t364;
t432 = rSges(6,3) + qJ(5);
t363 = cos(qJ(4));
t360 = sin(qJ(4));
t402 = t364 * t360;
t338 = t361 * t363 - t402;
t327 = t338 * t362;
t403 = t361 * t360;
t371 = t364 * t363 + t403;
t328 = t371 * t362;
t431 = t441 * t327 + t444 * t328 + t440 * t365;
t329 = t338 * t365;
t330 = t371 * t365;
t430 = t441 * t329 + t444 * t330 - t440 * t362;
t429 = t442 * t327 + t445 * t328 + t441 * t365;
t428 = t442 * t329 + t445 * t330 - t441 * t362;
t427 = t445 * t327 + t446 * t328 + t444 * t365;
t426 = t445 * t329 + t446 * t330 - t444 * t362;
t425 = t444 * t338 - t441 * t371;
t424 = t445 * t338 - t442 * t371;
t423 = t446 * t338 - t445 * t371;
t422 = t433 * t362 - t443 * t365;
t421 = t443 * t362 + t433 * t365;
t420 = t450 * t361 + t449 * t364;
t419 = t435 * t361 + t434 * t364;
t418 = t438 * t361 + t436 * t364;
t417 = -t439 * t361 + t437 * t364;
t411 = pkin(3) * t364;
t409 = t363 * pkin(4);
t369 = pkin(4) * t403 + t409 * t364;
t401 = t328 * rSges(6,1) + t327 * rSges(6,2) + t369 * t362 + t432 * t365;
t400 = t330 * rSges(6,1) + t329 * rSges(6,2) - t432 * t362 + t369 * t365;
t399 = t338 * rSges(6,1) - rSges(6,2) * t371 - pkin(4) * t402 + t409 * t361;
t385 = pkin(2) * t364 + qJ(3) * t361;
t334 = t385 * t362;
t355 = t362 * pkin(1) - t365 * pkin(6);
t398 = -t334 - t355;
t397 = qJD(2) * t362;
t396 = qJD(2) * t365;
t395 = qJD(3) * t361;
t394 = qJD(2) - qJD(4);
t335 = t385 * t365;
t341 = qJD(1) * (t365 * pkin(1) + t362 * pkin(6));
t393 = qJD(1) * t335 + t362 * t395 + t341;
t339 = t365 * pkin(7) + t362 * t411;
t392 = -t339 + t398;
t350 = t361 * pkin(2) - t364 * qJ(3);
t389 = qJD(2) * (-t361 * rSges(4,1) + t364 * rSges(4,3) - t350);
t388 = -qJD(3) * t364 + t334 * t397 + t335 * t396;
t387 = rSges(3,1) * t364 - rSges(3,2) * t361;
t386 = rSges(4,1) * t364 + rSges(4,3) * t361;
t384 = qJD(2) * (-t361 * pkin(3) - t350);
t340 = -t362 * pkin(7) + t365 * t411;
t370 = t339 * t397 + t340 * t396 + t388;
t357 = t365 * t395;
t368 = t365 * t384 + t357;
t367 = qJD(1) * t340 + t362 * t384 + t393;
t354 = t365 * rSges(2,1) - t362 * rSges(2,2);
t353 = t362 * rSges(2,1) + t365 * rSges(2,2);
t352 = t361 * rSges(3,1) + t364 * rSges(3,2);
t343 = t394 * t365;
t342 = t394 * t362;
t326 = t362 * rSges(3,3) + t387 * t365;
t325 = t362 * rSges(4,2) + t386 * t365;
t324 = -t365 * rSges(3,3) + t387 * t362;
t323 = -t365 * rSges(4,2) + t386 * t362;
t307 = t338 * rSges(5,1) - rSges(5,2) * t371;
t297 = t330 * rSges(5,1) + t329 * rSges(5,2) - t362 * rSges(5,3);
t295 = t328 * rSges(5,1) + t327 * rSges(5,2) + t365 * rSges(5,3);
t281 = qJD(1) * t326 - t352 * t397 + t341;
t280 = -t352 * t396 + (-t324 - t355) * qJD(1);
t279 = (t324 * t362 + t326 * t365) * qJD(2);
t278 = qJD(1) * t325 + t362 * t389 + t393;
t277 = t357 + t365 * t389 + (-t323 + t398) * qJD(1);
t276 = (t323 * t362 + t325 * t365) * qJD(2) + t388;
t275 = qJD(1) * t297 - t342 * t307 + t367;
t274 = -t343 * t307 + (-t295 + t392) * qJD(1) + t368;
t273 = t342 * t295 + t343 * t297 + t370;
t272 = t400 * qJD(1) + qJD(5) * t365 - t399 * t342 + t367;
t271 = -qJD(5) * t362 - t399 * t343 + (t392 - t401) * qJD(1) + t368;
t270 = t401 * t342 + t400 * t343 + t370;
t1 = m(3) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(4) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + m(5) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(6) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + (-(t429 * t329 + t427 * t330 - t431 * t362) * t343 + (t428 * t329 + t426 * t330 - t430 * t362) * t342 + (t424 * t329 + t423 * t330 - t425 * t362) * qJD(1)) * t342 / 0.2e1 - (-(t429 * t327 + t427 * t328 + t431 * t365) * t343 + (t428 * t327 + t426 * t328 + t430 * t365) * t342 + (t424 * t327 + t423 * t328 + t425 * t365) * qJD(1)) * t343 / 0.2e1 + (m(2) * (t353 ^ 2 + t354 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t421 * t362 ^ 2 + (t417 * t365 + (t418 - t422) * t362) * t365) * qJD(2) + (t420 * t362 + t419 * t365) * qJD(1)) * t397 / 0.2e1 - ((t422 * t365 ^ 2 + (t418 * t362 + (t417 - t421) * t365) * t362) * qJD(2) + (t419 * t362 - t420 * t365) * qJD(1)) * t396 / 0.2e1 + (-(t427 * t338 - t371 * t429) * t343 + (t426 * t338 - t371 * t428) * t342 + ((t437 * t361 + t439 * t364) * t365 + (t436 * t361 - t438 * t364) * t362) * qJD(2) + (t423 * t338 + t434 * t361 - t435 * t364 - t424 * t371) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
