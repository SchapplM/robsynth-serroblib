% Calculate kinetic energy for
% S5RRPRP5
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:53:59
% EndTime: 2019-12-31 19:54:01
% DurationCPUTime: 1.80s
% Computational Cost: add. (983->178), mult. (1060->267), div. (0->0), fcn. (901->8), ass. (0->113)
t461 = Icges(5,4) - Icges(6,5);
t460 = Icges(5,1) + Icges(6,1);
t459 = Icges(5,2) + Icges(6,3);
t362 = qJ(2) + pkin(8);
t358 = qJ(4) + t362;
t354 = cos(t358);
t458 = t461 * t354;
t353 = sin(t358);
t457 = t461 * t353;
t456 = Icges(6,4) + Icges(5,5);
t455 = Icges(5,6) - Icges(6,6);
t454 = t459 * t353 - t458;
t453 = t460 * t354 - t457;
t452 = rSges(6,1) + pkin(4);
t451 = rSges(6,3) + qJ(5);
t450 = Icges(6,2) + Icges(5,3);
t449 = Icges(3,3) + Icges(4,3);
t365 = sin(qJ(1));
t367 = cos(qJ(1));
t448 = -t454 * t365 - t455 * t367;
t447 = -t455 * t365 + t454 * t367;
t446 = t453 * t365 - t456 * t367;
t445 = t456 * t365 + t453 * t367;
t444 = -t459 * t354 - t457;
t443 = t460 * t353 + t458;
t442 = -t455 * t353 + t456 * t354;
t356 = sin(t362);
t357 = cos(t362);
t364 = sin(qJ(2));
t366 = cos(qJ(2));
t441 = Icges(3,5) * t366 + Icges(4,5) * t357 - Icges(3,6) * t364 - Icges(4,6) * t356;
t440 = t451 * t353 + t452 * t354;
t439 = t441 * t365 - t449 * t367;
t438 = t449 * t365 + t441 * t367;
t437 = Icges(3,5) * t364 + Icges(4,5) * t356 + Icges(3,6) * t366 + Icges(4,6) * t357;
t405 = qJD(2) + qJD(4);
t344 = t405 * t365;
t345 = t405 * t367;
t436 = (t448 * t353 - t446 * t354) * t345 + (t447 * t353 + t445 * t354) * t344 + (t444 * t353 + t443 * t354) * qJD(1);
t435 = (-t442 * t365 + t450 * t367) * t345 + (t450 * t365 + t442 * t367) * t344 + (t456 * t353 + t455 * t354) * qJD(1);
t420 = Icges(4,4) * t356;
t339 = Icges(4,2) * t357 + t420;
t419 = Icges(4,4) * t357;
t340 = Icges(4,1) * t356 + t419;
t422 = Icges(3,4) * t364;
t347 = Icges(3,2) * t366 + t422;
t421 = Icges(3,4) * t366;
t348 = Icges(3,1) * t364 + t421;
t434 = -t339 * t356 + t340 * t357 - t347 * t364 + t348 * t366;
t388 = -Icges(4,2) * t356 + t419;
t314 = Icges(4,6) * t365 + t367 * t388;
t392 = Icges(4,1) * t357 - t420;
t316 = Icges(4,5) * t365 + t367 * t392;
t389 = -Icges(3,2) * t364 + t421;
t324 = Icges(3,6) * t365 + t367 * t389;
t393 = Icges(3,1) * t366 - t422;
t326 = Icges(3,5) * t365 + t367 * t393;
t433 = -t314 * t356 + t316 * t357 - t324 * t364 + t326 * t366;
t313 = -Icges(4,6) * t367 + t365 * t388;
t315 = -Icges(4,5) * t367 + t365 * t392;
t323 = -Icges(3,6) * t367 + t365 * t389;
t325 = -Icges(3,5) * t367 + t365 * t393;
t432 = t313 * t356 - t315 * t357 + t323 * t364 - t325 * t366;
t426 = pkin(2) * t364;
t424 = t366 * pkin(2);
t309 = -qJ(3) * t367 + t365 * t424;
t310 = qJ(3) * t365 + t367 * t424;
t406 = qJD(2) * t367;
t407 = qJD(2) * t365;
t414 = t309 * t407 + t310 * t406;
t413 = -rSges(6,2) * t367 + t440 * t365;
t412 = rSges(6,2) * t365 + t440 * t367;
t352 = pkin(1) * t365 - pkin(6) * t367;
t411 = -t309 - t352;
t410 = t452 * t353 - t451 * t354;
t409 = pkin(3) * t357;
t288 = -pkin(7) * t367 + t365 * t409;
t404 = -t288 + t411;
t289 = pkin(7) * t365 + t367 * t409;
t401 = t288 * t407 + t289 * t406 + t414;
t343 = qJD(1) * (pkin(1) * t367 + pkin(6) * t365);
t400 = qJD(1) * t310 - qJD(3) * t367 + t343;
t399 = rSges(3,1) * t366 - rSges(3,2) * t364;
t398 = rSges(4,1) * t357 - rSges(4,2) * t356;
t397 = rSges(5,1) * t354 - rSges(5,2) * t353;
t394 = qJD(2) * (-rSges(4,1) * t356 - rSges(4,2) * t357 - t426);
t375 = qJD(1) * t289 + t400;
t374 = (-pkin(3) * t356 - t426) * qJD(2);
t371 = qJD(5) * t353 + t374;
t359 = qJD(3) * t365;
t351 = rSges(2,1) * t367 - rSges(2,2) * t365;
t350 = rSges(2,1) * t365 + rSges(2,2) * t367;
t349 = rSges(3,1) * t364 + rSges(3,2) * t366;
t337 = rSges(5,1) * t353 + rSges(5,2) * t354;
t328 = rSges(3,3) * t365 + t367 * t399;
t327 = -rSges(3,3) * t367 + t365 * t399;
t318 = rSges(4,3) * t365 + t367 * t398;
t317 = -rSges(4,3) * t367 + t365 * t398;
t308 = rSges(5,3) * t365 + t367 * t397;
t306 = -rSges(5,3) * t367 + t365 * t397;
t284 = qJD(1) * t328 - t349 * t407 + t343;
t283 = -t349 * t406 + (-t327 - t352) * qJD(1);
t282 = (t327 * t365 + t328 * t367) * qJD(2);
t281 = qJD(1) * t318 + t365 * t394 + t400;
t280 = t359 + t367 * t394 + (-t317 + t411) * qJD(1);
t279 = (t317 * t365 + t318 * t367) * qJD(2) + t414;
t278 = qJD(1) * t308 - t337 * t344 + t365 * t374 + t375;
t277 = -t337 * t345 + t359 + t367 * t374 + (-t306 + t404) * qJD(1);
t276 = t306 * t344 + t308 * t345 + t401;
t275 = qJD(1) * t412 - t344 * t410 + t365 * t371 + t375;
t274 = t359 - t410 * t345 + t371 * t367 + (t404 - t413) * qJD(1);
t273 = -qJD(5) * t354 + t344 * t413 + t345 * t412 + t401;
t1 = m(3) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + m(4) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(5) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + m(6) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + (t435 * t365 + t436 * t367) * t344 / 0.2e1 - (t436 * t365 - t435 * t367) * t345 / 0.2e1 + (m(2) * (t350 ^ 2 + t351 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t438 * t365 ^ 2 + (t432 * t367 + (t433 - t439) * t365) * t367) * qJD(2) + (t437 * t365 + t434 * t367) * qJD(1)) * t407 / 0.2e1 - ((t439 * t367 ^ 2 + (t433 * t365 + (t432 - t438) * t367) * t365) * qJD(2) + (t434 * t365 - t437 * t367) * qJD(1)) * t406 / 0.2e1 + (-(t446 * t353 + t448 * t354) * t345 + (t445 * t353 - t447 * t354) * t344 + ((-t313 * t357 - t356 * t315 - t323 * t366 - t325 * t364) * t367 + (t314 * t357 + t316 * t356 + t324 * t366 + t326 * t364) * t365) * qJD(2) + (t357 * t339 + t356 * t340 + t366 * t347 + t364 * t348 + t443 * t353 - t444 * t354) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
