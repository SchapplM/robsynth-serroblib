% Calculate kinetic energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:50
% EndTime: 2019-12-31 20:54:52
% DurationCPUTime: 1.83s
% Computational Cost: add. (1001->180), mult. (1078->264), div. (0->0), fcn. (919->8), ass. (0->114)
t445 = Icges(5,4) - Icges(6,5);
t444 = Icges(5,1) + Icges(6,1);
t443 = Icges(5,2) + Icges(6,3);
t356 = qJ(2) + qJ(3);
t350 = pkin(8) + t356;
t348 = cos(t350);
t442 = t445 * t348;
t347 = sin(t350);
t441 = t445 * t347;
t440 = Icges(6,4) + Icges(5,5);
t439 = Icges(5,6) - Icges(6,6);
t438 = t443 * t347 - t442;
t437 = t444 * t348 - t441;
t436 = rSges(6,1) + pkin(4);
t435 = rSges(6,3) + qJ(5);
t358 = sin(qJ(1));
t360 = cos(qJ(1));
t434 = -t438 * t358 - t439 * t360;
t433 = -t439 * t358 + t438 * t360;
t432 = t437 * t358 - t440 * t360;
t431 = t440 * t358 + t437 * t360;
t430 = -t443 * t348 - t441;
t429 = t444 * t347 + t442;
t428 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t352 = sin(t356);
t353 = cos(t356);
t427 = Icges(4,5) * t353 - Icges(4,6) * t352 - t439 * t347 + t440 * t348;
t426 = t435 * t347 + t436 * t348;
t400 = pkin(3) * t353;
t282 = qJ(4) * t358 + t400 * t360;
t425 = qJD(1) * t282 - qJD(4) * t360;
t412 = Icges(4,4) * t353;
t380 = -Icges(4,2) * t352 + t412;
t306 = -Icges(4,6) * t360 + t380 * t358;
t307 = Icges(4,6) * t358 + t380 * t360;
t413 = Icges(4,4) * t352;
t384 = Icges(4,1) * t353 - t413;
t308 = -Icges(4,5) * t360 + t384 * t358;
t309 = Icges(4,5) * t358 + t384 * t360;
t333 = Icges(4,2) * t353 + t413;
t334 = Icges(4,1) * t352 + t412;
t395 = qJD(2) + qJD(3);
t338 = t395 * t358;
t339 = t395 * t360;
t424 = (t306 * t352 - t308 * t353 + t434 * t347 - t432 * t348) * t339 + (-t307 * t352 + t309 * t353 + t433 * t347 + t431 * t348) * t338 + (-t333 * t352 + t334 * t353 + t430 * t347 + t429 * t348) * qJD(1);
t423 = (-t427 * t358 + t428 * t360) * t339 + (t428 * t358 + t427 * t360) * t338 + (Icges(4,5) * t352 + Icges(4,6) * t353 + t440 * t347 + t439 * t348) * qJD(1);
t419 = pkin(3) * t352;
t359 = cos(qJ(2));
t417 = t359 * pkin(2);
t357 = sin(qJ(2));
t415 = Icges(3,4) * t357;
t414 = Icges(3,4) * t359;
t302 = -pkin(7) * t360 + t417 * t358;
t303 = pkin(7) * t358 + t417 * t360;
t397 = qJD(2) * t360;
t398 = qJD(2) * t358;
t407 = t302 * t398 + t303 * t397;
t406 = -rSges(6,2) * t360 + t426 * t358;
t405 = rSges(6,2) * t358 + t426 * t360;
t337 = qJD(1) * (pkin(1) * t360 + pkin(6) * t358);
t404 = qJD(1) * t303 + t337;
t346 = pkin(1) * t358 - pkin(6) * t360;
t403 = -t302 - t346;
t402 = t436 * t347 - t435 * t348;
t401 = qJD(4) * t358 - t339 * t419;
t394 = pkin(2) * qJD(2) * t357;
t281 = -qJ(4) * t360 + t400 * t358;
t393 = t338 * t281 + t407;
t392 = -t281 + t403;
t391 = t360 * t394;
t390 = rSges(3,1) * t359 - rSges(3,2) * t357;
t389 = rSges(4,1) * t353 - rSges(4,2) * t352;
t388 = rSges(5,1) * t348 - rSges(5,2) * t347;
t385 = Icges(3,1) * t359 - t415;
t381 = -Icges(3,2) * t357 + t414;
t377 = Icges(3,5) * t359 - Icges(3,6) * t357;
t316 = -Icges(3,6) * t360 + t381 * t358;
t318 = -Icges(3,5) * t360 + t385 * t358;
t373 = t316 * t357 - t318 * t359;
t317 = Icges(3,6) * t358 + t381 * t360;
t319 = Icges(3,5) * t358 + t385 * t360;
t372 = -t317 * t357 + t319 * t359;
t341 = Icges(3,2) * t359 + t415;
t342 = Icges(3,1) * t357 + t414;
t371 = -t341 * t357 + t342 * t359;
t370 = qJD(5) * t347 - t394;
t369 = -t358 * t394 + t404;
t345 = rSges(2,1) * t360 - rSges(2,2) * t358;
t344 = rSges(2,1) * t358 + rSges(2,2) * t360;
t343 = rSges(3,1) * t357 + rSges(3,2) * t359;
t340 = Icges(3,5) * t357 + Icges(3,6) * t359;
t335 = rSges(4,1) * t352 + rSges(4,2) * t353;
t330 = rSges(5,1) * t347 + rSges(5,2) * t348;
t321 = rSges(3,3) * t358 + t390 * t360;
t320 = -rSges(3,3) * t360 + t390 * t358;
t315 = Icges(3,3) * t358 + t377 * t360;
t314 = -Icges(3,3) * t360 + t377 * t358;
t311 = rSges(4,3) * t358 + t389 * t360;
t310 = -rSges(4,3) * t360 + t389 * t358;
t300 = rSges(5,3) * t358 + t388 * t360;
t298 = -rSges(5,3) * t360 + t388 * t358;
t279 = qJD(1) * t321 - t343 * t398 + t337;
t278 = -t343 * t397 + (-t320 - t346) * qJD(1);
t277 = (t320 * t358 + t321 * t360) * qJD(2);
t275 = qJD(1) * t311 - t335 * t338 + t369;
t274 = -t391 - t335 * t339 + (-t310 + t403) * qJD(1);
t273 = t310 * t338 + t311 * t339 + t407;
t272 = qJD(1) * t300 + (-t330 - t419) * t338 + t369 + t425;
t271 = -t391 - t330 * t339 + (-t298 + t392) * qJD(1) + t401;
t270 = t370 * t358 + t405 * qJD(1) + (-t402 - t419) * t338 + t404 + t425;
t269 = t370 * t360 - t402 * t339 + (t392 - t406) * qJD(1) + t401;
t268 = t298 * t338 - (-t282 - t300) * t339 + t393;
t267 = -qJD(5) * t348 + t406 * t338 - (-t282 - t405) * t339 + t393;
t1 = m(3) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + ((t358 * t340 + t371 * t360) * qJD(1) + (t358 ^ 2 * t315 + (t373 * t360 + (-t314 + t372) * t358) * t360) * qJD(2)) * t398 / 0.2e1 - ((-t360 * t340 + t371 * t358) * qJD(1) + (t360 ^ 2 * t314 + (t372 * t358 + (-t315 + t373) * t360) * t358) * qJD(2)) * t397 / 0.2e1 + m(4) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(5) * (t268 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + m(6) * (t267 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + (m(2) * (t344 ^ 2 + t345 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (t423 * t358 + t424 * t360) * t338 / 0.2e1 - (t424 * t358 - t423 * t360) * t339 / 0.2e1 + (((t317 * t359 + t319 * t357) * t358 - (t316 * t359 + t318 * t357) * t360) * qJD(2) - (t306 * t353 + t308 * t352 + t432 * t347 + t434 * t348) * t339 + (t307 * t353 + t309 * t352 + t431 * t347 - t433 * t348) * t338 + (t353 * t333 + t352 * t334 + t359 * t341 + t357 * t342 + t429 * t347 - t430 * t348) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
