% Calculate kinetic energy for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:29
% EndTime: 2019-12-05 16:43:30
% DurationCPUTime: 1.24s
% Computational Cost: add. (923->131), mult. (753->207), div. (0->0), fcn. (636->6), ass. (0->87)
t383 = Icges(5,4) + Icges(6,4);
t382 = Icges(5,1) + Icges(6,1);
t381 = Icges(5,2) + Icges(6,2);
t313 = qJ(3) + qJ(4);
t309 = cos(t313);
t380 = t383 * t309;
t308 = sin(t313);
t379 = t383 * t308;
t378 = Icges(5,5) + Icges(6,5);
t377 = Icges(5,6) + Icges(6,6);
t376 = -t381 * t308 + t380;
t375 = t382 * t309 - t379;
t374 = rSges(6,1) + pkin(4);
t373 = Icges(5,3) + Icges(6,3);
t312 = pkin(8) + qJ(2);
t306 = sin(t312);
t307 = cos(t312);
t372 = t376 * t306 - t377 * t307;
t371 = t377 * t306 + t376 * t307;
t370 = t375 * t306 - t378 * t307;
t369 = t378 * t306 + t375 * t307;
t368 = t381 * t309 + t379;
t367 = t382 * t308 + t380;
t366 = -t377 * t308 + t378 * t309;
t365 = rSges(6,3) + qJ(5);
t364 = -rSges(6,2) * t308 + t374 * t309;
t343 = qJD(3) + qJD(4);
t287 = t343 * t306;
t288 = t343 * t307;
t363 = (t372 * t308 - t370 * t309) * t288 + (-t371 * t308 + t369 * t309) * t287 + (-t368 * t308 + t367 * t309) * qJD(2);
t362 = (-t366 * t306 + t373 * t307) * t288 + (t373 * t306 + t366 * t307) * t287 + (t378 * t308 + t377 * t309) * qJD(2);
t315 = cos(qJ(3));
t358 = t315 * pkin(3);
t314 = sin(qJ(3));
t356 = Icges(4,4) * t314;
t355 = Icges(4,4) * t315;
t350 = t364 * t306 - t365 * t307;
t349 = t365 * t306 + t364 * t307;
t260 = -pkin(7) * t307 + t358 * t306;
t291 = pkin(2) * t306 - pkin(6) * t307;
t348 = -t260 - t291;
t345 = qJD(3) * t306;
t344 = qJD(3) * t307;
t342 = pkin(3) * qJD(3) * t314;
t261 = pkin(7) * t306 + t358 * t307;
t341 = t260 * t345 + t261 * t344 + qJD(1);
t340 = rSges(6,2) * t309 + t374 * t308;
t339 = t307 * t342;
t338 = rSges(4,1) * t315 - rSges(4,2) * t314;
t337 = rSges(5,1) * t309 - rSges(5,2) * t308;
t335 = Icges(4,1) * t315 - t356;
t332 = -Icges(4,2) * t314 + t355;
t329 = Icges(4,5) * t315 - Icges(4,6) * t314;
t280 = -Icges(4,6) * t307 + t332 * t306;
t282 = -Icges(4,5) * t307 + t335 * t306;
t326 = t280 * t314 - t282 * t315;
t281 = Icges(4,6) * t306 + t332 * t307;
t283 = Icges(4,5) * t306 + t335 * t307;
t325 = -t281 * t314 + t283 * t315;
t302 = Icges(4,2) * t315 + t356;
t303 = Icges(4,1) * t314 + t355;
t324 = -t302 * t314 + t303 * t315;
t286 = qJD(2) * (pkin(2) * t307 + pkin(6) * t306);
t323 = qJD(2) * t261 - t306 * t342 + t286;
t318 = qJD(1) ^ 2;
t317 = qJD(2) ^ 2;
t304 = rSges(4,1) * t314 + rSges(4,2) * t315;
t301 = Icges(4,5) * t314 + Icges(4,6) * t315;
t299 = rSges(5,1) * t308 + rSges(5,2) * t309;
t290 = rSges(3,1) * t307 - rSges(3,2) * t306;
t289 = rSges(3,1) * t306 + rSges(3,2) * t307;
t285 = rSges(4,3) * t306 + t338 * t307;
t284 = -rSges(4,3) * t307 + t338 * t306;
t279 = Icges(4,3) * t306 + t329 * t307;
t278 = -Icges(4,3) * t307 + t329 * t306;
t277 = rSges(5,3) * t306 + t337 * t307;
t275 = -rSges(5,3) * t307 + t337 * t306;
t254 = qJD(2) * t285 - t304 * t345 + t286;
t253 = -t304 * t344 + (-t284 - t291) * qJD(2);
t252 = qJD(1) + (t284 * t306 + t285 * t307) * qJD(3);
t251 = qJD(2) * t277 - t287 * t299 + t323;
t250 = -t339 - t288 * t299 + (-t275 + t348) * qJD(2);
t249 = t275 * t287 + t277 * t288 + t341;
t248 = t349 * qJD(2) - qJD(5) * t307 - t340 * t287 + t323;
t247 = -t339 + qJD(5) * t306 - t340 * t288 + (t348 - t350) * qJD(2);
t246 = t350 * t287 + t349 * t288 + t341;
t1 = m(2) * t318 / 0.2e1 + m(3) * (t318 + (t289 ^ 2 + t290 ^ 2) * t317) / 0.2e1 + t317 * Icges(3,3) / 0.2e1 + m(4) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + ((t306 * t301 + t324 * t307) * qJD(2) + (t306 ^ 2 * t279 + (t326 * t307 + (-t278 + t325) * t306) * t307) * qJD(3)) * t345 / 0.2e1 - ((-t307 * t301 + t324 * t306) * qJD(2) + (t307 ^ 2 * t278 + (t325 * t306 + (-t279 + t326) * t307) * t306) * qJD(3)) * t344 / 0.2e1 + m(5) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(6) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + (t362 * t306 + t363 * t307) * t287 / 0.2e1 - (t363 * t306 - t362 * t307) * t288 / 0.2e1 + (((t281 * t315 + t283 * t314) * t306 - (t280 * t315 + t282 * t314) * t307) * qJD(3) - (t370 * t308 + t372 * t309) * t288 + (t369 * t308 + t371 * t309) * t287 + (t315 * t302 + t314 * t303 + t367 * t308 + t368 * t309) * qJD(2)) * qJD(2) / 0.2e1;
T = t1;
