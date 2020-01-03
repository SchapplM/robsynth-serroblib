% Calculate kinetic energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:37
% EndTime: 2019-12-31 18:08:38
% DurationCPUTime: 1.20s
% Computational Cost: add. (859->135), mult. (748->197), div. (0->0), fcn. (617->8), ass. (0->85)
t393 = Icges(5,4) - Icges(6,5);
t392 = Icges(5,1) + Icges(6,1);
t391 = Icges(5,2) + Icges(6,3);
t311 = qJ(3) + pkin(8);
t309 = cos(t311);
t390 = t393 * t309;
t307 = sin(t311);
t389 = t393 * t307;
t388 = Icges(6,4) + Icges(5,5);
t387 = Icges(5,6) - Icges(6,6);
t386 = t391 * t307 - t390;
t385 = t392 * t309 - t389;
t384 = rSges(6,1) + pkin(4);
t383 = rSges(6,3) + qJ(5);
t312 = qJ(1) + pkin(7);
t308 = sin(t312);
t310 = cos(t312);
t382 = t386 * t308 + t387 * t310;
t381 = -t387 * t308 + t386 * t310;
t380 = -t385 * t308 + t388 * t310;
t379 = t388 * t308 + t385 * t310;
t378 = -t391 * t309 - t389;
t377 = t392 * t307 + t390;
t376 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t314 = sin(qJ(3));
t316 = cos(qJ(3));
t375 = Icges(4,5) * t316 - Icges(4,6) * t314 - t387 * t307 + t388 * t309;
t374 = t383 * t307 + t384 * t309;
t373 = t375 * t308 - t376 * t310;
t372 = t376 * t308 + t375 * t310;
t371 = Icges(4,5) * t314 + Icges(4,6) * t316 + t388 * t307 + t387 * t309;
t359 = Icges(4,4) * t314;
t299 = Icges(4,2) * t316 + t359;
t358 = Icges(4,4) * t316;
t300 = Icges(4,1) * t314 + t358;
t370 = -t299 * t314 + t300 * t316 + t378 * t307 + t377 * t309;
t335 = -Icges(4,2) * t314 + t358;
t278 = Icges(4,6) * t308 + t335 * t310;
t338 = Icges(4,1) * t316 - t359;
t280 = Icges(4,5) * t308 + t338 * t310;
t369 = -t278 * t314 + t280 * t316 + t381 * t307 + t379 * t309;
t277 = -Icges(4,6) * t310 + t335 * t308;
t279 = -Icges(4,5) * t310 + t338 * t308;
t368 = t277 * t314 - t279 * t316 - t382 * t307 + t380 * t309;
t315 = sin(qJ(1));
t364 = pkin(1) * t315;
t363 = pkin(3) * t314;
t361 = pkin(3) * t316;
t353 = -rSges(6,2) * t310 + t374 * t308;
t352 = rSges(6,2) * t308 + t374 * t310;
t317 = cos(qJ(1));
t306 = qJD(1) * t317 * pkin(1);
t351 = qJD(1) * (pkin(2) * t310 + pkin(6) * t308) + t306;
t350 = qJD(3) * t308;
t349 = qJD(3) * t310;
t257 = -qJ(4) * t310 + t361 * t308;
t258 = qJ(4) * t308 + t361 * t310;
t348 = t257 * t350 + t258 * t349 + qJD(2);
t345 = -pkin(2) * t308 + pkin(6) * t310 - t364;
t344 = -t257 + t345;
t343 = rSges(4,1) * t316 - rSges(4,2) * t314;
t342 = rSges(5,1) * t309 - rSges(5,2) * t307;
t339 = qJD(3) * (-rSges(5,1) * t307 - rSges(5,2) * t309 - t363);
t320 = qJD(1) * t258 - qJD(4) * t310 + t351;
t319 = (t383 * t309 - t363) * qJD(3) + (-t384 * qJD(3) + qJD(5)) * t307;
t304 = qJD(4) * t308;
t303 = rSges(2,1) * t317 - rSges(2,2) * t315;
t302 = rSges(2,1) * t315 + rSges(2,2) * t317;
t301 = rSges(4,1) * t314 + rSges(4,2) * t316;
t286 = t306 + qJD(1) * (rSges(3,1) * t310 - rSges(3,2) * t308);
t285 = (-rSges(3,1) * t308 - rSges(3,2) * t310 - t364) * qJD(1);
t282 = rSges(4,3) * t308 + t343 * t310;
t281 = -rSges(4,3) * t310 + t343 * t308;
t274 = rSges(5,3) * t308 + t342 * t310;
t272 = -rSges(5,3) * t310 + t342 * t308;
t253 = qJD(1) * t282 - t301 * t350 + t351;
t252 = -t301 * t349 + (-t281 + t345) * qJD(1);
t251 = qJD(2) + (t281 * t308 + t282 * t310) * qJD(3);
t250 = qJD(1) * t274 + t308 * t339 + t320;
t249 = t304 + t310 * t339 + (-t272 + t344) * qJD(1);
t248 = (t272 * t308 + t274 * t310) * qJD(3) + t348;
t247 = t352 * qJD(1) + t319 * t308 + t320;
t246 = t304 + t319 * t310 + (t344 - t353) * qJD(1);
t245 = -qJD(5) * t309 + (t353 * t308 + t352 * t310) * qJD(3) + t348;
t1 = m(3) * (qJD(2) ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + m(4) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + m(5) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(6) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + (m(2) * (t302 ^ 2 + t303 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t277 * t316 - t279 * t314 + t380 * t307 + t382 * t309) * t310 + (t278 * t316 + t280 * t314 + t379 * t307 - t381 * t309) * t308) * qJD(3) + (t316 * t299 + t314 * t300 + t377 * t307 - t378 * t309) * qJD(1)) * qJD(1) / 0.2e1 + ((t372 * t308 ^ 2 + (t368 * t310 + (t369 - t373) * t308) * t310) * qJD(3) + (t308 * t371 + t310 * t370) * qJD(1)) * t350 / 0.2e1 - ((t373 * t310 ^ 2 + (t369 * t308 + (t368 - t372) * t310) * t308) * qJD(3) + (t308 * t370 - t310 * t371) * qJD(1)) * t349 / 0.2e1;
T = t1;
