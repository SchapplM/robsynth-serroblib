% Calculate kinetic energy for
% S5RPRRP3
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:31
% EndTime: 2022-01-23 09:29:33
% DurationCPUTime: 1.27s
% Computational Cost: add. (935->139), mult. (779->215), div. (0->0), fcn. (648->8), ass. (0->91)
t390 = Icges(5,4) + Icges(6,4);
t389 = Icges(5,1) + Icges(6,1);
t388 = Icges(5,2) + Icges(6,2);
t315 = qJ(3) + qJ(4);
t311 = cos(t315);
t387 = t390 * t311;
t310 = sin(t315);
t386 = t390 * t310;
t385 = Icges(5,5) + Icges(6,5);
t384 = Icges(5,6) + Icges(6,6);
t383 = -t388 * t310 + t387;
t382 = t389 * t311 - t386;
t381 = rSges(6,1) + pkin(4);
t380 = Icges(5,3) + Icges(6,3);
t314 = qJ(1) + pkin(8);
t308 = sin(t314);
t309 = cos(t314);
t379 = t383 * t308 - t384 * t309;
t378 = t384 * t308 + t383 * t309;
t377 = t382 * t308 - t385 * t309;
t376 = t385 * t308 + t382 * t309;
t375 = t388 * t311 + t386;
t374 = t389 * t310 + t387;
t373 = -t384 * t310 + t385 * t311;
t372 = rSges(6,3) + qJ(5);
t371 = -rSges(6,2) * t310 + t381 * t311;
t348 = qJD(3) + qJD(4);
t288 = t348 * t308;
t289 = t348 * t309;
t370 = (t379 * t310 - t377 * t311) * t289 + (-t378 * t310 + t376 * t311) * t288 + (-t375 * t310 + t374 * t311) * qJD(1);
t369 = (-t373 * t308 + t380 * t309) * t289 + (t380 * t308 + t373 * t309) * t288 + (t385 * t310 + t384 * t311) * qJD(1);
t317 = sin(qJ(1));
t365 = pkin(1) * t317;
t318 = cos(qJ(3));
t363 = t318 * pkin(3);
t316 = sin(qJ(3));
t361 = Icges(4,4) * t316;
t360 = Icges(4,4) * t318;
t355 = t371 * t308 - t372 * t309;
t354 = t372 * t308 + t371 * t309;
t319 = cos(qJ(1));
t307 = qJD(1) * t319 * pkin(1);
t353 = qJD(1) * (pkin(2) * t309 + pkin(6) * t308) + t307;
t350 = qJD(3) * t308;
t349 = qJD(3) * t309;
t347 = pkin(3) * qJD(3) * t316;
t259 = -pkin(7) * t309 + t363 * t308;
t260 = pkin(7) * t308 + t363 * t309;
t346 = t259 * t350 + t260 * t349 + qJD(2);
t345 = -pkin(2) * t308 + pkin(6) * t309 - t365;
t344 = rSges(6,2) * t311 + t381 * t310;
t343 = t309 * t347;
t342 = -t259 + t345;
t341 = rSges(4,1) * t318 - rSges(4,2) * t316;
t340 = rSges(5,1) * t311 - rSges(5,2) * t310;
t338 = Icges(4,1) * t318 - t361;
t335 = -Icges(4,2) * t316 + t360;
t332 = Icges(4,5) * t318 - Icges(4,6) * t316;
t279 = -Icges(4,6) * t309 + t335 * t308;
t281 = -Icges(4,5) * t309 + t338 * t308;
t329 = t279 * t316 - t281 * t318;
t280 = Icges(4,6) * t308 + t335 * t309;
t282 = Icges(4,5) * t308 + t338 * t309;
t328 = -t280 * t316 + t282 * t318;
t301 = Icges(4,2) * t318 + t361;
t302 = Icges(4,1) * t316 + t360;
t327 = -t301 * t316 + t302 * t318;
t326 = qJD(1) * t260 - t308 * t347 + t353;
t305 = rSges(2,1) * t319 - rSges(2,2) * t317;
t304 = rSges(2,1) * t317 + rSges(2,2) * t319;
t303 = rSges(4,1) * t316 + rSges(4,2) * t318;
t300 = Icges(4,5) * t316 + Icges(4,6) * t318;
t298 = rSges(5,1) * t310 + rSges(5,2) * t311;
t286 = t307 + qJD(1) * (rSges(3,1) * t309 - rSges(3,2) * t308);
t285 = (-rSges(3,1) * t308 - rSges(3,2) * t309 - t365) * qJD(1);
t284 = rSges(4,3) * t308 + t341 * t309;
t283 = -rSges(4,3) * t309 + t341 * t308;
t278 = Icges(4,3) * t308 + t332 * t309;
t277 = -Icges(4,3) * t309 + t332 * t308;
t276 = rSges(5,3) * t308 + t340 * t309;
t274 = -rSges(5,3) * t309 + t340 * t308;
t253 = qJD(1) * t284 - t303 * t350 + t353;
t252 = -t303 * t349 + (-t283 + t345) * qJD(1);
t251 = qJD(2) + (t283 * t308 + t284 * t309) * qJD(3);
t250 = qJD(1) * t276 - t288 * t298 + t326;
t249 = -t343 - t289 * t298 + (-t274 + t342) * qJD(1);
t248 = t274 * t288 + t276 * t289 + t346;
t247 = t354 * qJD(1) - qJD(5) * t309 - t344 * t288 + t326;
t246 = -t343 + qJD(5) * t308 - t344 * t289 + (t342 - t355) * qJD(1);
t245 = t355 * t288 + t354 * t289 + t346;
t1 = m(3) * (qJD(2) ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + m(4) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + ((t308 * t300 + t327 * t309) * qJD(1) + (t308 ^ 2 * t278 + (t329 * t309 + (-t277 + t328) * t308) * t309) * qJD(3)) * t350 / 0.2e1 - ((-t309 * t300 + t327 * t308) * qJD(1) + (t309 ^ 2 * t277 + (t328 * t308 + (-t278 + t329) * t309) * t308) * qJD(3)) * t349 / 0.2e1 + m(5) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(6) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + (t369 * t308 + t370 * t309) * t288 / 0.2e1 - (t370 * t308 - t369 * t309) * t289 / 0.2e1 + (m(2) * (t304 ^ 2 + t305 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + (((t280 * t318 + t282 * t316) * t308 - (t279 * t318 + t316 * t281) * t309) * qJD(3) - (t377 * t310 + t379 * t311) * t289 + (t376 * t310 + t378 * t311) * t288 + (t318 * t301 + t316 * t302 + t374 * t310 + t375 * t311) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
