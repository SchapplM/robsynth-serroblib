% Calculate kinetic energy for
% S5RPRRP7
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:32
% EndTime: 2019-12-31 18:44:33
% DurationCPUTime: 1.21s
% Computational Cost: add. (1000->164), mult. (1223->260), div. (0->0), fcn. (1214->8), ass. (0->92)
t370 = Icges(5,1) + Icges(6,1);
t369 = -Icges(5,4) + Icges(6,5);
t368 = Icges(6,4) + Icges(5,5);
t367 = Icges(5,2) + Icges(6,3);
t366 = -Icges(6,6) + Icges(5,6);
t365 = -Icges(5,3) - Icges(6,2);
t364 = rSges(6,1) + pkin(4);
t363 = rSges(6,3) + qJ(5);
t316 = qJ(1) + pkin(8);
t314 = sin(t316);
t315 = cos(t316);
t320 = cos(qJ(4));
t317 = sin(qJ(4));
t321 = cos(qJ(3));
t344 = t317 * t321;
t285 = t314 * t344 + t315 * t320;
t343 = t320 * t321;
t286 = t314 * t343 - t315 * t317;
t318 = sin(qJ(3));
t346 = t314 * t318;
t362 = t367 * t285 + t369 * t286 - t366 * t346;
t287 = -t314 * t320 + t315 * t344;
t288 = t314 * t317 + t315 * t343;
t345 = t315 * t318;
t361 = t367 * t287 + t369 * t288 - t366 * t345;
t360 = -t366 * t285 + t368 * t286 - t365 * t346;
t359 = -t366 * t287 + t368 * t288 - t365 * t345;
t358 = t369 * t285 + t370 * t286 + t368 * t346;
t357 = t369 * t287 + t370 * t288 + t368 * t345;
t356 = t366 * t321 + (t367 * t317 + t369 * t320) * t318;
t355 = t365 * t321 + (-t366 * t317 + t368 * t320) * t318;
t354 = -t368 * t321 + (t369 * t317 + t370 * t320) * t318;
t319 = sin(qJ(1));
t349 = t319 * pkin(1);
t348 = Icges(4,4) * t318;
t347 = Icges(4,4) * t321;
t342 = rSges(6,2) * t346 + t363 * t285 + t364 * t286;
t341 = rSges(6,2) * t345 + t363 * t287 + t364 * t288;
t340 = -t321 * rSges(6,2) + (t363 * t317 + t364 * t320) * t318;
t322 = cos(qJ(1));
t313 = qJD(1) * t322 * pkin(1);
t339 = qJD(1) * (t315 * pkin(2) + t314 * pkin(6)) + t313;
t338 = qJD(3) * t314;
t337 = qJD(3) * t315;
t336 = qJD(4) * t318;
t333 = pkin(3) * t321 + pkin(7) * t318;
t297 = t333 * t314;
t298 = t333 * t315;
t335 = t297 * t338 + t298 * t337 + qJD(2);
t334 = -t314 * pkin(2) + t315 * pkin(6) - t349;
t332 = rSges(4,1) * t321 - rSges(4,2) * t318;
t331 = Icges(4,1) * t321 - t348;
t330 = -Icges(4,2) * t318 + t347;
t329 = Icges(4,5) * t321 - Icges(4,6) * t318;
t274 = -Icges(4,6) * t315 + t314 * t330;
t276 = -Icges(4,5) * t315 + t314 * t331;
t328 = t274 * t318 - t276 * t321;
t275 = Icges(4,6) * t314 + t315 * t330;
t277 = Icges(4,5) * t314 + t315 * t331;
t327 = -t275 * t318 + t277 * t321;
t305 = Icges(4,2) * t321 + t348;
t306 = Icges(4,1) * t318 + t347;
t326 = -t305 * t318 + t306 * t321;
t310 = t318 * pkin(3) - t321 * pkin(7);
t325 = qJD(1) * t298 - t310 * t338 + t339;
t324 = (-t297 + t334) * qJD(1) - t310 * t337;
t312 = -qJD(4) * t321 + qJD(1);
t309 = t322 * rSges(2,1) - t319 * rSges(2,2);
t308 = t319 * rSges(2,1) + t322 * rSges(2,2);
t307 = t318 * rSges(4,1) + t321 * rSges(4,2);
t304 = Icges(4,5) * t318 + Icges(4,6) * t321;
t300 = t314 * t336 - t337;
t299 = t315 * t336 + t338;
t296 = -t321 * rSges(5,3) + (rSges(5,1) * t320 - rSges(5,2) * t317) * t318;
t283 = t313 + qJD(1) * (t315 * rSges(3,1) - t314 * rSges(3,2));
t282 = (-t314 * rSges(3,1) - t315 * rSges(3,2) - t349) * qJD(1);
t279 = t314 * rSges(4,3) + t315 * t332;
t278 = -t315 * rSges(4,3) + t314 * t332;
t273 = Icges(4,3) * t314 + t315 * t329;
t272 = -Icges(4,3) * t315 + t314 * t329;
t269 = t288 * rSges(5,1) - t287 * rSges(5,2) + rSges(5,3) * t345;
t267 = t286 * rSges(5,1) - t285 * rSges(5,2) + rSges(5,3) * t346;
t253 = qJD(1) * t279 - t307 * t338 + t339;
t252 = -t307 * t337 + (-t278 + t334) * qJD(1);
t251 = qJD(2) + (t278 * t314 + t279 * t315) * qJD(3);
t250 = t312 * t269 - t299 * t296 + t325;
t249 = -t312 * t267 + t300 * t296 + t324;
t248 = t299 * t267 - t300 * t269 + t335;
t247 = qJD(5) * t285 - t299 * t340 + t312 * t341 + t325;
t246 = qJD(5) * t287 + t300 * t340 - t312 * t342 + t324;
t245 = qJD(5) * t318 * t317 + t299 * t342 - t300 * t341 + t335;
t1 = m(3) * (qJD(2) ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(4) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + ((t314 * t304 + t315 * t326) * qJD(1) + (t314 ^ 2 * t273 + (t328 * t315 + (-t272 + t327) * t314) * t315) * qJD(3)) * t338 / 0.2e1 - ((-t315 * t304 + t314 * t326) * qJD(1) + (t315 ^ 2 * t272 + (t327 * t314 + (-t273 + t328) * t315) * t314) * qJD(3)) * t337 / 0.2e1 + qJD(1) * ((t321 * t305 + t318 * t306) * qJD(1) + ((t321 * t275 + t318 * t277) * t314 - (t321 * t274 + t318 * t276) * t315) * qJD(3)) / 0.2e1 + m(5) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(6) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + ((t356 * t287 + t354 * t288 + t355 * t345) * t312 + (t362 * t287 + t358 * t288 + t360 * t345) * t300 + (t361 * t287 + t357 * t288 + t359 * t345) * t299) * t299 / 0.2e1 + ((t356 * t285 + t354 * t286 + t355 * t346) * t312 + (t362 * t285 + t358 * t286 + t360 * t346) * t300 + (t361 * t285 + t357 * t286 + t359 * t346) * t299) * t300 / 0.2e1 + ((-t359 * t299 - t360 * t300 - t355 * t312) * t321 + ((t356 * t317 + t354 * t320) * t312 + (t362 * t317 + t358 * t320) * t300 + (t361 * t317 + t357 * t320) * t299) * t318) * t312 / 0.2e1 + (m(2) * (t308 ^ 2 + t309 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
