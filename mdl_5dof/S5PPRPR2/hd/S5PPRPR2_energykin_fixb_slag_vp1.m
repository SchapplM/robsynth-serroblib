% Calculate kinetic energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:58
% EndTime: 2019-12-05 15:02:59
% DurationCPUTime: 0.83s
% Computational Cost: add. (549->108), mult. (744->180), div. (0->0), fcn. (718->6), ass. (0->65)
t368 = Icges(5,1) + Icges(4,3);
t318 = pkin(8) + qJ(3);
t316 = sin(t318);
t317 = cos(t318);
t367 = (-Icges(5,4) + Icges(4,5)) * t317 + (Icges(5,5) - Icges(4,6)) * t316;
t320 = cos(pkin(7));
t359 = t320 ^ 2;
t319 = sin(pkin(7));
t360 = t319 ^ 2;
t343 = t359 + t360;
t361 = qJD(3) * t343;
t366 = t319 * t320;
t365 = t367 * t319 - t368 * t320;
t364 = t368 * t319 + t367 * t320;
t357 = t317 * t319;
t356 = t317 * t320;
t321 = sin(qJ(5));
t355 = t319 * t321;
t322 = cos(qJ(5));
t354 = t319 * t322;
t353 = t320 * t321;
t352 = t320 * t322;
t315 = qJD(2) * t319;
t347 = qJD(4) * t316;
t351 = t320 * t347 + t315;
t350 = qJD(2) * t320;
t349 = qJD(3) * t319;
t348 = qJD(3) * t320;
t346 = qJD(5) * t316;
t345 = qJD(5) * t317;
t344 = qJD(1) + (pkin(3) * t317 + qJ(4) * t316) * t361;
t310 = t316 * pkin(3) - t317 * qJ(4);
t340 = qJD(3) * (t316 * rSges(5,2) + t317 * rSges(5,3) - t310);
t339 = t319 * t347 - t350;
t338 = qJD(3) * (-pkin(6) * t316 - t310);
t324 = qJD(1) ^ 2;
t312 = t316 * rSges(4,1) + t317 * rSges(4,2);
t309 = t319 * t345 - t348;
t308 = t320 * t345 + t349;
t307 = t316 * t355 - t352;
t306 = t316 * t354 + t353;
t305 = t316 * t353 + t354;
t304 = t316 * t352 - t355;
t303 = -t312 * t349 - t350;
t302 = -t312 * t348 + t315;
t287 = t316 * rSges(6,3) + (-rSges(6,1) * t321 - rSges(6,2) * t322) * t317;
t286 = Icges(6,5) * t316 + (-Icges(6,1) * t321 - Icges(6,4) * t322) * t317;
t285 = Icges(6,6) * t316 + (-Icges(6,4) * t321 - Icges(6,2) * t322) * t317;
t284 = Icges(6,3) * t316 + (-Icges(6,5) * t321 - Icges(6,6) * t322) * t317;
t283 = t307 * rSges(6,1) + t306 * rSges(6,2) + rSges(6,3) * t357;
t282 = t305 * rSges(6,1) + t304 * rSges(6,2) + rSges(6,3) * t356;
t281 = Icges(6,1) * t307 + Icges(6,4) * t306 + Icges(6,5) * t357;
t280 = Icges(6,1) * t305 + Icges(6,4) * t304 + Icges(6,5) * t356;
t279 = Icges(6,4) * t307 + Icges(6,2) * t306 + Icges(6,6) * t357;
t278 = Icges(6,4) * t305 + Icges(6,2) * t304 + Icges(6,6) * t356;
t277 = Icges(6,5) * t307 + Icges(6,6) * t306 + Icges(6,3) * t357;
t276 = Icges(6,5) * t305 + Icges(6,6) * t304 + Icges(6,3) * t356;
t275 = t319 * t340 + t339;
t274 = t320 * t340 + t351;
t273 = qJD(1) + (rSges(4,1) * t317 - rSges(4,2) * t316) * t361;
t272 = -qJD(4) * t317 + t344 + (-rSges(5,2) * t317 + rSges(5,3) * t316) * t361;
t271 = t282 * t346 - t308 * t287 + t319 * t338 + t339;
t270 = -t283 * t346 + t309 * t287 + t320 * t338 + t351;
t269 = -t309 * t282 + t308 * t283 + (pkin(6) * t361 - qJD(4)) * t317 + t344;
t1 = m(2) * t324 / 0.2e1 + m(3) * (t343 * qJD(2) ^ 2 + t324) / 0.2e1 + m(4) * (t273 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(5) * (t272 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(6) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + t308 * ((t276 * t356 + t304 * t278 + t305 * t280) * t308 + (t277 * t356 + t304 * t279 + t305 * t281) * t309 + (t284 * t356 + t304 * t285 + t305 * t286) * t346) / 0.2e1 + t309 * ((t276 * t357 + t306 * t278 + t307 * t280) * t308 + (t277 * t357 + t306 * t279 + t307 * t281) * t309 + (t284 * t357 + t306 * t285 + t307 * t286) * t346) / 0.2e1 + ((t276 * t308 + t277 * t309 + t284 * t346) * t316 + ((-t278 * t322 - t280 * t321) * t308 + (-t279 * t322 - t281 * t321) * t309 + (-t285 * t322 - t286 * t321) * t346) * t317) * t346 / 0.2e1 + ((t364 * t360 - t365 * t366) * t319 / 0.2e1 - (t365 * t359 - t364 * t366) * t320 / 0.2e1) * qJD(3) ^ 2;
T = t1;
