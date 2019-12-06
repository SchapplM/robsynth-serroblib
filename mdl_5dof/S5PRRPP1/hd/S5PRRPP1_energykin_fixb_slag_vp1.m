% Calculate kinetic energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:10
% EndTime: 2019-12-05 16:06:11
% DurationCPUTime: 1.19s
% Computational Cost: add. (847->127), mult. (722->189), div. (0->0), fcn. (605->6), ass. (0->81)
t386 = Icges(5,4) - Icges(6,5);
t385 = Icges(5,1) + Icges(6,1);
t384 = Icges(5,2) + Icges(6,3);
t310 = qJ(3) + pkin(8);
t308 = cos(t310);
t383 = t386 * t308;
t306 = sin(t310);
t382 = t386 * t306;
t381 = Icges(6,4) + Icges(5,5);
t380 = Icges(5,6) - Icges(6,6);
t379 = t384 * t306 - t383;
t378 = t385 * t308 - t382;
t377 = rSges(6,1) + pkin(4);
t376 = rSges(6,3) + qJ(5);
t309 = pkin(7) + qJ(2);
t305 = sin(t309);
t307 = cos(t309);
t375 = t379 * t305 + t380 * t307;
t374 = -t380 * t305 + t379 * t307;
t373 = -t378 * t305 + t381 * t307;
t372 = t381 * t305 + t378 * t307;
t371 = -t384 * t308 - t382;
t370 = t385 * t306 + t383;
t369 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t312 = sin(qJ(3));
t313 = cos(qJ(3));
t368 = Icges(4,5) * t313 - Icges(4,6) * t312 - t380 * t306 + t381 * t308;
t367 = t376 * t306 + t377 * t308;
t366 = t368 * t305 - t369 * t307;
t365 = t369 * t305 + t368 * t307;
t364 = Icges(4,5) * t312 + Icges(4,6) * t313 + t381 * t306 + t380 * t308;
t354 = Icges(4,4) * t312;
t300 = Icges(4,2) * t313 + t354;
t353 = Icges(4,4) * t313;
t301 = Icges(4,1) * t312 + t353;
t363 = -t300 * t312 + t301 * t313 + t371 * t306 + t370 * t308;
t331 = -Icges(4,2) * t312 + t353;
t279 = Icges(4,6) * t305 + t331 * t307;
t334 = Icges(4,1) * t313 - t354;
t281 = Icges(4,5) * t305 + t334 * t307;
t362 = -t279 * t312 + t281 * t313 + t374 * t306 + t372 * t308;
t278 = -Icges(4,6) * t307 + t331 * t305;
t280 = -Icges(4,5) * t307 + t334 * t305;
t361 = t278 * t312 - t280 * t313 - t375 * t306 + t373 * t308;
t358 = pkin(3) * t312;
t356 = pkin(3) * t313;
t258 = -qJ(4) * t307 + t356 * t305;
t298 = pkin(2) * t305 - pkin(6) * t307;
t348 = -t258 - t298;
t347 = -t307 * rSges(6,2) + t367 * t305;
t346 = t305 * rSges(6,2) + t367 * t307;
t345 = qJD(3) * t305;
t344 = qJD(3) * t307;
t259 = qJ(4) * t305 + t356 * t307;
t343 = t258 * t345 + t259 * t344 + qJD(1);
t286 = qJD(2) * (pkin(2) * t307 + pkin(6) * t305);
t340 = qJD(2) * t259 - qJD(4) * t307 + t286;
t339 = rSges(4,1) * t313 - rSges(4,2) * t312;
t338 = rSges(5,1) * t308 - rSges(5,2) * t306;
t335 = qJD(3) * (-rSges(5,1) * t306 - rSges(5,2) * t308 - t358);
t316 = (t376 * t308 - t358) * qJD(3) + (-t377 * qJD(3) + qJD(5)) * t306;
t315 = qJD(1) ^ 2;
t314 = qJD(2) ^ 2;
t303 = qJD(4) * t305;
t302 = rSges(4,1) * t312 + rSges(4,2) * t313;
t297 = rSges(3,1) * t307 - rSges(3,2) * t305;
t293 = rSges(3,1) * t305 + rSges(3,2) * t307;
t283 = t305 * rSges(4,3) + t339 * t307;
t282 = -t307 * rSges(4,3) + t339 * t305;
t275 = t305 * rSges(5,3) + t338 * t307;
t273 = -t307 * rSges(5,3) + t338 * t305;
t254 = qJD(2) * t283 - t302 * t345 + t286;
t253 = -t302 * t344 + (-t282 - t298) * qJD(2);
t252 = qJD(1) + (t282 * t305 + t283 * t307) * qJD(3);
t251 = qJD(2) * t275 + t305 * t335 + t340;
t250 = t303 + t307 * t335 + (-t273 + t348) * qJD(2);
t249 = (t273 * t305 + t275 * t307) * qJD(3) + t343;
t248 = t346 * qJD(2) + t316 * t305 + t340;
t247 = t303 + t316 * t307 + (-t347 + t348) * qJD(2);
t246 = -qJD(5) * t308 + (t347 * t305 + t346 * t307) * qJD(3) + t343;
t1 = m(2) * t315 / 0.2e1 + m(3) * (t315 + (t293 ^ 2 + t297 ^ 2) * t314) / 0.2e1 + t314 * Icges(3,3) / 0.2e1 + m(4) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(5) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(6) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + (((-t313 * t278 - t312 * t280 + t373 * t306 + t375 * t308) * t307 + (t313 * t279 + t312 * t281 + t372 * t306 - t374 * t308) * t305) * qJD(3) + (t313 * t300 + t312 * t301 + t370 * t306 - t371 * t308) * qJD(2)) * qJD(2) / 0.2e1 + ((t365 * t305 ^ 2 + (t361 * t307 + (t362 - t366) * t305) * t307) * qJD(3) + (t364 * t305 + t363 * t307) * qJD(2)) * t345 / 0.2e1 - ((t366 * t307 ^ 2 + (t362 * t305 + (t361 - t365) * t307) * t305) * qJD(3) + (t363 * t305 - t364 * t307) * qJD(2)) * t344 / 0.2e1;
T = t1;
