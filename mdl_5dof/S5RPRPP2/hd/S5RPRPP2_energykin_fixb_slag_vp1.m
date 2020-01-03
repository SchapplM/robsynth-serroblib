% Calculate kinetic energy for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:39
% EndTime: 2019-12-31 18:10:40
% DurationCPUTime: 1.14s
% Computational Cost: add. (643->111), mult. (741->169), div. (0->0), fcn. (610->6), ass. (0->73)
t386 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t385 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t384 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t310 = cos(qJ(3));
t383 = t386 * t310;
t308 = sin(qJ(3));
t382 = t386 * t308;
t381 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t380 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t379 = t384 * t308 - t383;
t378 = t385 * t310 - t382;
t377 = rSges(6,1) + pkin(4);
t376 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t307 = qJ(1) + pkin(7);
t305 = sin(t307);
t306 = cos(t307);
t375 = t379 * t305 + t380 * t306;
t374 = -t380 * t305 + t379 * t306;
t373 = -t378 * t305 + t381 * t306;
t372 = t381 * t305 + t378 * t306;
t371 = -t384 * t310 - t382;
t370 = t385 * t308 + t383;
t369 = -t380 * t308 + t381 * t310;
t368 = rSges(6,3) + qJ(5);
t367 = rSges(6,2) * t308 + t377 * t310;
t366 = t369 * t305 - t376 * t306;
t365 = t376 * t305 + t369 * t306;
t364 = -t381 * t308 - t380 * t310;
t363 = t371 * t308 + t370 * t310;
t362 = t374 * t308 + t372 * t310;
t361 = -t375 * t308 + t373 * t310;
t309 = sin(qJ(1));
t357 = pkin(1) * t309;
t348 = t367 * t305 + t368 * t306;
t347 = -t368 * t305 + t367 * t306;
t311 = cos(qJ(1));
t304 = qJD(1) * t311 * pkin(1);
t346 = qJD(1) * (pkin(2) * t306 + pkin(6) * t305) + t304;
t345 = qJD(3) * t305;
t344 = qJD(3) * t306;
t343 = qJD(4) * t308;
t340 = -pkin(2) * t305 + pkin(6) * t306 - t357;
t298 = pkin(3) * t308 - qJ(4) * t310;
t339 = qJD(3) * (-rSges(5,1) * t308 + rSges(5,3) * t310 - t298);
t333 = pkin(3) * t310 + qJ(4) * t308;
t282 = t333 * t306;
t338 = qJD(1) * t282 + t305 * t343 + t346;
t281 = t333 * t305;
t337 = -t281 + t340;
t336 = rSges(4,1) * t310 - rSges(4,2) * t308;
t335 = rSges(5,1) * t310 + rSges(5,3) * t308;
t314 = -qJD(4) * t310 + t281 * t345 + t282 * t344 + qJD(2);
t313 = qJD(3) * (rSges(6,2) * t310 - t377 * t308 - t298);
t303 = rSges(2,1) * t311 - rSges(2,2) * t309;
t302 = rSges(2,1) * t309 + rSges(2,2) * t311;
t301 = rSges(4,1) * t308 + rSges(4,2) * t310;
t297 = t306 * t343;
t279 = t304 + qJD(1) * (rSges(3,1) * t306 - rSges(3,2) * t305);
t278 = (-rSges(3,1) * t305 - rSges(3,2) * t306 - t357) * qJD(1);
t277 = rSges(4,3) * t305 + t336 * t306;
t276 = rSges(5,2) * t305 + t335 * t306;
t274 = -rSges(4,3) * t306 + t336 * t305;
t273 = -rSges(5,2) * t306 + t335 * t305;
t251 = qJD(1) * t277 - t301 * t345 + t346;
t250 = -t301 * t344 + (-t274 + t340) * qJD(1);
t249 = qJD(2) + (t274 * t305 + t277 * t306) * qJD(3);
t248 = qJD(1) * t276 + t305 * t339 + t338;
t247 = t297 + t306 * t339 + (-t273 + t337) * qJD(1);
t246 = (t273 * t305 + t276 * t306) * qJD(3) + t314;
t245 = t347 * qJD(1) + qJD(5) * t306 + t305 * t313 + t338;
t244 = -qJD(5) * t305 + t297 + t306 * t313 + (t337 - t348) * qJD(1);
t243 = (t348 * t305 + t347 * t306) * qJD(3) + t314;
t1 = m(3) * (qJD(2) ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(4) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(5) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(6) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + (m(2) * (t302 ^ 2 + t303 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + (((t373 * t308 + t375 * t310) * t306 + (t372 * t308 - t374 * t310) * t305) * qJD(3) + (t370 * t308 - t371 * t310) * qJD(1)) * qJD(1) / 0.2e1 + ((t365 * t305 ^ 2 + (t361 * t306 + (t362 - t366) * t305) * t306) * qJD(3) + (-t364 * t305 + t363 * t306) * qJD(1)) * t345 / 0.2e1 - ((t366 * t306 ^ 2 + (t362 * t305 + (t361 - t365) * t306) * t305) * qJD(3) + (t363 * t305 + t364 * t306) * qJD(1)) * t344 / 0.2e1;
T = t1;
