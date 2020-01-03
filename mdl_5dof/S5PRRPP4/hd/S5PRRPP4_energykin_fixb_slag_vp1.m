% Calculate kinetic energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:45
% EndTime: 2019-12-31 17:40:46
% DurationCPUTime: 1.05s
% Computational Cost: add. (631->103), mult. (715->161), div. (0->0), fcn. (598->4), ass. (0->69)
t379 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t378 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t377 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t306 = sin(qJ(3));
t376 = t379 * t306;
t307 = cos(qJ(3));
t375 = t379 * t307;
t374 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t373 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t372 = t377 * t306 - t375;
t371 = t378 * t307 - t376;
t370 = rSges(6,1) + pkin(4);
t369 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t305 = pkin(7) + qJ(2);
t303 = sin(t305);
t304 = cos(t305);
t368 = t372 * t303 + t373 * t304;
t367 = -t373 * t303 + t372 * t304;
t366 = -t371 * t303 + t374 * t304;
t365 = t374 * t303 + t371 * t304;
t364 = -t377 * t307 - t376;
t363 = t378 * t306 + t375;
t362 = -t373 * t306 + t374 * t307;
t361 = rSges(6,3) + qJ(5);
t360 = rSges(6,2) * t306 + t370 * t307;
t359 = t362 * t303 - t369 * t304;
t358 = t369 * t303 + t362 * t304;
t357 = -t374 * t306 - t373 * t307;
t356 = t364 * t306 + t363 * t307;
t355 = t367 * t306 + t365 * t307;
t354 = -t368 * t306 + t366 * t307;
t343 = t360 * t303 + t361 * t304;
t342 = -t361 * t303 + t360 * t304;
t330 = pkin(3) * t307 + qJ(4) * t306;
t280 = t330 * t303;
t287 = t303 * pkin(2) - t304 * pkin(6);
t341 = -t280 - t287;
t340 = qJD(3) * t303;
t339 = qJD(3) * t304;
t338 = qJD(4) * t306;
t281 = t330 * t304;
t284 = qJD(2) * (t304 * pkin(2) + t303 * pkin(6));
t337 = qJD(2) * t281 + t303 * t338 + t284;
t299 = t306 * pkin(3) - t307 * qJ(4);
t334 = qJD(3) * (-t306 * rSges(5,1) + t307 * rSges(5,3) - t299);
t333 = rSges(4,1) * t307 - rSges(4,2) * t306;
t332 = rSges(5,1) * t307 + rSges(5,3) * t306;
t311 = -qJD(4) * t307 + t280 * t340 + t281 * t339 + qJD(1);
t310 = qJD(3) * (t307 * rSges(6,2) - t370 * t306 - t299);
t309 = qJD(1) ^ 2;
t308 = qJD(2) ^ 2;
t302 = t306 * rSges(4,1) + t307 * rSges(4,2);
t298 = t304 * t338;
t286 = t304 * rSges(3,1) - t303 * rSges(3,2);
t285 = t303 * rSges(3,1) + t304 * rSges(3,2);
t278 = t303 * rSges(4,3) + t333 * t304;
t277 = t303 * rSges(5,2) + t332 * t304;
t275 = -t304 * rSges(4,3) + t333 * t303;
t274 = -t304 * rSges(5,2) + t332 * t303;
t252 = qJD(2) * t278 - t302 * t340 + t284;
t251 = -t302 * t339 + (-t275 - t287) * qJD(2);
t250 = qJD(1) + (t275 * t303 + t278 * t304) * qJD(3);
t249 = qJD(2) * t277 + t303 * t334 + t337;
t248 = t298 + t304 * t334 + (-t274 + t341) * qJD(2);
t247 = (t274 * t303 + t277 * t304) * qJD(3) + t311;
t246 = t342 * qJD(2) + qJD(5) * t304 + t303 * t310 + t337;
t245 = -qJD(5) * t303 + t298 + t304 * t310 + (t341 - t343) * qJD(2);
t244 = (t343 * t303 + t342 * t304) * qJD(3) + t311;
t1 = m(2) * t309 / 0.2e1 + m(3) * (t309 + (t285 ^ 2 + t286 ^ 2) * t308) / 0.2e1 + t308 * Icges(3,3) / 0.2e1 + m(4) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(5) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + m(6) * (t244 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + (((t366 * t306 + t368 * t307) * t304 + (t365 * t306 - t367 * t307) * t303) * qJD(3) + (t363 * t306 - t364 * t307) * qJD(2)) * qJD(2) / 0.2e1 + ((t358 * t303 ^ 2 + (t354 * t304 + (t355 - t359) * t303) * t304) * qJD(3) + (-t357 * t303 + t356 * t304) * qJD(2)) * t340 / 0.2e1 - ((t359 * t304 ^ 2 + (t355 * t303 + (t354 - t358) * t304) * t303) * qJD(3) + (t356 * t303 + t357 * t304) * qJD(2)) * t339 / 0.2e1;
T = t1;
