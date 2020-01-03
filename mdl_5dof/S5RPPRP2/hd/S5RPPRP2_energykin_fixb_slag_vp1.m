% Calculate kinetic energy for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:04
% EndTime: 2019-12-31 17:49:05
% DurationCPUTime: 0.94s
% Computational Cost: add. (657->107), mult. (542->162), div. (0->0), fcn. (435->8), ass. (0->68)
t369 = Icges(5,4) - Icges(6,5);
t368 = Icges(5,1) + Icges(6,1);
t367 = Icges(5,2) + Icges(6,3);
t297 = pkin(8) + qJ(4);
t295 = cos(t297);
t366 = t369 * t295;
t293 = sin(t297);
t365 = t369 * t293;
t364 = Icges(6,4) + Icges(5,5);
t363 = Icges(5,6) - Icges(6,6);
t362 = t367 * t293 - t366;
t361 = t368 * t295 - t365;
t360 = rSges(6,1) + pkin(4);
t359 = rSges(6,3) + qJ(5);
t358 = Icges(6,2) + Icges(5,3);
t298 = qJ(1) + pkin(7);
t294 = sin(t298);
t296 = cos(t298);
t357 = t362 * t294 + t363 * t296;
t356 = -t363 * t294 + t362 * t296;
t355 = -t361 * t294 + t364 * t296;
t354 = t364 * t294 + t361 * t296;
t353 = -t367 * t295 - t365;
t352 = t368 * t293 + t366;
t351 = -t363 * t293 + t364 * t295;
t350 = t359 * t293 + t360 * t295;
t349 = t351 * t294 - t358 * t296;
t348 = t358 * t294 + t351 * t296;
t347 = t364 * t293 + t363 * t295;
t346 = t353 * t293 + t352 * t295;
t345 = t356 * t293 + t354 * t295;
t344 = -t357 * t293 + t355 * t295;
t302 = sin(qJ(1));
t340 = t302 * pkin(1);
t300 = cos(pkin(8));
t338 = pkin(3) * t300;
t332 = -rSges(6,2) * t296 + t350 * t294;
t331 = rSges(6,2) * t294 + t350 * t296;
t303 = cos(qJ(1));
t292 = qJD(1) * t303 * pkin(1);
t330 = qJD(1) * (pkin(2) * t296 + qJ(3) * t294) + t292;
t329 = qJD(4) * t294;
t328 = qJD(4) * t296;
t325 = -pkin(2) * t294 + qJ(3) * t296 - t340;
t324 = pkin(6) * t296 - t338 * t294 + t325;
t299 = sin(pkin(8));
t323 = rSges(4,1) * t300 - rSges(4,2) * t299;
t322 = rSges(5,1) * t295 - rSges(5,2) * t293;
t307 = -qJD(3) * t296 + qJD(1) * (pkin(6) * t294 + t338 * t296) + t330;
t306 = t359 * qJD(4) * t295 + (-t360 * qJD(4) + qJD(5)) * t293;
t304 = qJD(2) ^ 2;
t290 = qJD(3) * t294;
t289 = rSges(2,1) * t303 - rSges(2,2) * t302;
t288 = rSges(2,1) * t302 + rSges(2,2) * t303;
t286 = rSges(5,1) * t293 + rSges(5,2) * t295;
t276 = t292 + qJD(1) * (rSges(3,1) * t296 - rSges(3,2) * t294);
t275 = (-rSges(3,1) * t294 - rSges(3,2) * t296 - t340) * qJD(1);
t272 = rSges(5,3) * t294 + t322 * t296;
t270 = -rSges(5,3) * t296 + t322 * t294;
t254 = qJD(1) * t294 * rSges(4,3) + (qJD(1) * t323 - qJD(3)) * t296 + t330;
t253 = t290 + (t296 * rSges(4,3) - t323 * t294 + t325) * qJD(1);
t252 = qJD(2) + (t270 * t294 + t272 * t296) * qJD(4);
t251 = qJD(1) * t272 - t286 * t329 + t307;
t250 = -t286 * t328 + t290 + (-t270 + t324) * qJD(1);
t249 = -qJD(5) * t295 + qJD(2) + (t332 * t294 + t331 * t296) * qJD(4);
t248 = t331 * qJD(1) + t306 * t294 + t307;
t247 = t290 + t306 * t296 + (t324 - t332) * qJD(1);
t1 = m(3) * (t275 ^ 2 + t276 ^ 2 + t304) / 0.2e1 + m(4) * (t253 ^ 2 + t254 ^ 2 + t304) / 0.2e1 + m(5) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(6) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + (((t355 * t293 + t357 * t295) * t296 + (t354 * t293 - t356 * t295) * t294) * qJD(4) + (t352 * t293 - t353 * t295) * qJD(1)) * qJD(1) / 0.2e1 + ((t348 * t294 ^ 2 + (t344 * t296 + (t345 - t349) * t294) * t296) * qJD(4) + (t347 * t294 + t346 * t296) * qJD(1)) * t329 / 0.2e1 - ((t349 * t296 ^ 2 + (t345 * t294 + (t344 - t348) * t296) * t294) * qJD(4) + (t346 * t294 - t347 * t296) * qJD(1)) * t328 / 0.2e1 + (m(2) * (t288 ^ 2 + t289 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t300 ^ 2 + (Icges(4,1) * t299 + 0.2e1 * Icges(4,4) * t300) * t299) * qJD(1) ^ 2 / 0.2e1;
T = t1;
