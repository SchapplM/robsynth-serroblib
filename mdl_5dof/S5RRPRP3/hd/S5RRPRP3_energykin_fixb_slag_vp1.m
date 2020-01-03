% Calculate kinetic energy for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:54
% EndTime: 2019-12-31 19:50:55
% DurationCPUTime: 0.90s
% Computational Cost: add. (683->104), mult. (540->165), div. (0->0), fcn. (435->8), ass. (0->69)
t370 = Icges(5,4) - Icges(6,5);
t369 = Icges(5,1) + Icges(6,1);
t368 = Icges(5,2) + Icges(6,3);
t296 = pkin(8) + qJ(4);
t292 = cos(t296);
t367 = t370 * t292;
t291 = sin(t296);
t366 = t370 * t291;
t365 = Icges(6,4) + Icges(5,5);
t364 = Icges(5,6) - Icges(6,6);
t363 = t368 * t291 - t367;
t362 = t369 * t292 - t366;
t361 = rSges(6,1) + pkin(4);
t360 = rSges(6,3) + qJ(5);
t359 = Icges(6,2) + Icges(5,3);
t298 = qJ(1) + qJ(2);
t293 = sin(t298);
t294 = cos(t298);
t358 = t363 * t293 + t364 * t294;
t357 = -t364 * t293 + t363 * t294;
t356 = -t362 * t293 + t365 * t294;
t355 = t365 * t293 + t362 * t294;
t354 = -t368 * t292 - t366;
t353 = t369 * t291 + t367;
t352 = -t364 * t291 + t365 * t292;
t351 = t360 * t291 + t361 * t292;
t350 = t352 * t293 - t359 * t294;
t349 = t359 * t293 + t352 * t294;
t348 = t365 * t291 + t364 * t292;
t347 = t354 * t291 + t353 * t292;
t346 = t357 * t291 + t355 * t292;
t345 = -t358 * t291 + t356 * t292;
t300 = cos(pkin(8));
t339 = t300 * pkin(3);
t338 = pkin(1) * qJD(1);
t285 = t293 * pkin(2) - t294 * qJ(3);
t332 = pkin(7) * t294 - t339 * t293 - t285;
t331 = -t294 * rSges(6,2) + t351 * t293;
t330 = t293 * rSges(6,2) + t351 * t294;
t303 = cos(qJ(1));
t290 = t303 * t338;
t297 = qJD(1) + qJD(2);
t329 = t297 * (t294 * pkin(2) + t293 * qJ(3)) + t290;
t328 = qJD(4) * t293;
t327 = qJD(4) * t294;
t302 = sin(qJ(1));
t326 = t302 * t338;
t323 = qJD(3) * t293 - t326;
t299 = sin(pkin(8));
t322 = rSges(4,1) * t300 - rSges(4,2) * t299;
t321 = rSges(5,1) * t292 - rSges(5,2) * t291;
t306 = -qJD(3) * t294 + t297 * (pkin(7) * t293 + t339 * t294) + t329;
t305 = t360 * qJD(4) * t292 + (-t361 * qJD(4) + qJD(5)) * t291;
t287 = t303 * rSges(2,1) - t302 * rSges(2,2);
t286 = t302 * rSges(2,1) + t303 * rSges(2,2);
t284 = t291 * rSges(5,1) + t292 * rSges(5,2);
t272 = t290 + t297 * (t294 * rSges(3,1) - t293 * rSges(3,2));
t271 = -t326 - t297 * (t293 * rSges(3,1) + t294 * rSges(3,2));
t270 = t293 * rSges(5,3) + t321 * t294;
t268 = -t294 * rSges(5,3) + t321 * t293;
t252 = t297 * t293 * rSges(4,3) + (t297 * t322 - qJD(3)) * t294 + t329;
t251 = (t294 * rSges(4,3) - t322 * t293 - t285) * t297 + t323;
t250 = (t268 * t293 + t270 * t294) * qJD(4);
t249 = t297 * t270 - t284 * t328 + t306;
t248 = -t284 * t327 + (-t268 + t332) * t297 + t323;
t247 = -qJD(5) * t292 + (t331 * t293 + t330 * t294) * qJD(4);
t246 = t305 * t293 + t330 * t297 + t306;
t245 = t305 * t294 + (-t331 + t332) * t297 + t323;
t1 = m(3) * (t271 ^ 2 + t272 ^ 2) / 0.2e1 + m(4) * (t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(5) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(6) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t300 ^ 2 + (Icges(4,1) * t299 + 0.2e1 * Icges(4,4) * t300) * t299) * t297 ^ 2 / 0.2e1 + ((t353 * t291 - t354 * t292) * t297 + ((t356 * t291 + t358 * t292) * t294 + (t355 * t291 - t357 * t292) * t293) * qJD(4)) * t297 / 0.2e1 + (m(2) * (t286 ^ 2 + t287 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t348 * t293 + t347 * t294) * t297 + (t349 * t293 ^ 2 + (t345 * t294 + (t346 - t350) * t293) * t294) * qJD(4)) * t328 / 0.2e1 - ((t347 * t293 - t348 * t294) * t297 + (t350 * t294 ^ 2 + (t346 * t293 + (t345 - t349) * t294) * t293) * qJD(4)) * t327 / 0.2e1;
T = t1;
