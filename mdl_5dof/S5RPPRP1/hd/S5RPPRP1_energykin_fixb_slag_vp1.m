% Calculate kinetic energy for
% S5RPPRP1
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:07
% EndTime: 2020-01-03 11:25:08
% DurationCPUTime: 0.96s
% Computational Cost: add. (790->130), mult. (1002->203), div. (0->0), fcn. (1007->8), ass. (0->68)
t346 = Icges(5,1) + Icges(6,1);
t345 = Icges(5,4) + Icges(6,4);
t344 = -Icges(6,5) - Icges(5,5);
t343 = Icges(5,2) + Icges(6,2);
t342 = -Icges(6,6) - Icges(5,6);
t341 = Icges(5,3) + Icges(6,3);
t296 = sin(pkin(8));
t295 = qJ(1) + pkin(7);
t293 = sin(t295);
t294 = cos(t295);
t301 = cos(qJ(4));
t297 = cos(pkin(8));
t299 = sin(qJ(4));
t321 = t297 * t299;
t280 = -t293 * t321 - t294 * t301;
t320 = t297 * t301;
t322 = t294 * t299;
t281 = t293 * t320 - t322;
t282 = -t293 * t301 + t294 * t321;
t324 = t293 * t299;
t283 = -t294 * t320 - t324;
t323 = t294 * t296;
t325 = t293 * t296;
t332 = (t342 * t282 + t344 * t283 + t341 * t323) * t294 + (-t342 * t280 - t344 * t281 + t341 * t325) * t293;
t340 = t332 * t296;
t339 = t343 * t280 + t345 * t281 - t342 * t325;
t338 = t343 * t282 + t345 * t283 + t342 * t323;
t337 = t345 * t280 + t346 * t281 - t344 * t325;
t336 = -t345 * t282 - t346 * t283 - t344 * t323;
t335 = t341 * t297 + (-t342 * t299 + t344 * t301) * t296;
t334 = t342 * t297 + (-t343 * t299 + t345 * t301) * t296;
t333 = t344 * t297 + (-t345 * t299 + t346 * t301) * t296;
t328 = t301 * pkin(4);
t331 = qJ(5) * t296 + t297 * t328;
t326 = pkin(1) * qJD(1);
t319 = t281 * rSges(6,1) + t280 * rSges(6,2) + rSges(6,3) * t325 - pkin(4) * t322 + t331 * t293;
t318 = t283 * rSges(6,1) + t282 * rSges(6,2) - rSges(6,3) * t323 - pkin(4) * t324 - t331 * t294;
t300 = sin(qJ(1));
t291 = t300 * t326;
t317 = qJD(1) * (t293 * pkin(2) - t294 * qJ(3)) + t291;
t316 = qJD(4) * t296;
t315 = qJD(5) * t296;
t311 = pkin(3) * t297 + pkin(6) * t296;
t314 = qJD(1) * t311 * t293 + t317;
t310 = rSges(4,1) * t297 - rSges(4,2) * t296;
t309 = -(-t297 * rSges(5,3) + (rSges(5,1) * t301 - rSges(5,2) * t299) * t296) * t316 - qJD(3);
t286 = -t294 * pkin(2) - t293 * qJ(3);
t302 = cos(qJ(1));
t292 = t302 * t326;
t306 = t292 + (t311 * t294 - t286) * qJD(1);
t305 = -qJD(3) + ((qJ(5) + rSges(6,3)) * t297 + (-rSges(6,1) * t301 + rSges(6,2) * t299 - t328) * t296) * t316;
t303 = qJD(2) ^ 2;
t289 = -qJD(4) * t297 + qJD(1);
t288 = -t302 * rSges(2,1) + t300 * rSges(2,2);
t287 = t300 * rSges(2,1) + t302 * rSges(2,2);
t270 = t292 - qJD(1) * (-t294 * rSges(3,1) + t293 * rSges(3,2));
t269 = t291 + qJD(1) * (t293 * rSges(3,1) + t294 * rSges(3,2));
t267 = t283 * rSges(5,1) + t282 * rSges(5,2) - rSges(5,3) * t323;
t265 = t281 * rSges(5,1) + t280 * rSges(5,2) + rSges(5,3) * t325;
t249 = -qJD(3) * t294 + t292 + (t293 * rSges(4,3) + t294 * t310 - t286) * qJD(1);
t248 = -qJD(1) * t294 * rSges(4,3) + (qJD(1) * t310 - qJD(3)) * t293 + t317;
t247 = qJD(2) + (t265 * t294 + t267 * t293) * t316;
t246 = -t289 * t267 + t294 * t309 + t306;
t245 = t289 * t265 + t293 * t309 + t314;
t244 = -t289 * t318 + t293 * t315 + t294 * t305 + t306;
t243 = t289 * t319 + t293 * t305 - t294 * t315 + t314;
t242 = -qJD(5) * t297 + qJD(2) + (t293 * t318 + t294 * t319) * t316;
t1 = m(3) * (t269 ^ 2 + t270 ^ 2 + t303) / 0.2e1 + m(4) * (t248 ^ 2 + t249 ^ 2 + t303) / 0.2e1 + m(5) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + m(6) * (t242 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + ((-t332 * t297 + ((t338 * t299 + t336 * t301) * t294 + (-t339 * t299 + t337 * t301) * t293) * t296) * t316 + (t335 * t297 + (-t334 * t299 + t333 * t301) * t296) * t289) * t289 / 0.2e1 + (((-t338 * t280 + t336 * t281) * t294 + (t339 * t280 + t337 * t281 + t340) * t293) * t316 + (t334 * t280 + t333 * t281 - t335 * t325) * t289) * t293 * t316 / 0.2e1 - (((-t338 * t282 + t336 * t283 - t340) * t294 + (t339 * t282 + t337 * t283) * t293) * t316 + (t334 * t282 + t333 * t283 + t335 * t323) * t289) * t294 * t316 / 0.2e1 + (m(2) * (t287 ^ 2 + t288 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t297 ^ 2 + (Icges(4,1) * t296 + 0.2e1 * Icges(4,4) * t297) * t296) * qJD(1) ^ 2 / 0.2e1;
T = t1;
