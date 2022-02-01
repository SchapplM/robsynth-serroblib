% Calculate kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:11
% EndTime: 2022-01-20 09:12:12
% DurationCPUTime: 0.61s
% Computational Cost: add. (609->129), mult. (605->212), div. (0->0), fcn. (578->10), ass. (0->64)
t292 = cos(pkin(9));
t293 = cos(pkin(8));
t313 = t292 * t293;
t324 = pkin(4) * t313;
t295 = sin(qJ(1));
t321 = t295 * pkin(1);
t289 = qJ(1) + pkin(7);
t285 = sin(t289);
t291 = sin(pkin(8));
t319 = t285 * t291;
t318 = t285 * t293;
t287 = cos(t289);
t290 = sin(pkin(9));
t317 = t287 * t290;
t316 = t287 * t291;
t315 = t287 * t293;
t314 = t290 * t293;
t296 = cos(qJ(1));
t283 = qJD(1) * t296 * pkin(1);
t311 = qJD(1) * (t287 * pkin(2) + t285 * qJ(3)) + t283;
t281 = qJD(3) * t285;
t308 = qJD(4) * t291;
t310 = t287 * t308 + t281;
t309 = qJD(1) * t285;
t307 = qJD(5) * t291;
t305 = -t285 * pkin(2) + t287 * qJ(3) - t321;
t279 = -qJD(4) * t293 + qJD(2);
t301 = pkin(3) * t293 + qJ(4) * t291;
t304 = qJD(1) * t301 * t287 + t285 * t308 + t311;
t303 = -t301 * t285 + t305;
t302 = rSges(4,1) * t293 - rSges(4,2) * t291;
t288 = pkin(9) + qJ(5);
t284 = sin(t288);
t286 = cos(t288);
t265 = -t284 * t318 - t287 * t286;
t266 = -t287 * t284 + t286 * t318;
t267 = -t284 * t315 + t285 * t286;
t268 = t285 * t284 + t286 * t315;
t300 = (Icges(6,5) * t266 + Icges(6,6) * t265 + Icges(6,3) * t319) * t285 + (Icges(6,5) * t268 + Icges(6,6) * t267 + Icges(6,3) * t316) * t287;
t299 = t300 * t291;
t297 = qJD(2) ^ 2;
t280 = -qJD(5) * t293 + qJD(1);
t278 = t296 * rSges(2,1) - t295 * rSges(2,2);
t277 = t295 * rSges(2,1) + t296 * rSges(2,2);
t270 = t283 + qJD(1) * (t287 * rSges(3,1) - t285 * rSges(3,2));
t269 = (-t285 * rSges(3,1) - t287 * rSges(3,2) - t321) * qJD(1);
t264 = -t293 * rSges(6,3) + (rSges(6,1) * t286 - rSges(6,2) * t284) * t291;
t263 = -Icges(6,5) * t293 + (Icges(6,1) * t286 - Icges(6,4) * t284) * t291;
t262 = -Icges(6,6) * t293 + (Icges(6,4) * t286 - Icges(6,2) * t284) * t291;
t261 = -Icges(6,3) * t293 + (Icges(6,5) * t286 - Icges(6,6) * t284) * t291;
t260 = rSges(4,3) * t309 + (qJD(1) * t302 - qJD(3)) * t287 + t311;
t259 = t281 + (t287 * rSges(4,3) - t302 * t285 + t305) * qJD(1);
t258 = t268 * rSges(6,1) + t267 * rSges(6,2) + rSges(6,3) * t316;
t257 = t266 * rSges(6,1) + t265 * rSges(6,2) + rSges(6,3) * t319;
t256 = Icges(6,1) * t268 + Icges(6,4) * t267 + Icges(6,5) * t316;
t255 = Icges(6,1) * t266 + Icges(6,4) * t265 + Icges(6,5) * t319;
t254 = Icges(6,4) * t268 + Icges(6,2) * t267 + Icges(6,6) * t316;
t253 = Icges(6,4) * t266 + Icges(6,2) * t265 + Icges(6,6) * t319;
t250 = -qJD(3) * t287 + qJD(1) * ((t285 * t290 + t287 * t313) * rSges(5,1) + (t285 * t292 - t287 * t314) * rSges(5,2) + rSges(5,3) * t316) + t304;
t249 = (-(t285 * t313 - t317) * rSges(5,1) - (-t285 * t314 - t287 * t292) * rSges(5,2) - rSges(5,3) * t319 + t303) * qJD(1) + t310;
t248 = (t257 * t287 - t258 * t285) * t307 + t279;
t247 = t290 * pkin(4) * t309 + t280 * t258 + (-qJD(3) + qJD(1) * t324 + (pkin(6) * qJD(1) - qJD(5) * t264) * t291) * t287 + t304;
t246 = t285 * t264 * t307 - t280 * t257 + (pkin(4) * t317 + (-pkin(6) * t291 - t324) * t285 + t303) * qJD(1) + t310;
t1 = m(3) * (t269 ^ 2 + t270 ^ 2 + t297) / 0.2e1 + m(4) * (t259 ^ 2 + t260 ^ 2 + t297) / 0.2e1 + m(5) * (t249 ^ 2 + t250 ^ 2 + t279 ^ 2) / 0.2e1 + m(6) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t280 * ((-t293 * t261 + (-t262 * t284 + t263 * t286) * t291) * t280 + (((-t254 * t284 + t256 * t286) * t287 + (-t253 * t284 + t255 * t286) * t285) * t291 - t300 * t293) * t307) / 0.2e1 + (t287 * ((t261 * t316 + t267 * t262 + t268 * t263) * t280 + ((t267 * t253 + t268 * t255) * t285 + (t267 * t254 + t268 * t256 + t299) * t287) * t307) + t285 * ((t261 * t319 + t265 * t262 + t266 * t263) * t280 + ((t265 * t254 + t266 * t256) * t287 + (t265 * t253 + t266 * t255 + t299) * t285) * t307)) * t307 / 0.2e1 + (m(2) * (t277 ^ 2 + t278 ^ 2) + Icges(2,3) + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t293 ^ 2 + ((Icges(4,1) + Icges(5,1) * t292 ^ 2 + (-0.2e1 * Icges(5,4) * t292 + Icges(5,2) * t290) * t290) * t291 + 0.2e1 * (-Icges(5,5) * t292 + Icges(5,6) * t290 + Icges(4,4)) * t293) * t291) * qJD(1) ^ 2 / 0.2e1;
T = t1;
