% Calculate kinetic energy for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:56
% EndTime: 2019-12-31 18:06:56
% DurationCPUTime: 0.73s
% Computational Cost: add. (570->160), mult. (765->263), div. (0->0), fcn. (716->8), ass. (0->85)
t305 = sin(qJ(1));
t307 = cos(qJ(1));
t336 = qJD(1) * t307 * qJ(3) + qJD(3) * t305;
t301 = sin(pkin(8));
t334 = pkin(3) * t301;
t300 = pkin(8) + qJ(4);
t295 = sin(t300);
t333 = Icges(5,4) * t295;
t296 = cos(t300);
t332 = Icges(5,4) * t296;
t331 = t296 * t305;
t330 = t296 * t307;
t304 = sin(qJ(5));
t329 = t304 * t305;
t328 = t304 * t307;
t306 = cos(qJ(5));
t327 = t305 * t306;
t326 = t306 * t307;
t299 = qJD(2) * t305;
t324 = qJD(3) * t307 + t299;
t323 = qJD(4) * t305;
t322 = qJD(4) * t307;
t321 = qJD(5) * t296;
t289 = qJD(1) * (pkin(1) * t307 + qJ(2) * t305);
t320 = -qJD(2) * t307 + t289;
t319 = qJD(1) * (pkin(6) * t307 + t305 * t334) + t289 + t336;
t318 = pkin(4) * t295 - pkin(7) * t296;
t290 = pkin(1) * t305 - qJ(2) * t307;
t317 = t307 * t334 - t290 + (-pkin(6) - qJ(3)) * t305;
t302 = cos(pkin(8));
t316 = rSges(4,1) * t301 + rSges(4,2) * t302;
t315 = rSges(5,1) * t295 + rSges(5,2) * t296;
t314 = Icges(5,1) * t295 + t332;
t313 = Icges(5,2) * t296 + t333;
t312 = Icges(5,5) * t295 + Icges(5,6) * t296;
t268 = Icges(5,6) * t307 + t305 * t313;
t270 = Icges(5,5) * t307 + t305 * t314;
t311 = -t268 * t296 - t270 * t295;
t269 = Icges(5,6) * t305 - t307 * t313;
t271 = Icges(5,5) * t305 - t307 * t314;
t310 = t269 * t296 + t271 * t295;
t285 = -Icges(5,2) * t295 + t332;
t286 = Icges(5,1) * t296 - t333;
t309 = t285 * t296 + t286 * t295;
t293 = qJD(5) * t295 + qJD(1);
t292 = rSges(2,1) * t307 - rSges(2,2) * t305;
t291 = rSges(2,1) * t305 + rSges(2,2) * t307;
t288 = pkin(4) * t296 + pkin(7) * t295;
t287 = rSges(5,1) * t296 - rSges(5,2) * t295;
t284 = Icges(5,5) * t296 - Icges(5,6) * t295;
t283 = -t305 * t321 + t322;
t282 = t307 * t321 + t323;
t281 = -t295 * t326 + t329;
t280 = t295 * t328 + t327;
t279 = t295 * t327 + t328;
t278 = -t295 * t329 + t326;
t276 = t318 * t307;
t275 = t318 * t305;
t273 = rSges(5,3) * t305 - t307 * t315;
t272 = rSges(5,3) * t307 + t305 * t315;
t267 = Icges(5,3) * t305 - t307 * t312;
t266 = Icges(5,3) * t307 + t305 * t312;
t265 = rSges(6,3) * t295 + (rSges(6,1) * t306 - rSges(6,2) * t304) * t296;
t264 = Icges(6,5) * t295 + (Icges(6,1) * t306 - Icges(6,4) * t304) * t296;
t263 = Icges(6,6) * t295 + (Icges(6,4) * t306 - Icges(6,2) * t304) * t296;
t262 = Icges(6,3) * t295 + (Icges(6,5) * t306 - Icges(6,6) * t304) * t296;
t261 = qJD(1) * (-rSges(3,2) * t307 + rSges(3,3) * t305) + t320;
t260 = t299 + (rSges(3,2) * t305 + rSges(3,3) * t307 - t290) * qJD(1);
t259 = rSges(6,1) * t281 + rSges(6,2) * t280 + rSges(6,3) * t330;
t258 = rSges(6,1) * t279 + rSges(6,2) * t278 - rSges(6,3) * t331;
t257 = Icges(6,1) * t281 + Icges(6,4) * t280 + Icges(6,5) * t330;
t256 = Icges(6,1) * t279 + Icges(6,4) * t278 - Icges(6,5) * t331;
t255 = Icges(6,4) * t281 + Icges(6,2) * t280 + Icges(6,6) * t330;
t254 = Icges(6,4) * t279 + Icges(6,2) * t278 - Icges(6,6) * t331;
t253 = Icges(6,5) * t281 + Icges(6,6) * t280 + Icges(6,3) * t330;
t252 = Icges(6,5) * t279 + Icges(6,6) * t278 - Icges(6,3) * t331;
t251 = qJD(1) * (rSges(4,3) * t307 + t305 * t316) + t320 + t336;
t250 = (-t290 + t316 * t307 + (-rSges(4,3) - qJ(3)) * t305) * qJD(1) + t324;
t249 = (-t272 * t305 + t273 * t307) * qJD(4);
t248 = qJD(1) * t272 + (-qJD(4) * t287 - qJD(2)) * t307 + t319;
t247 = t287 * t323 + (-t273 + t317) * qJD(1) + t324;
t246 = -t258 * t282 + t259 * t283 + (-t275 * t305 - t276 * t307) * qJD(4);
t245 = qJD(1) * t275 + t258 * t293 - t265 * t283 + (-qJD(4) * t288 - qJD(2)) * t307 + t319;
t244 = t288 * t323 - t259 * t293 + t265 * t282 + (t276 + t317) * qJD(1) + t324;
t1 = m(3) * (t260 ^ 2 + t261 ^ 2) / 0.2e1 + m(4) * (t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(5) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + ((t307 * t284 + t305 * t309) * qJD(1) + (t307 ^ 2 * t266 + (t310 * t305 + (t267 - t311) * t307) * t305) * qJD(4)) * t322 / 0.2e1 + ((t305 * t284 - t307 * t309) * qJD(1) + (t305 ^ 2 * t267 + (t311 * t307 + (t266 - t310) * t305) * t307) * qJD(4)) * t323 / 0.2e1 + qJD(1) * ((-t295 * t285 + t296 * t286) * qJD(1) + ((-t268 * t295 + t270 * t296) * t307 + (-t269 * t295 + t271 * t296) * t305) * qJD(4)) / 0.2e1 + m(6) * (t244 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + t283 * ((-t252 * t331 + t278 * t254 + t279 * t256) * t283 + (-t253 * t331 + t255 * t278 + t257 * t279) * t282 + (-t262 * t331 + t263 * t278 + t264 * t279) * t293) / 0.2e1 + t282 * ((t252 * t330 + t254 * t280 + t281 * t256) * t283 + (t253 * t330 + t280 * t255 + t281 * t257) * t282 + (t262 * t330 + t263 * t280 + t264 * t281) * t293) / 0.2e1 + t293 * ((t252 * t283 + t253 * t282 + t262 * t293) * t295 + ((-t254 * t304 + t256 * t306) * t283 + (-t255 * t304 + t257 * t306) * t282 + (-t263 * t304 + t264 * t306) * t293) * t296) / 0.2e1 + (m(2) * (t291 ^ 2 + t292 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1) * t302 ^ 2 + (-0.2e1 * Icges(4,4) * t302 + Icges(4,2) * t301) * t301) * qJD(1) ^ 2 / 0.2e1;
T = t1;
