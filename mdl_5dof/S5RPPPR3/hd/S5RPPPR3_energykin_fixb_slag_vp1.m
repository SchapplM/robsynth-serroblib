% Calculate kinetic energy for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:51
% EndTime: 2019-12-31 17:43:51
% DurationCPUTime: 0.33s
% Computational Cost: add. (427->110), mult. (561->172), div. (0->0), fcn. (552->8), ass. (0->57)
t290 = cos(pkin(8));
t311 = pkin(4) * t290;
t292 = sin(qJ(1));
t310 = t292 * pkin(1);
t294 = cos(qJ(1));
t285 = qJD(1) * t294 * pkin(1);
t288 = qJ(1) + pkin(7);
t286 = sin(t288);
t287 = cos(t288);
t309 = qJD(1) * (t287 * pkin(2) + t286 * qJ(3)) + t285;
t284 = qJD(3) * t286;
t289 = sin(pkin(8));
t307 = qJD(4) * t289;
t308 = t287 * t307 + t284;
t306 = qJD(5) * t286;
t305 = qJD(5) * t287;
t304 = t286 * qJD(1);
t303 = -t286 * pkin(2) + t287 * qJ(3) - t310;
t283 = -qJD(4) * t290 + qJD(2);
t298 = pkin(3) * t290 + qJ(4) * t289;
t302 = qJD(1) * t298 * t287 + t286 * t307 + t309;
t301 = -t298 * t286 + t303;
t300 = rSges(4,1) * t290 - rSges(4,2) * t289;
t299 = rSges(5,1) * t290 + rSges(5,3) * t289;
t291 = sin(qJ(5));
t293 = cos(qJ(5));
t278 = t289 * t293 - t290 * t291;
t297 = t289 * t291 + t290 * t293;
t295 = qJD(2) ^ 2;
t282 = t294 * rSges(2,1) - t292 * rSges(2,2);
t281 = t292 * rSges(2,1) + t294 * rSges(2,2);
t272 = t285 + qJD(1) * (t287 * rSges(3,1) - t286 * rSges(3,2));
t271 = (-t286 * rSges(3,1) - t287 * rSges(3,2) - t310) * qJD(1);
t270 = t297 * t287;
t269 = t278 * t287;
t268 = t297 * t286;
t267 = t278 * t286;
t266 = t278 * rSges(6,1) - rSges(6,2) * t297;
t265 = Icges(6,1) * t278 - Icges(6,4) * t297;
t264 = Icges(6,4) * t278 - Icges(6,2) * t297;
t263 = Icges(6,5) * t278 - Icges(6,6) * t297;
t262 = t270 * rSges(6,1) + t269 * rSges(6,2) - t286 * rSges(6,3);
t261 = t268 * rSges(6,1) + t267 * rSges(6,2) + t287 * rSges(6,3);
t260 = Icges(6,1) * t270 + Icges(6,4) * t269 - Icges(6,5) * t286;
t259 = Icges(6,1) * t268 + Icges(6,4) * t267 + Icges(6,5) * t287;
t258 = Icges(6,4) * t270 + Icges(6,2) * t269 - Icges(6,6) * t286;
t257 = Icges(6,4) * t268 + Icges(6,2) * t267 + Icges(6,6) * t287;
t256 = Icges(6,5) * t270 + Icges(6,6) * t269 - Icges(6,3) * t286;
t255 = Icges(6,5) * t268 + Icges(6,6) * t267 + Icges(6,3) * t287;
t254 = rSges(4,3) * t304 + (qJD(1) * t300 - qJD(3)) * t287 + t309;
t253 = t284 + (t287 * rSges(4,3) - t300 * t286 + t303) * qJD(1);
t252 = rSges(5,2) * t304 + (qJD(1) * t299 - qJD(3)) * t287 + t302;
t251 = (t287 * rSges(5,2) - t299 * t286 + t301) * qJD(1) + t308;
t250 = (-t261 * t286 - t262 * t287) * qJD(5) + t283;
t249 = t266 * t306 - qJD(3) * t287 + (-t286 * pkin(6) + t287 * t311 + t262) * qJD(1) + t302;
t248 = t266 * t305 + (-t287 * pkin(6) - t286 * t311 - t261 + t301) * qJD(1) + t308;
t1 = m(3) * (t271 ^ 2 + t272 ^ 2 + t295) / 0.2e1 + m(4) * (t253 ^ 2 + t254 ^ 2 + t295) / 0.2e1 + m(5) * (t251 ^ 2 + t252 ^ 2 + t283 ^ 2) / 0.2e1 + m(6) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 - ((-t286 * t263 + t269 * t264 + t270 * t265) * qJD(1) + (-(-t286 * t256 + t269 * t258 + t270 * t260) * t286 + (-t286 * t255 + t269 * t257 + t270 * t259) * t287) * qJD(5)) * t306 / 0.2e1 + ((t287 * t263 + t267 * t264 + t268 * t265) * qJD(1) + (-(t287 * t256 + t267 * t258 + t268 * t260) * t286 + (t287 * t255 + t267 * t257 + t268 * t259) * t287) * qJD(5)) * t305 / 0.2e1 + qJD(1) * ((-t264 * t297 + t278 * t265) * qJD(1) + (-(-t258 * t297 + t278 * t260) * t286 + (-t257 * t297 + t278 * t259) * t287) * qJD(5)) / 0.2e1 + (m(2) * (t281 ^ 2 + t282 ^ 2) + Icges(2,3) + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t290 ^ 2 + ((Icges(4,1) + Icges(5,1)) * t289 + 0.2e1 * (Icges(4,4) - Icges(5,5)) * t290) * t289) * qJD(1) ^ 2 / 0.2e1;
T = t1;
