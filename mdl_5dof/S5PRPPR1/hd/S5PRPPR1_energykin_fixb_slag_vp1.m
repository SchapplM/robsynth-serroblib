% Calculate kinetic energy for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:45
% EndTime: 2019-12-05 15:21:45
% DurationCPUTime: 0.52s
% Computational Cost: add. (597->119), mult. (579->201), div. (0->0), fcn. (566->8), ass. (0->58)
t290 = cos(pkin(9));
t291 = cos(pkin(8));
t308 = t290 * t291;
t286 = pkin(9) + qJ(5);
t282 = sin(t286);
t284 = cos(t286);
t289 = sin(pkin(8));
t317 = -t289 * (qJD(2) * pkin(6) - qJD(5) * (-t291 * rSges(6,3) + (rSges(6,1) * t284 - rSges(6,2) * t282) * t289)) - qJD(2) * pkin(4) * t308;
t287 = pkin(7) + qJ(2);
t283 = sin(t287);
t314 = t283 * t289;
t313 = t283 * t291;
t285 = cos(t287);
t288 = sin(pkin(9));
t312 = t285 * t288;
t311 = t285 * t289;
t310 = t285 * t291;
t309 = t288 * t291;
t273 = t283 * pkin(2) - t285 * qJ(3);
t298 = pkin(3) * t291 + qJ(4) * t289;
t307 = -t298 * t283 - t273;
t280 = qJD(3) * t283;
t304 = qJD(4) * t289;
t306 = t285 * t304 + t280;
t305 = qJD(2) * t283;
t303 = qJD(5) * t289;
t272 = qJD(2) * (t285 * pkin(2) + t283 * qJ(3));
t302 = qJD(2) * t298 * t285 + t283 * t304 + t272;
t279 = -qJD(4) * t291 + qJD(1);
t299 = rSges(4,1) * t291 - rSges(4,2) * t289;
t266 = -t282 * t313 - t285 * t284;
t267 = -t285 * t282 + t284 * t313;
t268 = -t282 * t310 + t283 * t284;
t269 = t283 * t282 + t284 * t310;
t297 = (Icges(6,5) * t267 + Icges(6,6) * t266 + Icges(6,3) * t314) * t283 + (Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t311) * t285;
t296 = t297 * t289;
t294 = qJD(1) ^ 2;
t293 = qJD(2) ^ 2;
t278 = -qJD(5) * t291 + qJD(2);
t275 = t285 * rSges(3,1) - t283 * rSges(3,2);
t274 = t283 * rSges(3,1) + t285 * rSges(3,2);
t264 = -Icges(6,5) * t291 + (Icges(6,1) * t284 - Icges(6,4) * t282) * t289;
t263 = -Icges(6,6) * t291 + (Icges(6,4) * t284 - Icges(6,2) * t282) * t289;
t262 = -Icges(6,3) * t291 + (Icges(6,5) * t284 - Icges(6,6) * t282) * t289;
t261 = rSges(4,3) * t305 + t272 + (qJD(2) * t299 - qJD(3)) * t285;
t260 = t280 + (t285 * rSges(4,3) - t283 * t299 - t273) * qJD(2);
t259 = t269 * rSges(6,1) + t268 * rSges(6,2) + rSges(6,3) * t311;
t258 = t267 * rSges(6,1) + t266 * rSges(6,2) + rSges(6,3) * t314;
t257 = Icges(6,1) * t269 + Icges(6,4) * t268 + Icges(6,5) * t311;
t256 = Icges(6,1) * t267 + Icges(6,4) * t266 + Icges(6,5) * t314;
t255 = Icges(6,4) * t269 + Icges(6,2) * t268 + Icges(6,6) * t311;
t254 = Icges(6,4) * t267 + Icges(6,2) * t266 + Icges(6,6) * t314;
t251 = -qJD(3) * t285 + qJD(2) * ((t283 * t288 + t285 * t308) * rSges(5,1) + (t283 * t290 - t285 * t309) * rSges(5,2) + rSges(5,3) * t311) + t302;
t250 = (-(t283 * t308 - t312) * rSges(5,1) - (-t283 * t309 - t285 * t290) * rSges(5,2) - rSges(5,3) * t314 + t307) * qJD(2) + t306;
t249 = (t258 * t285 - t259 * t283) * t303 + t279;
t248 = t288 * pkin(4) * t305 + t278 * t259 + (-qJD(3) - t317) * t285 + t302;
t247 = -t278 * t258 + (pkin(4) * t312 + t307) * qJD(2) + t317 * t283 + t306;
t1 = m(2) * t294 / 0.2e1 + m(3) * (t294 + (t274 ^ 2 + t275 ^ 2) * t293) / 0.2e1 + m(4) * (t260 ^ 2 + t261 ^ 2 + t294) / 0.2e1 + m(5) * (t250 ^ 2 + t251 ^ 2 + t279 ^ 2) / 0.2e1 + m(6) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t278 * ((-t291 * t262 + (-t263 * t282 + t264 * t284) * t289) * t278 + (((-t255 * t282 + t257 * t284) * t285 + (-t254 * t282 + t256 * t284) * t283) * t289 - t297 * t291) * t303) / 0.2e1 + (t285 * ((t262 * t311 + t268 * t263 + t269 * t264) * t278 + ((t268 * t254 + t269 * t256) * t283 + (t268 * t255 + t269 * t257 + t296) * t285) * t303) + t283 * ((t262 * t314 + t266 * t263 + t267 * t264) * t278 + ((t266 * t255 + t267 * t257) * t285 + (t266 * t254 + t267 * t256 + t296) * t283) * t303)) * t303 / 0.2e1 + (Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t291 ^ 2 + ((Icges(4,1) + Icges(5,1) * t290 ^ 2 + (-0.2e1 * Icges(5,4) * t290 + Icges(5,2) * t288) * t288) * t289 + 0.2e1 * (-Icges(5,5) * t290 + Icges(5,6) * t288 + Icges(4,4)) * t291) * t289) * t293 / 0.2e1;
T = t1;
