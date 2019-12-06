% Calculate kinetic energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:45
% EndTime: 2019-12-05 16:39:46
% DurationCPUTime: 0.88s
% Computational Cost: add. (675->91), mult. (479->143), div. (0->0), fcn. (386->6), ass. (0->62)
t334 = Icges(5,4) + Icges(6,4);
t333 = Icges(5,1) + Icges(6,1);
t332 = Icges(5,2) + Icges(6,2);
t273 = cos(qJ(4));
t331 = t334 * t273;
t272 = sin(qJ(4));
t330 = t334 * t272;
t329 = Icges(5,5) + Icges(6,5);
t328 = Icges(5,6) + Icges(6,6);
t327 = -t332 * t272 + t331;
t326 = t333 * t273 - t330;
t325 = rSges(6,1) + pkin(4);
t324 = Icges(5,3) + Icges(6,3);
t269 = pkin(8) + qJ(2);
t268 = qJ(3) + t269;
t263 = sin(t268);
t264 = cos(t268);
t323 = t327 * t263 - t328 * t264;
t322 = t328 * t263 + t327 * t264;
t321 = -t326 * t263 + t329 * t264;
t320 = t329 * t263 + t326 * t264;
t319 = t332 * t273 + t330;
t318 = t333 * t272 + t331;
t317 = -t328 * t272 + t329 * t273;
t316 = rSges(6,3) + qJ(5);
t315 = -rSges(6,2) * t272 + t325 * t273;
t314 = t317 * t263 - t324 * t264;
t313 = t324 * t263 + t317 * t264;
t312 = t329 * t272 + t328 * t273;
t311 = -t319 * t272 + t318 * t273;
t310 = -t322 * t272 + t320 * t273;
t309 = t323 * t272 + t321 * t273;
t303 = pkin(2) * qJD(2);
t298 = t315 * t263 - t316 * t264;
t297 = t316 * t263 + t315 * t264;
t267 = cos(t269);
t262 = t267 * t303;
t270 = qJD(2) + qJD(3);
t296 = t270 * (pkin(3) * t264 + pkin(7) * t263) + t262;
t295 = qJD(4) * t263;
t294 = qJD(4) * t264;
t266 = sin(t269);
t293 = t266 * t303;
t290 = rSges(5,1) * t273 - rSges(5,2) * t272;
t288 = qJD(4) * (-rSges(6,2) * t273 - t325 * t272);
t275 = qJD(1) ^ 2;
t274 = qJD(2) ^ 2;
t261 = rSges(5,1) * t272 + rSges(5,2) * t273;
t253 = rSges(3,1) * t267 - rSges(3,2) * t266;
t252 = rSges(3,1) * t266 + rSges(3,2) * t267;
t251 = pkin(3) * t263 - pkin(7) * t264;
t249 = t262 + t270 * (rSges(4,1) * t264 - rSges(4,2) * t263);
t248 = -t293 - t270 * (rSges(4,1) * t263 + rSges(4,2) * t264);
t247 = rSges(5,3) * t263 + t290 * t264;
t245 = -rSges(5,3) * t264 + t290 * t263;
t229 = qJD(1) + (t245 * t263 + t247 * t264) * qJD(4);
t228 = t247 * t270 - t261 * t295 + t296;
t227 = -t293 - t261 * t294 + (-t245 - t251) * t270;
t226 = -qJD(5) * t264 + t263 * t288 + t297 * t270 + t296;
t225 = -t293 + qJD(5) * t263 + t264 * t288 + (-t251 - t298) * t270;
t224 = qJD(1) + (t298 * t263 + t297 * t264) * qJD(4);
t1 = m(2) * t275 / 0.2e1 + m(3) * (t275 + (t252 ^ 2 + t253 ^ 2) * t274) / 0.2e1 + t274 * Icges(3,3) / 0.2e1 + m(4) * (t248 ^ 2 + t249 ^ 2 + t275) / 0.2e1 + t270 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + m(6) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + ((t318 * t272 + t319 * t273) * t270 + ((t321 * t272 - t323 * t273) * t264 + (t320 * t272 + t322 * t273) * t263) * qJD(4)) * t270 / 0.2e1 + ((t312 * t263 + t311 * t264) * t270 + (t313 * t263 ^ 2 + (t309 * t264 + (t310 - t314) * t263) * t264) * qJD(4)) * t295 / 0.2e1 - ((t311 * t263 - t312 * t264) * t270 + (t314 * t264 ^ 2 + (t310 * t263 + (t309 - t313) * t264) * t263) * qJD(4)) * t294 / 0.2e1;
T = t1;
