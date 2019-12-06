% Calculate kinetic energy for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:21
% EndTime: 2019-12-05 18:01:22
% DurationCPUTime: 0.91s
% Computational Cost: add. (687->99), mult. (505->151), div. (0->0), fcn. (398->8), ass. (0->66)
t342 = Icges(5,4) + Icges(6,4);
t341 = Icges(5,1) + Icges(6,1);
t340 = Icges(5,2) + Icges(6,2);
t277 = cos(qJ(4));
t339 = t342 * t277;
t275 = sin(qJ(4));
t338 = t342 * t275;
t337 = Icges(5,5) + Icges(6,5);
t336 = Icges(5,6) + Icges(6,6);
t335 = -t275 * t340 + t339;
t334 = t277 * t341 - t338;
t333 = rSges(6,1) + pkin(4);
t332 = Icges(5,3) + Icges(6,3);
t273 = qJ(1) + pkin(8);
t271 = qJ(3) + t273;
t266 = sin(t271);
t267 = cos(t271);
t331 = -t335 * t266 + t336 * t267;
t330 = t336 * t266 + t335 * t267;
t329 = -t334 * t266 + t337 * t267;
t328 = t337 * t266 + t334 * t267;
t327 = t277 * t340 + t338;
t326 = t275 * t341 + t339;
t325 = -t336 * t275 + t337 * t277;
t324 = rSges(6,3) + qJ(5);
t323 = -rSges(6,2) * t275 + t333 * t277;
t322 = -t325 * t266 + t332 * t267;
t321 = t332 * t266 + t325 * t267;
t320 = t337 * t275 + t336 * t277;
t319 = t327 * t275 - t326 * t277;
t318 = t330 * t275 - t328 * t277;
t317 = -t331 * t275 + t329 * t277;
t276 = sin(qJ(1));
t312 = pkin(1) * t276;
t278 = cos(qJ(1));
t311 = pkin(1) * t278;
t304 = -t323 * t266 + t324 * t267;
t303 = t324 * t266 + t323 * t267;
t302 = qJD(4) * t266;
t301 = qJD(4) * t267;
t298 = rSges(6,2) * t277 + t333 * t275;
t297 = rSges(5,1) * t277 - rSges(5,2) * t275;
t269 = sin(t273);
t283 = (-pkin(2) * t269 - t312) * qJD(1);
t270 = cos(t273);
t282 = (-pkin(2) * t270 - t311) * qJD(1);
t272 = qJD(1) + qJD(3);
t281 = t272 * (-pkin(3) * t266 + pkin(7) * t267) + t283;
t279 = qJD(2) ^ 2;
t265 = rSges(2,1) * t278 - rSges(2,2) * t276;
t264 = -rSges(2,1) * t276 - rSges(2,2) * t278;
t263 = rSges(5,1) * t275 + rSges(5,2) * t277;
t255 = pkin(3) * t267 + pkin(7) * t266;
t253 = (-rSges(3,1) * t270 + rSges(3,2) * t269 - t311) * qJD(1);
t252 = (-rSges(3,1) * t269 - rSges(3,2) * t270 - t312) * qJD(1);
t251 = rSges(5,3) * t266 + t297 * t267;
t249 = rSges(5,3) * t267 - t297 * t266;
t235 = -t272 * (rSges(4,1) * t267 - rSges(4,2) * t266) + t282;
t234 = t272 * (-rSges(4,1) * t266 - rSges(4,2) * t267) + t283;
t231 = qJD(2) + (-t249 * t266 + t251 * t267) * qJD(4);
t230 = t263 * t302 + (-t251 - t255) * t272 + t282;
t229 = t249 * t272 - t263 * t301 + t281;
t228 = qJD(5) * t267 + t298 * t302 + t282 + (-t255 - t303) * t272;
t227 = qJD(5) * t266 + t304 * t272 - t298 * t301 + t281;
t226 = qJD(2) + (-t304 * t266 + t303 * t267) * qJD(4);
t1 = m(3) * (t252 ^ 2 + t253 ^ 2 + t279) / 0.2e1 + m(4) * (t234 ^ 2 + t235 ^ 2 + t279) / 0.2e1 + t272 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + m(6) * (t226 ^ 2 + t227 ^ 2 + t228 ^ 2) / 0.2e1 + ((t326 * t275 + t327 * t277) * t272 + ((t329 * t275 + t331 * t277) * t267 + (t328 * t275 + t330 * t277) * t266) * qJD(4)) * t272 / 0.2e1 + ((t320 * t266 - t319 * t267) * t272 + (t321 * t266 ^ 2 + (t317 * t267 + (-t318 + t322) * t266) * t267) * qJD(4)) * t302 / 0.2e1 + ((t319 * t266 + t320 * t267) * t272 + (t322 * t267 ^ 2 + (t318 * t266 + (-t317 + t321) * t267) * t266) * qJD(4)) * t301 / 0.2e1 + (m(2) * (t264 ^ 2 + t265 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
