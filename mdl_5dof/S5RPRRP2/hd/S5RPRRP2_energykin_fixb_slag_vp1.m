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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:27:49
% EndTime: 2022-01-23 09:27:50
% DurationCPUTime: 0.88s
% Computational Cost: add. (687->99), mult. (505->152), div. (0->0), fcn. (398->8), ass. (0->66)
t341 = Icges(5,4) + Icges(6,4);
t340 = Icges(5,1) + Icges(6,1);
t339 = Icges(5,2) + Icges(6,2);
t277 = cos(qJ(4));
t338 = t341 * t277;
t275 = sin(qJ(4));
t337 = t341 * t275;
t336 = Icges(5,5) + Icges(6,5);
t335 = Icges(5,6) + Icges(6,6);
t334 = -t339 * t275 + t338;
t333 = t340 * t277 - t337;
t332 = rSges(6,1) + pkin(4);
t331 = Icges(5,3) + Icges(6,3);
t273 = qJ(1) + pkin(8);
t271 = qJ(3) + t273;
t265 = sin(t271);
t266 = cos(t271);
t330 = t334 * t265 - t335 * t266;
t329 = t335 * t265 + t334 * t266;
t328 = -t333 * t265 + t336 * t266;
t327 = t336 * t265 + t333 * t266;
t326 = t339 * t277 + t337;
t325 = t340 * t275 + t338;
t324 = -t335 * t275 + t336 * t277;
t323 = rSges(6,3) + qJ(5);
t322 = -rSges(6,2) * t275 + t332 * t277;
t321 = t324 * t265 - t331 * t266;
t320 = t331 * t265 + t324 * t266;
t319 = t336 * t275 + t335 * t277;
t318 = -t326 * t275 + t325 * t277;
t317 = -t329 * t275 + t327 * t277;
t316 = t330 * t275 + t328 * t277;
t276 = sin(qJ(1));
t311 = pkin(1) * t276;
t304 = t322 * t265 - t323 * t266;
t303 = t323 * t265 + t322 * t266;
t278 = cos(qJ(1));
t268 = qJD(1) * t278 * pkin(1);
t270 = cos(t273);
t302 = qJD(1) * pkin(2) * t270 + t268;
t301 = qJD(4) * t265;
t300 = qJD(4) * t266;
t272 = qJD(1) + qJD(3);
t299 = t272 * (pkin(3) * t266 + pkin(7) * t265) + t302;
t296 = rSges(5,1) * t277 - rSges(5,2) * t275;
t294 = qJD(4) * (-rSges(6,2) * t277 - t332 * t275);
t269 = sin(t273);
t281 = (-pkin(2) * t269 - t311) * qJD(1);
t279 = qJD(2) ^ 2;
t263 = rSges(2,1) * t278 - rSges(2,2) * t276;
t262 = rSges(2,1) * t276 + rSges(2,2) * t278;
t261 = rSges(5,1) * t275 + rSges(5,2) * t277;
t253 = pkin(3) * t265 - pkin(7) * t266;
t251 = t268 + qJD(1) * (rSges(3,1) * t270 - rSges(3,2) * t269);
t250 = (-rSges(3,1) * t269 - rSges(3,2) * t270 - t311) * qJD(1);
t249 = rSges(5,3) * t265 + t296 * t266;
t247 = -rSges(5,3) * t266 + t296 * t265;
t233 = t272 * (rSges(4,1) * t266 - rSges(4,2) * t265) + t302;
t232 = -t272 * (rSges(4,1) * t265 + rSges(4,2) * t266) + t281;
t229 = qJD(2) + (t247 * t265 + t249 * t266) * qJD(4);
t228 = t249 * t272 - t261 * t301 + t299;
t227 = -t261 * t300 + (-t247 - t253) * t272 + t281;
t226 = -qJD(5) * t266 + t265 * t294 + t303 * t272 + t299;
t225 = qJD(5) * t265 + t266 * t294 + t281 + (-t253 - t304) * t272;
t224 = qJD(2) + (t304 * t265 + t303 * t266) * qJD(4);
t1 = m(3) * (t250 ^ 2 + t251 ^ 2 + t279) / 0.2e1 + m(4) * (t232 ^ 2 + t233 ^ 2 + t279) / 0.2e1 + t272 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + m(6) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + ((t325 * t275 + t326 * t277) * t272 + ((t328 * t275 - t330 * t277) * t266 + (t327 * t275 + t329 * t277) * t265) * qJD(4)) * t272 / 0.2e1 + ((t319 * t265 + t318 * t266) * t272 + (t320 * t265 ^ 2 + (t316 * t266 + (t317 - t321) * t265) * t266) * qJD(4)) * t301 / 0.2e1 - ((t318 * t265 - t319 * t266) * t272 + (t321 * t266 ^ 2 + (t317 * t265 + (t316 - t320) * t266) * t265) * qJD(4)) * t300 / 0.2e1 + (m(2) * (t262 ^ 2 + t263 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
