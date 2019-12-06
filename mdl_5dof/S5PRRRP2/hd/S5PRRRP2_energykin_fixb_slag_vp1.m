% Calculate kinetic energy for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:36
% EndTime: 2019-12-05 16:41:37
% DurationCPUTime: 0.83s
% Computational Cost: add. (646->91), mult. (480->143), div. (0->0), fcn. (387->6), ass. (0->62)
t337 = Icges(5,4) - Icges(6,5);
t336 = Icges(5,1) + Icges(6,1);
t335 = Icges(5,2) + Icges(6,3);
t277 = cos(qJ(4));
t334 = t337 * t277;
t276 = sin(qJ(4));
t333 = t337 * t276;
t332 = Icges(6,4) + Icges(5,5);
t331 = Icges(5,6) - Icges(6,6);
t330 = t335 * t276 - t334;
t329 = t336 * t277 - t333;
t328 = rSges(6,1) + pkin(4);
t327 = rSges(6,3) + qJ(5);
t326 = Icges(6,2) + Icges(5,3);
t274 = pkin(8) + qJ(2);
t273 = qJ(3) + t274;
t269 = sin(t273);
t270 = cos(t273);
t325 = t330 * t269 + t331 * t270;
t324 = -t331 * t269 + t330 * t270;
t323 = -t329 * t269 + t332 * t270;
t322 = t332 * t269 + t329 * t270;
t321 = -t335 * t277 - t333;
t320 = t336 * t276 + t334;
t319 = -t331 * t276 + t332 * t277;
t318 = t327 * t276 + t328 * t277;
t317 = t319 * t269 - t270 * t326;
t316 = t269 * t326 + t319 * t270;
t315 = t332 * t276 + t331 * t277;
t314 = t276 * t321 + t277 * t320;
t313 = t276 * t324 + t277 * t322;
t312 = -t276 * t325 + t277 * t323;
t308 = pkin(2) * qJD(2);
t303 = -rSges(6,2) * t270 + t318 * t269;
t302 = rSges(6,2) * t269 + t318 * t270;
t272 = cos(t274);
t268 = t272 * t308;
t275 = qJD(2) + qJD(3);
t301 = t275 * (pkin(3) * t270 + pkin(7) * t269) + t268;
t300 = qJD(4) * t269;
t299 = qJD(4) * t270;
t271 = sin(t274);
t298 = t271 * t308;
t295 = rSges(5,1) * t277 - rSges(5,2) * t276;
t280 = t327 * qJD(4) * t277 + (-t328 * qJD(4) + qJD(5)) * t276;
t279 = qJD(1) ^ 2;
t278 = qJD(2) ^ 2;
t267 = rSges(5,1) * t276 + rSges(5,2) * t277;
t258 = rSges(3,1) * t272 - rSges(3,2) * t271;
t257 = rSges(3,1) * t271 + rSges(3,2) * t272;
t256 = pkin(3) * t269 - pkin(7) * t270;
t252 = t268 + t275 * (rSges(4,1) * t270 - rSges(4,2) * t269);
t251 = -t298 - t275 * (rSges(4,1) * t269 + rSges(4,2) * t270);
t250 = rSges(5,3) * t269 + t295 * t270;
t248 = -rSges(5,3) * t270 + t295 * t269;
t234 = qJD(1) + (t248 * t269 + t250 * t270) * qJD(4);
t233 = t250 * t275 - t267 * t300 + t301;
t232 = -t298 - t267 * t299 + (-t248 - t256) * t275;
t231 = t280 * t269 + t302 * t275 + t301;
t230 = -t298 + (-t256 - t303) * t275 + t280 * t270;
t229 = -qJD(5) * t277 + qJD(1) + (t303 * t269 + t302 * t270) * qJD(4);
t1 = m(2) * t279 / 0.2e1 + m(3) * (t279 + (t257 ^ 2 + t258 ^ 2) * t278) / 0.2e1 + t278 * Icges(3,3) / 0.2e1 + m(4) * (t251 ^ 2 + t252 ^ 2 + t279) / 0.2e1 + t275 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t232 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + m(6) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + ((t320 * t276 - t321 * t277) * t275 + ((t276 * t323 + t277 * t325) * t270 + (t276 * t322 - t277 * t324) * t269) * qJD(4)) * t275 / 0.2e1 + ((t315 * t269 + t314 * t270) * t275 + (t316 * t269 ^ 2 + (t312 * t270 + (t313 - t317) * t269) * t270) * qJD(4)) * t300 / 0.2e1 - ((t314 * t269 - t315 * t270) * t275 + (t317 * t270 ^ 2 + (t313 * t269 + (t312 - t316) * t270) * t269) * qJD(4)) * t299 / 0.2e1;
T = t1;
