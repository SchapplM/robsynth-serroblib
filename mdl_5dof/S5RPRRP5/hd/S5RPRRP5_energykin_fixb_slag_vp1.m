% Calculate kinetic energy for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:40
% EndTime: 2019-12-31 18:40:41
% DurationCPUTime: 0.91s
% Computational Cost: add. (658->99), mult. (506->152), div. (0->0), fcn. (399->8), ass. (0->66)
t344 = Icges(5,4) - Icges(6,5);
t343 = Icges(5,1) + Icges(6,1);
t342 = Icges(5,2) + Icges(6,3);
t281 = cos(qJ(4));
t341 = t344 * t281;
t279 = sin(qJ(4));
t340 = t344 * t279;
t339 = Icges(6,4) + Icges(5,5);
t338 = Icges(5,6) - Icges(6,6);
t337 = t342 * t279 - t341;
t336 = t343 * t281 - t340;
t335 = rSges(6,1) + pkin(4);
t334 = rSges(6,3) + qJ(5);
t333 = Icges(6,2) + Icges(5,3);
t278 = qJ(1) + pkin(8);
t276 = qJ(3) + t278;
t271 = sin(t276);
t272 = cos(t276);
t332 = t337 * t271 + t338 * t272;
t331 = -t338 * t271 + t337 * t272;
t330 = -t336 * t271 + t339 * t272;
t329 = t339 * t271 + t336 * t272;
t328 = -t342 * t281 - t340;
t327 = t343 * t279 + t341;
t326 = -t338 * t279 + t339 * t281;
t325 = t334 * t279 + t335 * t281;
t324 = t326 * t271 - t333 * t272;
t323 = t333 * t271 + t326 * t272;
t322 = t339 * t279 + t338 * t281;
t321 = t328 * t279 + t327 * t281;
t320 = t331 * t279 + t329 * t281;
t319 = -t332 * t279 + t330 * t281;
t280 = sin(qJ(1));
t314 = pkin(1) * t280;
t309 = -rSges(6,2) * t272 + t325 * t271;
t308 = rSges(6,2) * t271 + t325 * t272;
t282 = cos(qJ(1));
t273 = qJD(1) * t282 * pkin(1);
t275 = cos(t278);
t307 = qJD(1) * pkin(2) * t275 + t273;
t306 = qJD(4) * t271;
t305 = qJD(4) * t272;
t277 = qJD(1) + qJD(3);
t304 = t277 * (pkin(3) * t272 + pkin(7) * t271) + t307;
t301 = rSges(5,1) * t281 - rSges(5,2) * t279;
t274 = sin(t278);
t286 = (-pkin(2) * t274 - t314) * qJD(1);
t285 = t334 * qJD(4) * t281 + (-t335 * qJD(4) + qJD(5)) * t279;
t283 = qJD(2) ^ 2;
t269 = rSges(2,1) * t282 - rSges(2,2) * t280;
t268 = rSges(2,1) * t280 + rSges(2,2) * t282;
t267 = rSges(5,1) * t279 + rSges(5,2) * t281;
t258 = pkin(3) * t271 - pkin(7) * t272;
t254 = t273 + qJD(1) * (rSges(3,1) * t275 - rSges(3,2) * t274);
t253 = (-rSges(3,1) * t274 - rSges(3,2) * t275 - t314) * qJD(1);
t252 = rSges(5,3) * t271 + t301 * t272;
t250 = -rSges(5,3) * t272 + t301 * t271;
t236 = t277 * (rSges(4,1) * t272 - rSges(4,2) * t271) + t307;
t235 = -t277 * (rSges(4,1) * t271 + rSges(4,2) * t272) + t286;
t234 = qJD(2) + (t250 * t271 + t252 * t272) * qJD(4);
t233 = t252 * t277 - t267 * t306 + t304;
t232 = -t267 * t305 + (-t250 - t258) * t277 + t286;
t231 = t285 * t271 + t308 * t277 + t304;
t230 = t286 + (-t258 - t309) * t277 + t285 * t272;
t229 = -qJD(5) * t281 + qJD(2) + (t309 * t271 + t308 * t272) * qJD(4);
t1 = m(3) * (t253 ^ 2 + t254 ^ 2 + t283) / 0.2e1 + m(4) * (t235 ^ 2 + t236 ^ 2 + t283) / 0.2e1 + t277 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t232 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + m(6) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + ((t327 * t279 - t328 * t281) * t277 + ((t330 * t279 + t332 * t281) * t272 + (t329 * t279 - t331 * t281) * t271) * qJD(4)) * t277 / 0.2e1 + ((t322 * t271 + t321 * t272) * t277 + (t323 * t271 ^ 2 + (t319 * t272 + (t320 - t324) * t271) * t272) * qJD(4)) * t306 / 0.2e1 - ((t321 * t271 - t322 * t272) * t277 + (t324 * t272 ^ 2 + (t320 * t271 + (t319 - t323) * t272) * t271) * qJD(4)) * t305 / 0.2e1 + (m(2) * (t268 ^ 2 + t269 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
