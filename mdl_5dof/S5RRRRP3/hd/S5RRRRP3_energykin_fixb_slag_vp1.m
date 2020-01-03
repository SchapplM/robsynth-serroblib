% Calculate kinetic energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:05
% EndTime: 2019-12-31 21:49:06
% DurationCPUTime: 0.80s
% Computational Cost: add. (684->96), mult. (504->154), div. (0->0), fcn. (399->8), ass. (0->68)
t344 = Icges(5,4) - Icges(6,5);
t343 = Icges(5,1) + Icges(6,1);
t342 = Icges(5,2) + Icges(6,3);
t280 = cos(qJ(4));
t341 = t344 * t280;
t278 = sin(qJ(4));
t340 = t344 * t278;
t339 = Icges(6,4) + Icges(5,5);
t338 = Icges(5,6) - Icges(6,6);
t337 = t342 * t278 - t341;
t336 = t343 * t280 - t340;
t335 = rSges(6,1) + pkin(4);
t334 = rSges(6,3) + qJ(5);
t333 = Icges(6,2) + Icges(5,3);
t277 = qJ(1) + qJ(2);
t275 = qJ(3) + t277;
t269 = sin(t275);
t270 = cos(t275);
t332 = t337 * t269 + t338 * t270;
t331 = -t338 * t269 + t337 * t270;
t330 = -t336 * t269 + t339 * t270;
t329 = t339 * t269 + t336 * t270;
t328 = -t342 * t280 - t340;
t327 = t343 * t278 + t341;
t326 = -t338 * t278 + t339 * t280;
t325 = t334 * t278 + t335 * t280;
t324 = t326 * t269 - t333 * t270;
t323 = t333 * t269 + t326 * t270;
t322 = t339 * t278 + t338 * t280;
t321 = t328 * t278 + t327 * t280;
t320 = t331 * t278 + t329 * t280;
t319 = -t332 * t278 + t330 * t280;
t276 = qJD(1) + qJD(2);
t314 = pkin(2) * t276;
t313 = pkin(1) * qJD(1);
t308 = -rSges(6,2) * t270 + t325 * t269;
t307 = rSges(6,2) * t269 + t325 * t270;
t281 = cos(qJ(1));
t271 = t281 * t313;
t274 = cos(t277);
t306 = t274 * t314 + t271;
t305 = qJD(4) * t269;
t304 = qJD(4) * t270;
t279 = sin(qJ(1));
t303 = t279 * t313;
t272 = qJD(3) + t276;
t302 = t272 * (pkin(3) * t270 + pkin(8) * t269) + t306;
t299 = rSges(5,1) * t280 - rSges(5,2) * t278;
t273 = sin(t277);
t284 = -t273 * t314 - t303;
t283 = t334 * qJD(4) * t280 + (-t335 * qJD(4) + qJD(5)) * t278;
t268 = rSges(2,1) * t281 - rSges(2,2) * t279;
t267 = rSges(2,1) * t279 + rSges(2,2) * t281;
t266 = rSges(5,1) * t278 + rSges(5,2) * t280;
t256 = pkin(3) * t269 - pkin(8) * t270;
t252 = t271 + t276 * (rSges(3,1) * t274 - rSges(3,2) * t273);
t251 = -t303 - t276 * (rSges(3,1) * t273 + rSges(3,2) * t274);
t250 = rSges(5,3) * t269 + t299 * t270;
t248 = -rSges(5,3) * t270 + t299 * t269;
t234 = t272 * (rSges(4,1) * t270 - rSges(4,2) * t269) + t306;
t233 = -t272 * (rSges(4,1) * t269 + rSges(4,2) * t270) + t284;
t232 = (t248 * t269 + t250 * t270) * qJD(4);
t231 = t250 * t272 - t266 * t305 + t302;
t230 = -t266 * t304 + (-t248 - t256) * t272 + t284;
t229 = -qJD(5) * t280 + (t308 * t269 + t307 * t270) * qJD(4);
t228 = t283 * t269 + t307 * t272 + t302;
t227 = (-t256 - t308) * t272 + t283 * t270 + t284;
t1 = m(3) * (t251 ^ 2 + t252 ^ 2) / 0.2e1 + t276 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t233 ^ 2 + t234 ^ 2) / 0.2e1 + t272 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + m(6) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + ((t327 * t278 - t328 * t280) * t272 + ((t330 * t278 + t332 * t280) * t270 + (t329 * t278 - t331 * t280) * t269) * qJD(4)) * t272 / 0.2e1 + (m(2) * (t267 ^ 2 + t268 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t322 * t269 + t321 * t270) * t272 + (t323 * t269 ^ 2 + (t319 * t270 + (t320 - t324) * t269) * t270) * qJD(4)) * t305 / 0.2e1 - ((t321 * t269 - t322 * t270) * t272 + (t324 * t270 ^ 2 + (t320 * t269 + (t319 - t323) * t270) * t269) * qJD(4)) * t304 / 0.2e1;
T = t1;
