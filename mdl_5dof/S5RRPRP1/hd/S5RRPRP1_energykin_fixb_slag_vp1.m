% Calculate kinetic energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:53
% EndTime: 2020-01-03 11:58:53
% DurationCPUTime: 0.84s
% Computational Cost: add. (695->98), mult. (504->151), div. (0->0), fcn. (398->8), ass. (0->67)
t348 = Icges(5,4) + Icges(6,4);
t347 = Icges(5,1) + Icges(6,1);
t346 = Icges(5,2) + Icges(6,2);
t283 = cos(qJ(4));
t345 = t348 * t283;
t281 = sin(qJ(4));
t344 = t348 * t281;
t343 = Icges(5,5) + Icges(6,5);
t342 = Icges(5,6) + Icges(6,6);
t341 = t346 * t281 - t345;
t340 = t347 * t283 - t344;
t339 = rSges(6,1) + pkin(4);
t338 = -Icges(5,3) - Icges(6,3);
t279 = qJ(1) + qJ(2);
t274 = pkin(8) + t279;
t269 = sin(t274);
t270 = cos(t274);
t337 = t341 * t269 + t342 * t270;
t336 = -t342 * t269 + t341 * t270;
t335 = t340 * t269 - t343 * t270;
t334 = t343 * t269 + t340 * t270;
t333 = t346 * t283 + t344;
t332 = t347 * t281 + t345;
t331 = -t342 * t281 + t343 * t283;
t330 = rSges(6,3) + qJ(5);
t329 = -rSges(6,2) * t281 + t339 * t283;
t328 = t331 * t269 + t338 * t270;
t327 = t338 * t269 - t331 * t270;
t326 = -t343 * t281 - t342 * t283;
t325 = t333 * t281 - t332 * t283;
t324 = t336 * t281 + t334 * t283;
t323 = t337 * t281 + t335 * t283;
t278 = qJD(1) + qJD(2);
t317 = pkin(2) * t278;
t314 = pkin(1) * qJD(1);
t309 = t329 * t269 - t330 * t270;
t308 = t330 * t269 + t329 * t270;
t282 = sin(qJ(1));
t272 = t282 * t314;
t275 = sin(t279);
t307 = t275 * t317 + t272;
t284 = cos(qJ(1));
t273 = t284 * t314;
t276 = cos(t279);
t306 = t276 * t317 + t273;
t305 = qJD(4) * t269;
t304 = qJD(4) * t270;
t303 = t278 * (pkin(3) * t269 - pkin(7) * t270) + t307;
t300 = rSges(6,2) * t283 + t339 * t281;
t299 = rSges(5,1) * t283 - rSges(5,2) * t281;
t268 = -rSges(2,1) * t284 + rSges(2,2) * t282;
t267 = rSges(2,1) * t282 + rSges(2,2) * t284;
t266 = rSges(5,1) * t281 + rSges(5,2) * t283;
t256 = -pkin(3) * t270 - pkin(7) * t269;
t254 = t273 - t278 * (-rSges(3,1) * t276 + rSges(3,2) * t275);
t253 = t272 + t278 * (rSges(3,1) * t275 + rSges(3,2) * t276);
t252 = -rSges(5,3) * t269 - t299 * t270;
t250 = -rSges(5,3) * t270 + t299 * t269;
t236 = -t278 * (-rSges(4,1) * t270 + rSges(4,2) * t269) + t306;
t235 = t278 * (rSges(4,1) * t269 + rSges(4,2) * t270) + t307;
t232 = qJD(3) + (t250 * t269 - t252 * t270) * qJD(4);
t231 = -t266 * t305 + (-t252 - t256) * t278 + t306;
t230 = t250 * t278 + t266 * t304 + t303;
t229 = -qJD(5) * t270 - t300 * t305 + (-t256 + t308) * t278 + t306;
t228 = -qJD(5) * t269 + t309 * t278 + t300 * t304 + t303;
t227 = qJD(3) + (t309 * t269 + t308 * t270) * qJD(4);
t1 = m(3) * (t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + m(5) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + m(6) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t278 ^ 2 / 0.2e1 + ((t332 * t281 + t333 * t283) * t278 + ((-t335 * t281 + t337 * t283) * t270 + (t334 * t281 - t336 * t283) * t269) * qJD(4)) * t278 / 0.2e1 + (m(2) * (t267 ^ 2 + t268 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 - ((t326 * t269 + t325 * t270) * t278 + (t327 * t269 ^ 2 + (t323 * t270 + (-t324 + t328) * t269) * t270) * qJD(4)) * t305 / 0.2e1 - ((-t325 * t269 + t326 * t270) * t278 + (t328 * t270 ^ 2 + (t324 * t269 + (-t323 + t327) * t270) * t269) * qJD(4)) * t304 / 0.2e1;
T = t1;
