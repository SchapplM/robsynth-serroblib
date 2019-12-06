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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:21:58
% EndTime: 2019-12-05 18:21:59
% DurationCPUTime: 0.96s
% Computational Cost: add. (695->101), mult. (504->150), div. (0->0), fcn. (398->8), ass. (0->66)
t345 = Icges(5,4) + Icges(6,4);
t344 = Icges(5,1) + Icges(6,1);
t343 = Icges(5,2) + Icges(6,2);
t278 = cos(qJ(4));
t342 = t345 * t278;
t276 = sin(qJ(4));
t341 = t345 * t276;
t340 = Icges(5,5) + Icges(6,5);
t339 = Icges(5,6) + Icges(6,6);
t338 = -t343 * t276 + t342;
t337 = t278 * t344 - t341;
t336 = rSges(6,1) + pkin(4);
t335 = Icges(5,3) + Icges(6,3);
t274 = qJ(1) + qJ(2);
t269 = pkin(8) + t274;
t266 = sin(t269);
t267 = cos(t269);
t334 = -t338 * t266 + t339 * t267;
t333 = t339 * t266 + t338 * t267;
t332 = -t337 * t266 + t340 * t267;
t331 = t340 * t266 + t337 * t267;
t330 = t278 * t343 + t341;
t329 = t344 * t276 + t342;
t328 = -t339 * t276 + t340 * t278;
t327 = rSges(6,3) + qJ(5);
t326 = -rSges(6,2) * t276 + t336 * t278;
t325 = -t328 * t266 + t335 * t267;
t324 = t335 * t266 + t328 * t267;
t323 = t340 * t276 + t339 * t278;
t322 = t330 * t276 - t329 * t278;
t321 = t333 * t276 - t331 * t278;
t320 = -t334 * t276 + t332 * t278;
t270 = sin(t274);
t314 = pkin(2) * t270;
t271 = cos(t274);
t313 = pkin(2) * t271;
t310 = pkin(1) * qJD(1);
t305 = t326 * t266 - t327 * t267;
t304 = t327 * t266 + t326 * t267;
t303 = qJD(4) * t266;
t302 = qJD(4) * t267;
t277 = sin(qJ(1));
t301 = t277 * t310;
t279 = cos(qJ(1));
t300 = t279 * t310;
t297 = -pkin(3) * t267 - pkin(7) * t266 - t313;
t296 = rSges(6,2) * t278 + t336 * t276;
t273 = qJD(1) + qJD(2);
t295 = t273 * (-pkin(3) * t266 + pkin(7) * t267) - t301;
t294 = rSges(5,1) * t278 - rSges(5,2) * t276;
t265 = rSges(2,1) * t279 - rSges(2,2) * t277;
t264 = -rSges(2,1) * t277 - rSges(2,2) * t279;
t263 = rSges(5,1) * t276 + rSges(5,2) * t278;
t253 = -t300 - t273 * (rSges(3,1) * t271 - rSges(3,2) * t270);
t252 = -t301 + t273 * (-rSges(3,1) * t270 - rSges(3,2) * t271);
t251 = rSges(5,3) * t266 + t294 * t267;
t249 = rSges(5,3) * t267 - t294 * t266;
t235 = -t300 + (-rSges(4,1) * t267 + rSges(4,2) * t266 - t313) * t273;
t234 = -t301 + (-rSges(4,1) * t266 - rSges(4,2) * t267 - t314) * t273;
t231 = qJD(3) + (-t249 * t266 + t251 * t267) * qJD(4);
t230 = -t300 + t263 * t303 + (-t251 + t297) * t273;
t229 = -t263 * t302 + (t249 - t314) * t273 + t295;
t228 = -t300 + qJD(5) * t267 + t296 * t303 + (t297 - t304) * t273;
t227 = qJD(5) * t266 - t296 * t302 + (-t305 - t314) * t273 + t295;
t226 = qJD(3) + (t305 * t266 + t304 * t267) * qJD(4);
t1 = m(3) * (t252 ^ 2 + t253 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + m(5) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + m(6) * (t226 ^ 2 + t227 ^ 2 + t228 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t273 ^ 2 / 0.2e1 + ((t329 * t276 + t330 * t278) * t273 + ((t332 * t276 + t334 * t278) * t267 + (t331 * t276 + t333 * t278) * t266) * qJD(4)) * t273 / 0.2e1 + (m(2) * (t264 ^ 2 + t265 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t323 * t266 - t322 * t267) * t273 + (t324 * t266 ^ 2 + (t320 * t267 + (-t321 + t325) * t266) * t267) * qJD(4)) * t303 / 0.2e1 + ((t322 * t266 + t323 * t267) * t273 + (t325 * t267 ^ 2 + (t321 * t266 + (-t320 + t324) * t267) * t266) * qJD(4)) * t302 / 0.2e1;
T = t1;
