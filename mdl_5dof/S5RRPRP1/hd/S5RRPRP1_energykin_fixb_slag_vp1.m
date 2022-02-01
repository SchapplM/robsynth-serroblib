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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:19:31
% EndTime: 2022-01-20 10:19:32
% DurationCPUTime: 0.91s
% Computational Cost: add. (695->99), mult. (504->152), div. (0->0), fcn. (398->8), ass. (0->66)
t344 = Icges(5,4) + Icges(6,4);
t343 = Icges(5,1) + Icges(6,1);
t342 = Icges(5,2) + Icges(6,2);
t278 = cos(qJ(4));
t341 = t344 * t278;
t276 = sin(qJ(4));
t340 = t344 * t276;
t339 = Icges(5,5) + Icges(6,5);
t338 = Icges(5,6) + Icges(6,6);
t337 = -t342 * t276 + t341;
t336 = t343 * t278 - t340;
t335 = rSges(6,1) + pkin(4);
t334 = Icges(5,3) + Icges(6,3);
t274 = qJ(1) + qJ(2);
t269 = pkin(8) + t274;
t265 = sin(t269);
t266 = cos(t269);
t333 = t337 * t265 - t338 * t266;
t332 = t338 * t265 + t337 * t266;
t331 = -t336 * t265 + t339 * t266;
t330 = t339 * t265 + t336 * t266;
t329 = t342 * t278 + t340;
t328 = t343 * t276 + t341;
t327 = -t338 * t276 + t339 * t278;
t326 = rSges(6,3) + qJ(5);
t325 = -rSges(6,2) * t276 + t335 * t278;
t324 = t327 * t265 - t334 * t266;
t323 = t334 * t265 + t327 * t266;
t322 = t339 * t276 + t338 * t278;
t321 = -t329 * t276 + t328 * t278;
t320 = -t332 * t276 + t330 * t278;
t319 = t333 * t276 + t331 * t278;
t270 = sin(t274);
t313 = pkin(2) * t270;
t310 = pkin(1) * qJD(1);
t305 = t325 * t265 - t326 * t266;
t304 = t326 * t265 + t325 * t266;
t279 = cos(qJ(1));
t268 = t279 * t310;
t271 = cos(t274);
t273 = qJD(1) + qJD(2);
t303 = t273 * pkin(2) * t271 + t268;
t302 = qJD(4) * t265;
t301 = qJD(4) * t266;
t277 = sin(qJ(1));
t300 = t277 * t310;
t299 = t273 * (pkin(3) * t266 + pkin(7) * t265) + t303;
t296 = -pkin(3) * t265 + pkin(7) * t266 - t313;
t295 = rSges(5,1) * t278 - rSges(5,2) * t276;
t293 = qJD(4) * (-rSges(6,2) * t278 - t335 * t276);
t264 = rSges(2,1) * t279 - rSges(2,2) * t277;
t263 = rSges(2,1) * t277 + rSges(2,2) * t279;
t262 = rSges(5,1) * t276 + rSges(5,2) * t278;
t251 = t268 + t273 * (rSges(3,1) * t271 - rSges(3,2) * t270);
t250 = -t300 - t273 * (rSges(3,1) * t270 + rSges(3,2) * t271);
t249 = rSges(5,3) * t265 + t295 * t266;
t247 = -rSges(5,3) * t266 + t295 * t265;
t233 = t273 * (rSges(4,1) * t266 - rSges(4,2) * t265) + t303;
t232 = -t300 + (-rSges(4,1) * t265 - rSges(4,2) * t266 - t313) * t273;
t229 = qJD(3) + (t247 * t265 + t249 * t266) * qJD(4);
t228 = t249 * t273 - t262 * t302 + t299;
t227 = -t300 - t262 * t301 + (-t247 + t296) * t273;
t226 = -qJD(5) * t266 + t265 * t293 + t304 * t273 + t299;
t225 = -t300 + qJD(5) * t265 + t266 * t293 + (t296 - t305) * t273;
t224 = qJD(3) + (t305 * t265 + t304 * t266) * qJD(4);
t1 = m(3) * (t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t232 ^ 2 + t233 ^ 2) / 0.2e1 + m(5) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + m(6) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t273 ^ 2 / 0.2e1 + ((t328 * t276 + t329 * t278) * t273 + ((t331 * t276 - t333 * t278) * t266 + (t330 * t276 + t332 * t278) * t265) * qJD(4)) * t273 / 0.2e1 + (m(2) * (t263 ^ 2 + t264 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t322 * t265 + t321 * t266) * t273 + (t323 * t265 ^ 2 + (t319 * t266 + (t320 - t324) * t265) * t266) * qJD(4)) * t302 / 0.2e1 - ((t321 * t265 - t322 * t266) * t273 + (t324 * t266 ^ 2 + (t320 * t265 + (t319 - t323) * t266) * t265) * qJD(4)) * t301 / 0.2e1;
T = t1;
