% Calculate kinetic energy for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:21
% EndTime: 2019-12-31 16:57:22
% DurationCPUTime: 1.09s
% Computational Cost: add. (498->122), mult. (718->184), div. (0->0), fcn. (605->6), ass. (0->78)
t341 = Icges(4,4) - Icges(5,5);
t340 = Icges(4,1) + Icges(5,1);
t339 = Icges(4,2) + Icges(5,3);
t263 = qJ(2) + pkin(6);
t261 = cos(t263);
t338 = t341 * t261;
t260 = sin(t263);
t337 = t341 * t260;
t336 = Icges(5,4) + Icges(4,5);
t335 = Icges(4,6) - Icges(5,6);
t334 = t339 * t260 - t338;
t333 = t340 * t261 - t337;
t332 = rSges(5,1) + pkin(3);
t331 = rSges(5,3) + qJ(4);
t266 = sin(qJ(1));
t268 = cos(qJ(1));
t330 = t334 * t266 + t335 * t268;
t329 = -t335 * t266 + t334 * t268;
t328 = -t333 * t266 + t336 * t268;
t327 = t336 * t266 + t333 * t268;
t326 = -t339 * t261 - t337;
t325 = t340 * t260 + t338;
t324 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t265 = sin(qJ(2));
t267 = cos(qJ(2));
t323 = Icges(3,5) * t267 - Icges(3,6) * t265 - t335 * t260 + t336 * t261;
t322 = t331 * t260 + t332 * t261;
t321 = t323 * t266 - t324 * t268;
t320 = t324 * t266 + t323 * t268;
t319 = Icges(3,5) * t265 + Icges(3,6) * t267 + t336 * t260 + t335 * t261;
t308 = Icges(3,4) * t265;
t253 = Icges(3,2) * t267 + t308;
t307 = Icges(3,4) * t267;
t254 = Icges(3,1) * t265 + t307;
t318 = -t253 * t265 + t254 * t267 + t326 * t260 + t325 * t261;
t285 = -Icges(3,2) * t265 + t307;
t235 = Icges(3,6) * t266 + t285 * t268;
t288 = Icges(3,1) * t267 - t308;
t237 = Icges(3,5) * t266 + t288 * t268;
t317 = -t235 * t265 + t237 * t267 + t329 * t260 + t327 * t261;
t234 = -Icges(3,6) * t268 + t285 * t266;
t236 = -Icges(3,5) * t268 + t288 * t266;
t316 = t234 * t265 - t236 * t267 - t330 * t260 + t328 * t261;
t312 = pkin(2) * t265;
t310 = t267 * pkin(2);
t214 = -qJ(3) * t268 + t310 * t266;
t215 = qJ(3) * t266 + t310 * t268;
t297 = qJD(2) * t268;
t298 = qJD(2) * t266;
t302 = t214 * t298 + t215 * t297;
t258 = t266 * pkin(1) - t268 * pkin(5);
t301 = -t214 - t258;
t300 = -t268 * rSges(5,2) + t322 * t266;
t299 = t266 * rSges(5,2) + t322 * t268;
t251 = qJD(1) * (t268 * pkin(1) + t266 * pkin(5));
t294 = qJD(1) * t215 - qJD(3) * t268 + t251;
t293 = rSges(3,1) * t267 - rSges(3,2) * t265;
t292 = rSges(4,1) * t261 - rSges(4,2) * t260;
t289 = qJD(2) * (-t260 * rSges(4,1) - t261 * rSges(4,2) - t312);
t270 = (t331 * t261 - t312) * qJD(2) + (-t332 * qJD(2) + qJD(4)) * t260;
t262 = qJD(3) * t266;
t257 = t268 * rSges(2,1) - t266 * rSges(2,2);
t256 = t266 * rSges(2,1) + t268 * rSges(2,2);
t255 = t265 * rSges(3,1) + t267 * rSges(3,2);
t241 = t266 * rSges(3,3) + t293 * t268;
t240 = -t268 * rSges(3,3) + t293 * t266;
t231 = t266 * rSges(4,3) + t292 * t268;
t229 = -t268 * rSges(4,3) + t292 * t266;
t210 = qJD(1) * t241 - t255 * t298 + t251;
t209 = -t255 * t297 + (-t240 - t258) * qJD(1);
t208 = (t240 * t266 + t241 * t268) * qJD(2);
t207 = qJD(1) * t231 + t266 * t289 + t294;
t206 = t262 + t268 * t289 + (-t229 + t301) * qJD(1);
t205 = (t229 * t266 + t231 * t268) * qJD(2) + t302;
t204 = t299 * qJD(1) + t270 * t266 + t294;
t203 = t262 + t270 * t268 + (-t300 + t301) * qJD(1);
t202 = -qJD(4) * t261 + (t300 * t266 + t299 * t268) * qJD(2) + t302;
t1 = m(3) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(4) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(5) * (t202 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + (m(2) * (t256 ^ 2 + t257 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t267 * t234 - t265 * t236 + t328 * t260 + t330 * t261) * t268 + (t267 * t235 + t265 * t237 + t327 * t260 - t329 * t261) * t266) * qJD(2) + (t267 * t253 + t265 * t254 + t325 * t260 - t326 * t261) * qJD(1)) * qJD(1) / 0.2e1 + ((t320 * t266 ^ 2 + (t316 * t268 + (t317 - t321) * t266) * t268) * qJD(2) + (t319 * t266 + t318 * t268) * qJD(1)) * t298 / 0.2e1 - ((t321 * t268 ^ 2 + (t317 * t266 + (t316 - t320) * t268) * t266) * qJD(2) + (t318 * t266 - t319 * t268) * qJD(1)) * t297 / 0.2e1;
T = t1;
