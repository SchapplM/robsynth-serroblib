% Calculate kinetic energy for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:36
% EndTime: 2019-12-31 17:16:37
% DurationCPUTime: 1.12s
% Computational Cost: add. (522->127), mult. (742->204), div. (0->0), fcn. (629->6), ass. (0->84)
t337 = Icges(4,4) - Icges(5,5);
t336 = Icges(4,1) + Icges(5,1);
t335 = Icges(4,2) + Icges(5,3);
t266 = qJ(2) + qJ(3);
t265 = cos(t266);
t334 = t337 * t265;
t264 = sin(t266);
t333 = t337 * t264;
t332 = Icges(5,4) + Icges(4,5);
t331 = Icges(4,6) - Icges(5,6);
t330 = t335 * t264 - t334;
t329 = t336 * t265 - t333;
t328 = rSges(5,1) + pkin(3);
t327 = rSges(5,3) + qJ(4);
t326 = Icges(5,2) + Icges(4,3);
t268 = sin(qJ(1));
t270 = cos(qJ(1));
t325 = -t330 * t268 - t331 * t270;
t324 = -t331 * t268 + t330 * t270;
t323 = t329 * t268 - t332 * t270;
t322 = t332 * t268 + t329 * t270;
t321 = -t335 * t265 - t333;
t320 = t336 * t264 + t334;
t319 = -t331 * t264 + t332 * t265;
t318 = t327 * t264 + t328 * t265;
t295 = qJD(2) + qJD(3);
t254 = t295 * t268;
t255 = t295 * t270;
t317 = (t325 * t264 - t323 * t265) * t255 + (t324 * t264 + t322 * t265) * t254 + (t321 * t264 + t320 * t265) * qJD(1);
t316 = (-t319 * t268 + t326 * t270) * t255 + (t326 * t268 + t319 * t270) * t254 + (t332 * t264 + t331 * t265) * qJD(1);
t269 = cos(qJ(2));
t311 = t269 * pkin(2);
t267 = sin(qJ(2));
t309 = Icges(3,4) * t267;
t308 = Icges(3,4) * t269;
t216 = -pkin(6) * t270 + t311 * t268;
t217 = pkin(6) * t268 + t311 * t270;
t296 = qJD(2) * t270;
t297 = qJD(2) * t268;
t303 = t216 * t297 + t217 * t296;
t253 = qJD(1) * (t270 * pkin(1) + t268 * pkin(5));
t302 = qJD(1) * t217 + t253;
t262 = t268 * pkin(1) - t270 * pkin(5);
t301 = -t216 - t262;
t300 = -t270 * rSges(5,2) + t318 * t268;
t299 = t268 * rSges(5,2) + t318 * t270;
t298 = t328 * t264 - t327 * t265;
t294 = pkin(2) * qJD(2) * t267;
t293 = rSges(3,1) * t269 - rSges(3,2) * t267;
t292 = rSges(4,1) * t265 - rSges(4,2) * t264;
t289 = Icges(3,1) * t269 - t309;
t286 = -Icges(3,2) * t267 + t308;
t283 = Icges(3,5) * t269 - Icges(3,6) * t267;
t236 = -Icges(3,6) * t270 + t286 * t268;
t238 = -Icges(3,5) * t270 + t289 * t268;
t280 = t236 * t267 - t238 * t269;
t237 = Icges(3,6) * t268 + t286 * t270;
t239 = Icges(3,5) * t268 + t289 * t270;
t279 = -t237 * t267 + t239 * t269;
t257 = Icges(3,2) * t269 + t309;
t258 = Icges(3,1) * t267 + t308;
t278 = -t257 * t267 + t258 * t269;
t277 = qJD(4) * t264 - t294;
t261 = t270 * rSges(2,1) - t268 * rSges(2,2);
t260 = t268 * rSges(2,1) + t270 * rSges(2,2);
t259 = t267 * rSges(3,1) + t269 * rSges(3,2);
t256 = Icges(3,5) * t267 + Icges(3,6) * t269;
t252 = t264 * rSges(4,1) + t265 * rSges(4,2);
t241 = t268 * rSges(3,3) + t293 * t270;
t240 = -t270 * rSges(3,3) + t293 * t268;
t235 = Icges(3,3) * t268 + t283 * t270;
t234 = -Icges(3,3) * t270 + t283 * t268;
t233 = t268 * rSges(4,3) + t292 * t270;
t231 = -t270 * rSges(4,3) + t292 * t268;
t212 = qJD(1) * t241 - t259 * t297 + t253;
t211 = -t259 * t296 + (-t240 - t262) * qJD(1);
t210 = (t240 * t268 + t241 * t270) * qJD(2);
t209 = qJD(1) * t233 - t254 * t252 - t268 * t294 + t302;
t208 = -t270 * t294 - t255 * t252 + (-t231 + t301) * qJD(1);
t207 = t254 * t231 + t255 * t233 + t303;
t206 = t299 * qJD(1) - t298 * t254 + t277 * t268 + t302;
t205 = t277 * t270 - t298 * t255 + (-t300 + t301) * qJD(1);
t204 = -qJD(4) * t265 + t300 * t254 + t299 * t255 + t303;
t1 = m(3) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + ((t268 * t256 + t278 * t270) * qJD(1) + (t268 ^ 2 * t235 + (t280 * t270 + (-t234 + t279) * t268) * t270) * qJD(2)) * t297 / 0.2e1 - ((-t270 * t256 + t278 * t268) * qJD(1) + (t270 ^ 2 * t234 + (t279 * t268 + (-t235 + t280) * t270) * t268) * qJD(2)) * t296 / 0.2e1 + m(4) * (t207 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + m(5) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + (t316 * t268 + t317 * t270) * t254 / 0.2e1 - (t317 * t268 - t316 * t270) * t255 / 0.2e1 + (m(2) * (t260 ^ 2 + t261 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t269 * t237 + t267 * t239) * t268 - (t269 * t236 + t267 * t238) * t270) * qJD(2) - (t323 * t264 + t325 * t265) * t255 + (t322 * t264 - t324 * t265) * t254 + (t269 * t257 + t267 * t258 + t320 * t264 - t321 * t265) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
