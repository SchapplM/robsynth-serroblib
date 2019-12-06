% Calculate kinetic energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:48
% EndTime: 2019-12-05 18:31:49
% DurationCPUTime: 0.58s
% Computational Cost: add. (803->130), mult. (514->210), div. (0->0), fcn. (408->10), ass. (0->81)
t263 = qJ(1) + qJ(2);
t255 = pkin(9) + t263;
t252 = sin(t255);
t253 = cos(t255);
t262 = qJ(4) + qJ(5);
t256 = sin(t262);
t258 = cos(t262);
t293 = Icges(6,4) * t258;
t279 = -Icges(6,2) * t256 + t293;
t222 = Icges(6,6) * t253 - t279 * t252;
t223 = Icges(6,6) * t252 + t279 * t253;
t294 = Icges(6,4) * t256;
t281 = Icges(6,1) * t258 - t294;
t224 = Icges(6,5) * t253 - t281 * t252;
t225 = Icges(6,5) * t252 + t281 * t253;
t290 = qJD(4) + qJD(5);
t239 = t290 * t252;
t240 = t290 * t253;
t243 = Icges(6,2) * t258 + t294;
t244 = Icges(6,1) * t256 + t293;
t261 = qJD(1) + qJD(2);
t305 = (t243 * t256 - t244 * t258) * t261 + (t222 * t256 - t224 * t258) * t240 + (t223 * t256 - t225 * t258) * t239;
t257 = sin(t263);
t301 = pkin(2) * t257;
t259 = cos(t263);
t300 = pkin(2) * t259;
t266 = cos(qJ(4));
t299 = pkin(4) * t266;
t297 = pkin(1) * qJD(1);
t264 = sin(qJ(4));
t296 = Icges(5,4) * t264;
t295 = Icges(5,4) * t266;
t292 = qJD(4) * t252;
t291 = qJD(4) * t253;
t289 = pkin(4) * qJD(4) * t264;
t265 = sin(qJ(1));
t288 = t265 * t297;
t267 = cos(qJ(1));
t287 = t267 * t297;
t286 = -pkin(3) * t253 - pkin(7) * t252 - t300;
t285 = t261 * (-pkin(3) * t252 + pkin(7) * t253) - t288;
t284 = rSges(5,1) * t266 - rSges(5,2) * t264;
t283 = rSges(6,1) * t258 - rSges(6,2) * t256;
t282 = Icges(5,1) * t266 - t296;
t280 = -Icges(5,2) * t264 + t295;
t278 = Icges(5,5) * t266 - Icges(5,6) * t264;
t277 = Icges(6,5) * t258 - Icges(6,6) * t256;
t230 = Icges(5,6) * t253 - t280 * t252;
t232 = Icges(5,5) * t253 - t282 * t252;
t274 = -t230 * t264 + t232 * t266;
t231 = Icges(5,6) * t252 + t280 * t253;
t233 = Icges(5,5) * t252 + t282 * t253;
t273 = t231 * t264 - t233 * t266;
t247 = Icges(5,2) * t266 + t296;
t248 = Icges(5,1) * t264 + t295;
t271 = t247 * t264 - t248 * t266;
t270 = (Icges(6,3) * t253 - t277 * t252) * t240 + (Icges(6,3) * t252 + t277 * t253) * t239 + (Icges(6,5) * t256 + Icges(6,6) * t258) * t261;
t251 = rSges(2,1) * t267 - rSges(2,2) * t265;
t250 = -rSges(2,1) * t265 - rSges(2,2) * t267;
t249 = rSges(5,1) * t264 + rSges(5,2) * t266;
t246 = Icges(5,5) * t264 + Icges(5,6) * t266;
t245 = rSges(6,1) * t256 + rSges(6,2) * t258;
t237 = -t287 - t261 * (rSges(3,1) * t259 - rSges(3,2) * t257);
t236 = -t288 + t261 * (-rSges(3,1) * t257 - rSges(3,2) * t259);
t235 = rSges(5,3) * t252 + t284 * t253;
t234 = rSges(5,3) * t253 - t284 * t252;
t229 = Icges(5,3) * t252 + t278 * t253;
t228 = Icges(5,3) * t253 - t278 * t252;
t227 = rSges(6,3) * t252 + t283 * t253;
t226 = rSges(6,3) * t253 - t283 * t252;
t219 = -t287 + (-rSges(4,1) * t253 + rSges(4,2) * t252 - t300) * t261;
t218 = -t288 + (-rSges(4,1) * t252 - rSges(4,2) * t253 - t301) * t261;
t217 = pkin(8) * t252 + t299 * t253;
t216 = pkin(8) * t253 - t299 * t252;
t215 = qJD(3) + (-t234 * t252 + t235 * t253) * qJD(4);
t214 = -t287 + t249 * t292 + (-t235 + t286) * t261;
t213 = -t249 * t291 + (t234 - t301) * t261 + t285;
t212 = t252 * t289 - t287 + t239 * t245 + (-t217 - t227 + t286) * t261;
t211 = -t253 * t289 - t240 * t245 + (t216 + t226 - t301) * t261 + t285;
t210 = -t226 * t239 + t227 * t240 + qJD(3) + (-t216 * t252 + t217 * t253) * qJD(4);
t1 = m(3) * (t236 ^ 2 + t237 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + ((t253 * t246 + t271 * t252) * t261 + (t253 ^ 2 * t228 + (t273 * t252 + (t229 - t274) * t253) * t252) * qJD(4)) * t291 / 0.2e1 + ((t252 * t246 - t271 * t253) * t261 + (t252 ^ 2 * t229 + (t274 * t253 + (t228 - t273) * t252) * t253) * qJD(4)) * t292 / 0.2e1 + m(6) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + t240 * (t305 * t252 + t270 * t253) / 0.2e1 + t239 * (t270 * t252 - t305 * t253) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t261 ^ 2 / 0.2e1 + (((t230 * t266 + t232 * t264) * t253 + (t231 * t266 + t233 * t264) * t252) * qJD(4) + (t222 * t258 + t224 * t256) * t240 + (t223 * t258 + t225 * t256) * t239 + (t243 * t258 + t244 * t256 + t247 * t266 + t248 * t264) * t261) * t261 / 0.2e1 + (m(2) * (t250 ^ 2 + t251 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
