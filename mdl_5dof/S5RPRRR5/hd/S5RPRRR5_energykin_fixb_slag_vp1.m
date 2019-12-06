% Calculate kinetic energy for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:05
% EndTime: 2019-12-05 18:16:05
% DurationCPUTime: 0.61s
% Computational Cost: add. (795->128), mult. (515->211), div. (0->0), fcn. (408->10), ass. (0->81)
t261 = qJ(1) + pkin(9);
t257 = qJ(3) + t261;
t252 = sin(t257);
t253 = cos(t257);
t262 = qJ(4) + qJ(5);
t258 = sin(t262);
t259 = cos(t262);
t292 = Icges(6,4) * t259;
t282 = -Icges(6,2) * t258 + t292;
t222 = Icges(6,6) * t253 - t252 * t282;
t223 = Icges(6,6) * t252 + t253 * t282;
t293 = Icges(6,4) * t258;
t284 = Icges(6,1) * t259 - t293;
t224 = Icges(6,5) * t253 - t252 * t284;
t225 = Icges(6,5) * t252 + t253 * t284;
t289 = qJD(4) + qJD(5);
t239 = t289 * t252;
t240 = t289 * t253;
t243 = Icges(6,2) * t259 + t293;
t244 = Icges(6,1) * t258 + t292;
t260 = qJD(1) + qJD(3);
t302 = (t243 * t258 - t244 * t259) * t260 + (t222 * t258 - t224 * t259) * t240 + (t223 * t258 - t225 * t259) * t239;
t264 = sin(qJ(1));
t299 = pkin(1) * t264;
t266 = cos(qJ(1));
t298 = pkin(1) * t266;
t265 = cos(qJ(4));
t297 = pkin(4) * t265;
t263 = sin(qJ(4));
t295 = Icges(5,4) * t263;
t294 = Icges(5,4) * t265;
t291 = qJD(4) * t252;
t290 = qJD(4) * t253;
t288 = pkin(4) * qJD(4) * t263;
t287 = rSges(5,1) * t265 - rSges(5,2) * t263;
t286 = rSges(6,1) * t259 - rSges(6,2) * t258;
t285 = Icges(5,1) * t265 - t295;
t283 = -Icges(5,2) * t263 + t294;
t281 = Icges(5,5) * t265 - Icges(5,6) * t263;
t280 = Icges(6,5) * t259 - Icges(6,6) * t258;
t230 = Icges(5,6) * t253 - t252 * t283;
t232 = Icges(5,5) * t253 - t252 * t285;
t277 = -t230 * t263 + t232 * t265;
t231 = Icges(5,6) * t252 + t253 * t283;
t233 = Icges(5,5) * t252 + t253 * t285;
t276 = t231 * t263 - t233 * t265;
t247 = Icges(5,2) * t265 + t295;
t248 = Icges(5,1) * t263 + t294;
t274 = t247 * t263 - t248 * t265;
t255 = sin(t261);
t273 = (-pkin(2) * t255 - t299) * qJD(1);
t256 = cos(t261);
t272 = (-pkin(2) * t256 - t298) * qJD(1);
t271 = (Icges(6,3) * t253 - t252 * t280) * t240 + (Icges(6,3) * t252 + t253 * t280) * t239 + (Icges(6,5) * t258 + Icges(6,6) * t259) * t260;
t270 = t260 * (-pkin(3) * t252 + pkin(7) * t253) + t273;
t268 = qJD(2) ^ 2;
t251 = rSges(2,1) * t266 - rSges(2,2) * t264;
t250 = -rSges(2,1) * t264 - rSges(2,2) * t266;
t249 = rSges(5,1) * t263 + rSges(5,2) * t265;
t246 = Icges(5,5) * t263 + Icges(5,6) * t265;
t245 = rSges(6,1) * t258 + rSges(6,2) * t259;
t241 = pkin(3) * t253 + pkin(7) * t252;
t237 = (-rSges(3,1) * t256 + rSges(3,2) * t255 - t298) * qJD(1);
t236 = (-rSges(3,1) * t255 - rSges(3,2) * t256 - t299) * qJD(1);
t235 = rSges(5,3) * t252 + t253 * t287;
t234 = rSges(5,3) * t253 - t252 * t287;
t229 = Icges(5,3) * t252 + t253 * t281;
t228 = Icges(5,3) * t253 - t252 * t281;
t227 = rSges(6,3) * t252 + t253 * t286;
t226 = rSges(6,3) * t253 - t252 * t286;
t219 = -t260 * (rSges(4,1) * t253 - rSges(4,2) * t252) + t272;
t218 = t260 * (-rSges(4,1) * t252 - rSges(4,2) * t253) + t273;
t217 = pkin(8) * t252 + t253 * t297;
t216 = pkin(8) * t253 - t252 * t297;
t215 = qJD(2) + (-t234 * t252 + t235 * t253) * qJD(4);
t214 = t249 * t291 + (-t235 - t241) * t260 + t272;
t213 = t234 * t260 - t249 * t290 + t270;
t212 = t252 * t288 + t239 * t245 + t272 + (-t217 - t227 - t241) * t260;
t211 = -t253 * t288 - t240 * t245 + (t216 + t226) * t260 + t270;
t210 = -t226 * t239 + t227 * t240 + qJD(2) + (-t216 * t252 + t217 * t253) * qJD(4);
t1 = m(3) * (t236 ^ 2 + t237 ^ 2 + t268) / 0.2e1 + m(4) * (t218 ^ 2 + t219 ^ 2 + t268) / 0.2e1 + t260 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + ((t253 * t246 + t252 * t274) * t260 + (t253 ^ 2 * t228 + (t276 * t252 + (t229 - t277) * t253) * t252) * qJD(4)) * t290 / 0.2e1 + ((t252 * t246 - t253 * t274) * t260 + (t252 ^ 2 * t229 + (t277 * t253 + (t228 - t276) * t252) * t253) * qJD(4)) * t291 / 0.2e1 + m(6) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + t240 * (t302 * t252 + t271 * t253) / 0.2e1 + t239 * (t271 * t252 - t302 * t253) / 0.2e1 + (((t230 * t265 + t232 * t263) * t253 + (t231 * t265 + t233 * t263) * t252) * qJD(4) + (t222 * t259 + t224 * t258) * t240 + (t223 * t259 + t225 * t258) * t239 + (t259 * t243 + t258 * t244 + t265 * t247 + t263 * t248) * t260) * t260 / 0.2e1 + (m(2) * (t250 ^ 2 + t251 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
