% Calculate kinetic energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR11_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:31
% EndTime: 2019-12-31 18:05:32
% DurationCPUTime: 0.60s
% Computational Cost: add. (319->152), mult. (741->250), div. (0->0), fcn. (692->6), ass. (0->79)
t263 = sin(qJ(1));
t291 = pkin(6) * t263;
t262 = sin(qJ(4));
t290 = Icges(5,4) * t262;
t265 = cos(qJ(4));
t289 = Icges(5,4) * t265;
t261 = sin(qJ(5));
t288 = t261 * t263;
t266 = cos(qJ(1));
t287 = t261 * t266;
t264 = cos(qJ(5));
t286 = t263 * t264;
t285 = t263 * t265;
t284 = t264 * t266;
t283 = t265 * t266;
t260 = qJD(2) * t263;
t282 = qJD(3) * t266 + t260;
t281 = qJD(4) * t263;
t280 = qJD(4) * t266;
t279 = qJD(5) * t265;
t278 = -qJD(2) * t266 + qJD(1) * (pkin(1) * t266 + qJ(2) * t263);
t277 = pkin(4) * t262 - pkin(7) * t265;
t276 = rSges(5,1) * t262 + rSges(5,2) * t265;
t275 = Icges(5,1) * t262 + t289;
t274 = Icges(5,2) * t265 + t290;
t273 = Icges(5,5) * t262 + Icges(5,6) * t265;
t231 = Icges(5,6) * t266 + t263 * t274;
t234 = Icges(5,5) * t266 + t263 * t275;
t272 = t231 * t265 + t234 * t262;
t232 = -Icges(5,6) * t263 + t266 * t274;
t235 = -Icges(5,5) * t263 + t266 * t275;
t271 = -t232 * t265 - t235 * t262;
t249 = -Icges(5,2) * t262 + t289;
t250 = Icges(5,1) * t265 - t290;
t270 = t249 * t265 + t250 * t262;
t269 = qJD(1) * t266 * qJ(3) + qJD(3) * t263 + t278;
t251 = pkin(1) * t263 - qJ(2) * t266;
t268 = -pkin(6) * t266 - qJ(3) * t263 - t251;
t256 = qJD(5) * t262 + qJD(1);
t255 = pkin(4) * t265 + pkin(7) * t262;
t254 = rSges(2,1) * t266 - rSges(2,2) * t263;
t253 = rSges(5,1) * t265 - rSges(5,2) * t262;
t252 = rSges(2,1) * t263 + rSges(2,2) * t266;
t248 = Icges(5,5) * t265 - Icges(5,6) * t262;
t246 = -t263 * t279 + t280;
t245 = -t266 * t279 - t281;
t244 = t277 * t266;
t243 = t277 * t263;
t242 = t262 * t284 - t288;
t241 = -t262 * t287 - t286;
t240 = t262 * t286 + t287;
t239 = -t262 * t288 + t284;
t238 = -rSges(5,3) * t263 + t266 * t276;
t237 = rSges(6,3) * t262 + (rSges(6,1) * t264 - rSges(6,2) * t261) * t265;
t236 = rSges(5,3) * t266 + t263 * t276;
t233 = Icges(6,5) * t262 + (Icges(6,1) * t264 - Icges(6,4) * t261) * t265;
t230 = Icges(6,6) * t262 + (Icges(6,4) * t264 - Icges(6,2) * t261) * t265;
t229 = -Icges(5,3) * t263 + t266 * t273;
t228 = Icges(5,3) * t266 + t263 * t273;
t227 = Icges(6,3) * t262 + (Icges(6,5) * t264 - Icges(6,6) * t261) * t265;
t226 = qJD(1) * (-rSges(3,2) * t266 + rSges(3,3) * t263) + t278;
t225 = t260 + (rSges(3,2) * t263 + rSges(3,3) * t266 - t251) * qJD(1);
t224 = qJD(1) * (rSges(4,2) * t263 + rSges(4,3) * t266) + t269;
t223 = (t266 * rSges(4,2) - t251 + (-rSges(4,3) - qJ(3)) * t263) * qJD(1) + t282;
t222 = rSges(6,1) * t242 + rSges(6,2) * t241 - rSges(6,3) * t283;
t221 = rSges(6,1) * t240 + rSges(6,2) * t239 - rSges(6,3) * t285;
t220 = Icges(6,1) * t242 + Icges(6,4) * t241 - Icges(6,5) * t283;
t219 = Icges(6,1) * t240 + Icges(6,4) * t239 - Icges(6,5) * t285;
t218 = Icges(6,4) * t242 + Icges(6,2) * t241 - Icges(6,6) * t283;
t217 = Icges(6,4) * t240 + Icges(6,2) * t239 - Icges(6,6) * t285;
t216 = Icges(6,5) * t242 + Icges(6,6) * t241 - Icges(6,3) * t283;
t215 = Icges(6,5) * t240 + Icges(6,6) * t239 - Icges(6,3) * t285;
t214 = (-t236 * t263 - t238 * t266) * qJD(4);
t213 = t253 * t281 + (t238 - t291) * qJD(1) + t269;
t212 = t253 * t280 + (-t236 + t268) * qJD(1) + t282;
t211 = t255 * t281 + t222 * t256 - t237 * t245 + (t244 - t291) * qJD(1) + t269;
t210 = t255 * t280 - t221 * t256 + t237 * t246 + (-t243 + t268) * qJD(1) + t282;
t209 = t221 * t245 - t222 * t246 + (-t243 * t263 - t244 * t266) * qJD(4);
t1 = m(3) * (t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(4) * (t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(5) * (t212 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 - ((-t263 * t248 + t266 * t270) * qJD(1) + (t263 ^ 2 * t229 + (t272 * t266 + (-t228 + t271) * t263) * t266) * qJD(4)) * t281 / 0.2e1 + ((t266 * t248 + t263 * t270) * qJD(1) + (t266 ^ 2 * t228 + (t271 * t263 + (-t229 + t272) * t266) * t263) * qJD(4)) * t280 / 0.2e1 + qJD(1) * ((-t262 * t249 + t265 * t250) * qJD(1) + (-(-t232 * t262 + t235 * t265) * t263 + (-t231 * t262 + t234 * t265) * t266) * qJD(4)) / 0.2e1 + m(6) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + t245 * ((-t216 * t283 + t241 * t218 + t242 * t220) * t245 + (-t215 * t283 + t217 * t241 + t219 * t242) * t246 + (-t227 * t283 + t230 * t241 + t233 * t242) * t256) / 0.2e1 + t246 * ((-t216 * t285 + t218 * t239 + t220 * t240) * t245 + (-t215 * t285 + t239 * t217 + t240 * t219) * t246 + (-t227 * t285 + t230 * t239 + t233 * t240) * t256) / 0.2e1 + t256 * ((t215 * t246 + t216 * t245 + t227 * t256) * t262 + ((-t218 * t261 + t220 * t264) * t245 + (-t217 * t261 + t219 * t264) * t246 + (-t230 * t261 + t233 * t264) * t256) * t265) / 0.2e1 + (m(2) * (t252 ^ 2 + t254 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
