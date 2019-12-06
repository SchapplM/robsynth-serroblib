% Calculate kinetic energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:36
% EndTime: 2019-12-05 18:08:37
% DurationCPUTime: 1.06s
% Computational Cost: add. (593->192), mult. (1436->326), div. (0->0), fcn. (1566->8), ass. (0->95)
t262 = sin(qJ(3));
t288 = Icges(4,4) * t262;
t266 = cos(qJ(3));
t287 = Icges(4,4) * t266;
t261 = sin(qJ(4));
t286 = t261 * t262;
t263 = sin(qJ(1));
t285 = t262 * t263;
t265 = cos(qJ(4));
t284 = t262 * t265;
t267 = cos(qJ(1));
t283 = t262 * t267;
t282 = t263 * t266;
t281 = t266 * t267;
t279 = qJ(2) * qJD(1);
t280 = qJD(2) * t263 + t267 * t279;
t258 = qJD(3) * t263;
t277 = qJD(4) * t262;
t244 = t267 * t277 + t258;
t278 = qJD(3) * t267;
t276 = -qJD(2) * t267 + t263 * t279;
t245 = t263 * t277 - t278;
t254 = -qJD(4) * t266 + qJD(1);
t275 = rSges(4,1) * t266 - rSges(4,2) * t262;
t274 = Icges(4,1) * t266 - t288;
t273 = -Icges(4,2) * t262 + t287;
t272 = Icges(4,5) * t266 - Icges(4,6) * t262;
t229 = -Icges(4,6) * t267 + t263 * t273;
t232 = -Icges(4,5) * t267 + t263 * t274;
t271 = t229 * t262 - t232 * t266;
t230 = Icges(4,6) * t263 + t267 * t273;
t233 = Icges(4,5) * t263 + t267 * t274;
t270 = -t230 * t262 + t233 * t266;
t247 = Icges(4,2) * t266 + t288;
t248 = Icges(4,1) * t262 + t287;
t269 = -t247 * t262 + t248 * t266;
t264 = cos(qJ(5));
t260 = sin(qJ(5));
t251 = rSges(2,1) * t267 - rSges(2,2) * t263;
t250 = rSges(2,1) * t263 + rSges(2,2) * t267;
t249 = rSges(4,1) * t262 + rSges(4,2) * t266;
t246 = Icges(4,5) * t262 + Icges(4,6) * t266;
t243 = qJD(5) * t286 + t254;
t242 = t261 * t263 + t265 * t281;
t241 = t261 * t281 - t263 * t265;
t240 = -t261 * t267 + t265 * t282;
t239 = t261 * t282 + t265 * t267;
t238 = -t260 * t266 + t264 * t284;
t237 = -t260 * t284 - t264 * t266;
t236 = rSges(4,3) * t263 + t267 * t275;
t235 = -rSges(4,3) * t267 + t263 * t275;
t234 = -rSges(5,3) * t266 + (rSges(5,1) * t265 - rSges(5,2) * t261) * t262;
t231 = -Icges(5,5) * t266 + (Icges(5,1) * t265 - Icges(5,4) * t261) * t262;
t228 = -Icges(5,6) * t266 + (Icges(5,4) * t265 - Icges(5,2) * t261) * t262;
t227 = Icges(4,3) * t263 + t267 * t272;
t226 = -Icges(4,3) * t267 + t263 * t272;
t225 = -Icges(5,3) * t266 + (Icges(5,5) * t265 - Icges(5,6) * t261) * t262;
t224 = -qJD(1) * (rSges(3,1) * t263 - rSges(3,3) * t267) + t280;
t223 = qJD(1) * (rSges(3,1) * t267 + rSges(3,3) * t263) + t276;
t222 = t242 * t264 + t260 * t283;
t221 = -t242 * t260 + t264 * t283;
t220 = t240 * t264 + t260 * t285;
t219 = -t240 * t260 + t264 * t285;
t218 = qJD(5) * t239 + t245;
t217 = qJD(5) * t241 + t244;
t216 = rSges(5,1) * t242 - rSges(5,2) * t241 + rSges(5,3) * t283;
t215 = rSges(5,1) * t240 - rSges(5,2) * t239 + rSges(5,3) * t285;
t214 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t286;
t213 = Icges(5,1) * t242 - Icges(5,4) * t241 + Icges(5,5) * t283;
t212 = Icges(5,1) * t240 - Icges(5,4) * t239 + Icges(5,5) * t285;
t211 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t286;
t210 = Icges(5,4) * t242 - Icges(5,2) * t241 + Icges(5,6) * t283;
t209 = Icges(5,4) * t240 - Icges(5,2) * t239 + Icges(5,6) * t285;
t208 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t286;
t207 = Icges(5,5) * t242 - Icges(5,6) * t241 + Icges(5,3) * t283;
t206 = Icges(5,5) * t240 - Icges(5,6) * t239 + Icges(5,3) * t285;
t205 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t286;
t204 = -qJD(1) * t235 - t249 * t278 + t280;
t203 = qJD(1) * t236 - t249 * t258 + t276;
t202 = (t235 * t263 + t236 * t267) * qJD(3);
t201 = rSges(6,1) * t222 + rSges(6,2) * t221 + rSges(6,3) * t241;
t200 = rSges(6,1) * t220 + rSges(6,2) * t219 + rSges(6,3) * t239;
t199 = Icges(6,1) * t222 + Icges(6,4) * t221 + Icges(6,5) * t241;
t198 = Icges(6,1) * t220 + Icges(6,4) * t219 + Icges(6,5) * t239;
t197 = Icges(6,4) * t222 + Icges(6,2) * t221 + Icges(6,6) * t241;
t196 = Icges(6,4) * t220 + Icges(6,2) * t219 + Icges(6,6) * t239;
t195 = Icges(6,5) * t222 + Icges(6,6) * t221 + Icges(6,3) * t241;
t194 = Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t239;
t193 = -t215 * t254 + t234 * t245 + t280;
t192 = t216 * t254 - t234 * t244 + t276;
t191 = t215 * t244 - t216 * t245;
t190 = -t200 * t243 + t214 * t218 + t280;
t189 = t201 * t243 - t214 * t217 + t276;
t188 = t200 * t217 - t201 * t218;
t1 = m(3) * (t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + ((t246 * t263 + t267 * t269) * qJD(1) + (t227 * t263 ^ 2 + (t271 * t267 + (-t226 + t270) * t263) * t267) * qJD(3)) * t258 / 0.2e1 - ((-t246 * t267 + t263 * t269) * qJD(1) + (t226 * t267 ^ 2 + (t270 * t263 + (-t227 + t271) * t267) * t263) * qJD(3)) * t278 / 0.2e1 + qJD(1) * ((t247 * t266 + t248 * t262) * qJD(1) + ((t230 * t266 + t233 * t262) * t263 - (t229 * t266 + t232 * t262) * t267) * qJD(3)) / 0.2e1 + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + t244 * ((t207 * t283 - t210 * t241 + t213 * t242) * t244 + (t206 * t283 - t209 * t241 + t212 * t242) * t245 + (t225 * t283 - t228 * t241 + t231 * t242) * t254) / 0.2e1 + t245 * ((t207 * t285 - t210 * t239 + t213 * t240) * t244 + (t206 * t285 - t209 * t239 + t212 * t240) * t245 + (t225 * t285 - t228 * t239 + t231 * t240) * t254) / 0.2e1 + t254 * ((-t206 * t245 - t207 * t244 - t225 * t254) * t266 + ((-t210 * t261 + t213 * t265) * t244 + (-t209 * t261 + t212 * t265) * t245 + (-t228 * t261 + t231 * t265) * t254) * t262) / 0.2e1 + m(6) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + t217 * ((t195 * t241 + t197 * t221 + t199 * t222) * t217 + (t194 * t241 + t196 * t221 + t198 * t222) * t218 + (t205 * t241 + t208 * t221 + t211 * t222) * t243) / 0.2e1 + t218 * ((t195 * t239 + t197 * t219 + t199 * t220) * t217 + (t194 * t239 + t196 * t219 + t198 * t220) * t218 + (t205 * t239 + t208 * t219 + t211 * t220) * t243) / 0.2e1 + t243 * ((t195 * t286 + t197 * t237 + t199 * t238) * t217 + (t194 * t286 + t196 * t237 + t198 * t238) * t218 + (t205 * t286 + t208 * t237 + t211 * t238) * t243) / 0.2e1 + (m(2) * (t250 ^ 2 + t251 ^ 2) + Icges(2,3) + Icges(3,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
