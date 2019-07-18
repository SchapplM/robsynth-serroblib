% Calculate kinetic energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:19
% EndTime: 2019-07-18 13:28:19
% DurationCPUTime: 0.94s
% Computational Cost: add. (600->166), mult. (887->290), div. (0->0), fcn. (836->8), ass. (0->96)
t246 = sin(qJ(2));
t241 = t246 ^ 2;
t249 = cos(qJ(2));
t242 = t249 ^ 2;
t245 = sin(qJ(3));
t279 = Icges(4,4) * t245;
t248 = cos(qJ(3));
t278 = Icges(4,4) * t248;
t243 = qJ(3) + qJ(4);
t239 = sin(t243);
t277 = Icges(5,4) * t239;
t240 = cos(t243);
t276 = Icges(5,4) * t240;
t275 = t239 * t246;
t274 = t239 * t249;
t244 = sin(qJ(5));
t273 = t244 * t246;
t272 = t244 * t249;
t247 = cos(qJ(5));
t271 = t246 * t247;
t270 = t247 * t249;
t238 = qJD(3) * t246;
t228 = qJD(4) * t246 + t238;
t269 = qJD(2) * t248;
t268 = qJD(3) * t249;
t267 = qJD(5) * t239;
t229 = (-qJD(3) - qJD(4)) * t249;
t266 = rSges(4,1) * t248 - rSges(4,2) * t245;
t265 = rSges(5,1) * t240 - rSges(5,2) * t239;
t264 = Icges(4,1) * t248 - t279;
t263 = Icges(5,1) * t240 - t277;
t262 = -Icges(4,2) * t245 + t278;
t261 = -Icges(5,2) * t239 + t276;
t260 = Icges(4,5) * t248 - Icges(4,6) * t245;
t259 = Icges(5,5) * t240 - Icges(5,6) * t239;
t211 = -Icges(4,6) * t249 + t246 * t262;
t213 = -Icges(4,5) * t249 + t246 * t264;
t258 = t211 * t245 - t213 * t248;
t212 = Icges(4,6) * t246 + t249 * t262;
t214 = Icges(4,5) * t246 + t249 * t264;
t257 = -t212 * t245 + t214 * t248;
t231 = -Icges(4,2) * t248 - t279;
t232 = -Icges(4,1) * t245 - t278;
t256 = -t231 * t245 + t232 * t248;
t255 = (-t241 - t242) * t248 * qJD(3) * pkin(2);
t254 = qJD(1) + (-t238 * t245 + t249 * t269) * pkin(2);
t253 = qJD(2) * (-Icges(5,5) * t239 - Icges(5,6) * t240) - (-Icges(5,3) * t249 + t246 * t259) * t229 - (Icges(5,3) * t246 + t249 * t259) * t228;
t252 = (-t245 * t268 - t246 * t269) * pkin(2);
t203 = -Icges(5,6) * t249 + t246 * t261;
t204 = Icges(5,6) * t246 + t249 * t261;
t205 = -Icges(5,5) * t249 + t246 * t263;
t206 = Icges(5,5) * t246 + t249 * t263;
t224 = -Icges(5,2) * t240 - t277;
t225 = -Icges(5,1) * t239 - t276;
t251 = (-t204 * t239 + t206 * t240) * t228 - (-t224 * t239 + t225 * t240) * qJD(2) + (-t203 * t239 + t205 * t240) * t229;
t250 = qJD(2) ^ 2;
t235 = qJD(5) * t240 - qJD(2);
t234 = rSges(3,1) * t246 + rSges(3,2) * t249;
t233 = -rSges(4,1) * t245 - rSges(4,2) * t248;
t230 = -Icges(4,5) * t245 - Icges(4,6) * t248;
t227 = -rSges(5,1) * t239 - rSges(5,2) * t240;
t226 = qJD(1) + qJD(2) * (rSges(3,1) * t249 - rSges(3,2) * t246);
t222 = t240 * t270 + t273;
t221 = -t240 * t272 + t271;
t220 = t240 * t271 - t272;
t219 = -t240 * t273 - t270;
t218 = rSges(4,3) * t246 + t249 * t266;
t217 = -rSges(4,3) * t249 + t246 * t266;
t216 = t246 * t267 + t229;
t215 = t249 * t267 + t228;
t210 = Icges(4,3) * t246 + t249 * t260;
t209 = -Icges(4,3) * t249 + t246 * t260;
t208 = rSges(5,3) * t246 + t249 * t265;
t207 = -rSges(5,3) * t249 + t246 * t265;
t200 = rSges(6,3) * t240 + (-rSges(6,1) * t247 + rSges(6,2) * t244) * t239;
t199 = Icges(6,5) * t240 + (-Icges(6,1) * t247 + Icges(6,4) * t244) * t239;
t198 = Icges(6,6) * t240 + (-Icges(6,4) * t247 + Icges(6,2) * t244) * t239;
t197 = Icges(6,3) * t240 + (-Icges(6,5) * t247 + Icges(6,6) * t244) * t239;
t196 = -qJD(2) * t217 + t233 * t268;
t195 = qJD(2) * t218 + t233 * t238 + qJD(1);
t194 = rSges(6,1) * t222 + rSges(6,2) * t221 + rSges(6,3) * t274;
t193 = rSges(6,1) * t220 + rSges(6,2) * t219 + rSges(6,3) * t275;
t192 = Icges(6,1) * t222 + Icges(6,4) * t221 + Icges(6,5) * t274;
t191 = Icges(6,1) * t220 + Icges(6,4) * t219 + Icges(6,5) * t275;
t190 = Icges(6,4) * t222 + Icges(6,2) * t221 + Icges(6,6) * t274;
t189 = Icges(6,4) * t220 + Icges(6,2) * t219 + Icges(6,6) * t275;
t188 = Icges(6,5) * t222 + Icges(6,6) * t221 + Icges(6,3) * t274;
t187 = Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t275;
t186 = (-t217 * t246 - t218 * t249) * qJD(3);
t185 = -qJD(2) * t207 - t227 * t229 + t252;
t184 = qJD(2) * t208 + t227 * t228 + t254;
t183 = -t207 * t228 + t208 * t229 + t255;
t182 = t193 * t235 - t200 * t216 + t252;
t181 = -t194 * t235 + t200 * t215 + t254;
t180 = -t193 * t215 + t194 * t216 + t255;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t234 ^ 2 * t250 + t226 ^ 2) / 0.2e1 + t250 * Icges(3,3) / 0.2e1 + m(4) * (t186 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + (-(t246 * t230 + t249 * t256) * qJD(2) + (t241 * t210 + (t258 * t249 + (-t209 + t257) * t246) * t249) * qJD(3)) * t238 / 0.2e1 - (-(-t249 * t230 + t246 * t256) * qJD(2) + (t242 * t209 + (t257 * t246 + (-t210 + t258) * t249) * t246) * qJD(3)) * t268 / 0.2e1 + m(5) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + t228 * (-t253 * t246 + t251 * t249) / 0.2e1 + t229 * (t251 * t246 + t253 * t249) / 0.2e1 + m(6) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + t215 * ((t188 * t274 + t221 * t190 + t222 * t192) * t215 + (t197 * t274 + t198 * t221 + t199 * t222) * t235 + (t187 * t274 + t189 * t221 + t191 * t222) * t216) / 0.2e1 + t235 * ((t187 * t216 + t188 * t215 + t197 * t235) * t240 + ((t190 * t244 - t192 * t247) * t215 + (t198 * t244 - t199 * t247) * t235 + (t189 * t244 - t191 * t247) * t216) * t239) / 0.2e1 + t216 * ((t188 * t275 + t190 * t219 + t192 * t220) * t215 + (t197 * t275 + t198 * t219 + t199 * t220) * t235 + (t187 * t275 + t219 * t189 + t220 * t191) * t216) / 0.2e1 - (((-t212 * t248 - t214 * t245) * t246 - (-t211 * t248 - t213 * t245) * t249) * qJD(3) + (-t204 * t240 - t206 * t239) * t228 + (-t203 * t240 - t205 * t239) * t229 + (t240 * t224 + t239 * t225 + t248 * t231 + t245 * t232) * qJD(2)) * qJD(2) / 0.2e1;
T  = t1;
