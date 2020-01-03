% Calculate kinetic energy for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:34
% EndTime: 2019-12-31 16:33:34
% DurationCPUTime: 0.66s
% Computational Cost: add. (516->120), mult. (753->218), div. (0->0), fcn. (730->8), ass. (0->68)
t246 = sin(pkin(7));
t247 = cos(pkin(7));
t285 = t246 * t247;
t282 = t247 ^ 2;
t283 = t246 ^ 2;
t284 = t282 + t283;
t281 = qJD(2) ^ 2;
t251 = cos(qJ(2));
t280 = pkin(2) * t251;
t245 = qJ(2) + qJ(3);
t243 = sin(t245);
t278 = t243 * t246;
t277 = t243 * t247;
t248 = sin(qJ(4));
t276 = t246 * t248;
t250 = cos(qJ(4));
t275 = t246 * t250;
t274 = t247 * t248;
t273 = t247 * t250;
t242 = qJD(2) * t246;
t236 = qJD(3) * t246 + t242;
t272 = qJD(4) * t243;
t244 = cos(t245);
t271 = qJD(4) * t244;
t249 = sin(qJ(2));
t270 = pkin(2) * qJD(2) * t249;
t269 = qJD(1) + (-pkin(5) * t247 + t280 * t246) * t242 + qJD(2) * t247 * (pkin(5) * t246 + t280 * t247);
t237 = (-qJD(2) - qJD(3)) * t247;
t268 = t246 * t270;
t267 = t247 * t270;
t265 = rSges(4,1) * t244 - rSges(4,2) * t243;
t263 = Icges(4,1) * t244 - Icges(4,4) * t243;
t261 = Icges(4,4) * t244 - Icges(4,2) * t243;
t260 = Icges(3,5) * t251 - Icges(3,6) * t249;
t259 = Icges(4,5) * t244 - Icges(4,6) * t243;
t258 = (-Icges(4,3) * t247 + t259 * t246) * t237 + (Icges(4,3) * t246 + t259 * t247) * t236;
t254 = (-(Icges(4,6) * t246 + t261 * t247) * t243 + (Icges(4,5) * t246 + t263 * t247) * t244) * t236 + (-(-Icges(4,6) * t247 + t261 * t246) * t243 + (-Icges(4,5) * t247 + t263 * t246) * t244) * t237;
t239 = rSges(3,1) * t249 + rSges(3,2) * t251;
t235 = pkin(3) * t243 - pkin(6) * t244;
t234 = rSges(4,1) * t243 + rSges(4,2) * t244;
t233 = t244 * t273 + t276;
t232 = -t244 * t274 + t275;
t231 = t244 * t275 - t274;
t230 = -t244 * t276 - t273;
t229 = t246 * t272 + t237;
t228 = t247 * t272 + t236;
t223 = Icges(3,3) * t246 + t260 * t247;
t222 = -Icges(3,3) * t247 + t260 * t246;
t215 = -rSges(5,3) * t244 + (rSges(5,1) * t250 - rSges(5,2) * t248) * t243;
t214 = -Icges(5,5) * t244 + (Icges(5,1) * t250 - Icges(5,4) * t248) * t243;
t213 = -Icges(5,6) * t244 + (Icges(5,4) * t250 - Icges(5,2) * t248) * t243;
t212 = -Icges(5,3) * t244 + (Icges(5,5) * t250 - Icges(5,6) * t248) * t243;
t209 = t234 * t237 - t267;
t208 = -t234 * t236 - t268;
t207 = rSges(5,1) * t233 + rSges(5,2) * t232 + rSges(5,3) * t277;
t206 = rSges(5,1) * t231 + rSges(5,2) * t230 + rSges(5,3) * t278;
t205 = Icges(5,1) * t233 + Icges(5,4) * t232 + Icges(5,5) * t277;
t204 = Icges(5,1) * t231 + Icges(5,4) * t230 + Icges(5,5) * t278;
t203 = Icges(5,4) * t233 + Icges(5,2) * t232 + Icges(5,6) * t277;
t202 = Icges(5,4) * t231 + Icges(5,2) * t230 + Icges(5,6) * t278;
t201 = Icges(5,5) * t233 + Icges(5,6) * t232 + Icges(5,3) * t277;
t200 = Icges(5,5) * t231 + Icges(5,6) * t230 + Icges(5,3) * t278;
t199 = qJD(1) + t284 * qJD(2) * (rSges(3,1) * t251 - rSges(3,2) * t249);
t198 = t236 * (-rSges(4,3) * t247 + t265 * t246) - t237 * (rSges(4,3) * t246 + t265 * t247) + t269;
t197 = t206 * t271 + t215 * t229 + t235 * t237 - t267;
t196 = -t207 * t271 - t215 * t228 - t235 * t236 - t268;
t195 = t228 * t206 - t229 * t207 + t269 + (t246 * t236 - t247 * t237) * (pkin(3) * t244 + pkin(6) * t243);
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t284 * t281 * t239 ^ 2 + t199 ^ 2) / 0.2e1 + t281 * t246 * (-t222 * t285 + t283 * t223) / 0.2e1 - t281 * t247 * (t282 * t222 - t223 * t285) / 0.2e1 + m(4) * (t198 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + t236 * (t258 * t246 + t254 * t247) / 0.2e1 + t237 * (t254 * t246 - t258 * t247) / 0.2e1 + m(5) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + t228 * ((t201 * t277 + t232 * t203 + t233 * t205) * t228 + (t200 * t277 + t202 * t232 + t204 * t233) * t229 - (t212 * t277 + t213 * t232 + t214 * t233) * t271) / 0.2e1 + t229 * ((t201 * t278 + t203 * t230 + t205 * t231) * t228 + (t200 * t278 + t230 * t202 + t231 * t204) * t229 - (t212 * t278 + t213 * t230 + t214 * t231) * t271) / 0.2e1 - ((-t200 * t229 - t201 * t228 + t212 * t271) * t244 + ((-t203 * t248 + t205 * t250) * t228 + (-t202 * t248 + t204 * t250) * t229 - (-t213 * t248 + t214 * t250) * t271) * t243) * t271 / 0.2e1;
T = t1;
