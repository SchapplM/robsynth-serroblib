% Calculate kinetic energy for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR7_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:41
% EndTime: 2019-12-31 16:53:42
% DurationCPUTime: 0.65s
% Computational Cost: add. (555->145), mult. (730->248), div. (0->0), fcn. (698->8), ass. (0->80)
t246 = cos(pkin(7));
t276 = t246 * pkin(2);
t244 = pkin(7) + qJ(3);
t241 = sin(t244);
t275 = Icges(4,4) * t241;
t242 = cos(t244);
t274 = Icges(4,4) * t242;
t249 = sin(qJ(1));
t273 = t241 * t249;
t251 = cos(qJ(1));
t272 = t241 * t251;
t248 = sin(qJ(4));
t271 = t249 * t248;
t250 = cos(qJ(4));
t270 = t249 * t250;
t269 = t251 * t248;
t268 = t251 * t250;
t236 = t249 * pkin(1) - t251 * qJ(2);
t266 = pkin(5) * t251 - t276 * t249 - t236;
t265 = qJD(3) * t249;
t264 = qJD(3) * t251;
t263 = qJD(4) * t241;
t262 = pkin(3) * t242 + pkin(6) * t241;
t235 = qJD(1) * (t251 * pkin(1) + t249 * qJ(2));
t261 = -qJD(2) * t251 + qJD(1) * (pkin(5) * t249 + t276 * t251) + t235;
t245 = sin(pkin(7));
t260 = rSges(3,1) * t246 - rSges(3,2) * t245;
t259 = rSges(4,1) * t242 - rSges(4,2) * t241;
t258 = Icges(4,1) * t242 - t275;
t257 = -Icges(4,2) * t241 + t274;
t256 = Icges(4,5) * t242 - Icges(4,6) * t241;
t216 = -Icges(4,6) * t251 + t257 * t249;
t218 = -Icges(4,5) * t251 + t258 * t249;
t255 = t216 * t241 - t218 * t242;
t217 = Icges(4,6) * t249 + t257 * t251;
t219 = Icges(4,5) * t249 + t258 * t251;
t254 = -t217 * t241 + t219 * t242;
t231 = Icges(4,2) * t242 + t275;
t232 = Icges(4,1) * t241 + t274;
t253 = -t231 * t241 + t232 * t242;
t243 = qJD(2) * t249;
t239 = -qJD(4) * t242 + qJD(1);
t238 = t251 * rSges(2,1) - t249 * rSges(2,2);
t237 = t249 * rSges(2,1) + t251 * rSges(2,2);
t234 = t241 * pkin(3) - t242 * pkin(6);
t233 = t241 * rSges(4,1) + t242 * rSges(4,2);
t230 = Icges(4,5) * t241 + Icges(4,6) * t242;
t229 = t249 * t263 - t264;
t228 = t251 * t263 + t265;
t227 = t242 * t268 + t271;
t226 = -t242 * t269 + t270;
t225 = t242 * t270 - t269;
t224 = -t242 * t271 - t268;
t223 = t262 * t251;
t222 = t262 * t249;
t221 = t249 * rSges(4,3) + t259 * t251;
t220 = -t251 * rSges(4,3) + t259 * t249;
t215 = Icges(4,3) * t249 + t256 * t251;
t214 = -Icges(4,3) * t251 + t256 * t249;
t212 = -t242 * rSges(5,3) + (rSges(5,1) * t250 - rSges(5,2) * t248) * t241;
t211 = -Icges(5,5) * t242 + (Icges(5,1) * t250 - Icges(5,4) * t248) * t241;
t210 = -Icges(5,6) * t242 + (Icges(5,4) * t250 - Icges(5,2) * t248) * t241;
t209 = -Icges(5,3) * t242 + (Icges(5,5) * t250 - Icges(5,6) * t248) * t241;
t207 = qJD(1) * t249 * rSges(3,3) + t235 + (qJD(1) * t260 - qJD(2)) * t251;
t206 = t243 + (t251 * rSges(3,3) - t260 * t249 - t236) * qJD(1);
t205 = t227 * rSges(5,1) + t226 * rSges(5,2) + rSges(5,3) * t272;
t204 = t225 * rSges(5,1) + t224 * rSges(5,2) + rSges(5,3) * t273;
t203 = Icges(5,1) * t227 + Icges(5,4) * t226 + Icges(5,5) * t272;
t202 = Icges(5,1) * t225 + Icges(5,4) * t224 + Icges(5,5) * t273;
t201 = Icges(5,4) * t227 + Icges(5,2) * t226 + Icges(5,6) * t272;
t200 = Icges(5,4) * t225 + Icges(5,2) * t224 + Icges(5,6) * t273;
t199 = Icges(5,5) * t227 + Icges(5,6) * t226 + Icges(5,3) * t272;
t198 = Icges(5,5) * t225 + Icges(5,6) * t224 + Icges(5,3) * t273;
t197 = (t220 * t249 + t221 * t251) * qJD(3);
t196 = qJD(1) * t221 - t233 * t265 + t261;
t195 = -t233 * t264 + t243 + (-t220 + t266) * qJD(1);
t194 = t228 * t204 - t229 * t205 + (t222 * t249 + t223 * t251) * qJD(3);
t193 = qJD(1) * t223 + t239 * t205 - t228 * t212 - t234 * t265 + t261;
t192 = -t234 * t264 - t239 * t204 + t229 * t212 + t243 + (-t222 + t266) * qJD(1);
t1 = m(3) * (t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(4) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + ((t249 * t230 + t253 * t251) * qJD(1) + (t249 ^ 2 * t215 + (t255 * t251 + (-t214 + t254) * t249) * t251) * qJD(3)) * t265 / 0.2e1 - ((-t251 * t230 + t253 * t249) * qJD(1) + (t251 ^ 2 * t214 + (t254 * t249 + (-t215 + t255) * t251) * t249) * qJD(3)) * t264 / 0.2e1 + qJD(1) * ((t242 * t231 + t241 * t232) * qJD(1) + ((t242 * t217 + t241 * t219) * t249 - (t242 * t216 + t241 * t218) * t251) * qJD(3)) / 0.2e1 + m(5) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + t228 * ((t199 * t272 + t226 * t201 + t227 * t203) * t228 + (t198 * t272 + t226 * t200 + t227 * t202) * t229 + (t209 * t272 + t226 * t210 + t227 * t211) * t239) / 0.2e1 + t229 * ((t199 * t273 + t224 * t201 + t225 * t203) * t228 + (t198 * t273 + t224 * t200 + t225 * t202) * t229 + (t209 * t273 + t224 * t210 + t225 * t211) * t239) / 0.2e1 + t239 * ((-t198 * t229 - t199 * t228 - t209 * t239) * t242 + ((-t201 * t248 + t203 * t250) * t228 + (-t200 * t248 + t202 * t250) * t229 + (-t210 * t248 + t211 * t250) * t239) * t241) / 0.2e1 + (m(2) * (t237 ^ 2 + t238 ^ 2) + Icges(2,3) + Icges(3,2) * t246 ^ 2 + (Icges(3,1) * t245 + 0.2e1 * Icges(3,4) * t246) * t245) * qJD(1) ^ 2 / 0.2e1;
T = t1;
