% Calculate kinetic energy for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR9_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:07
% EndTime: 2019-12-31 16:56:08
% DurationCPUTime: 0.58s
% Computational Cost: add. (296->140), mult. (702->238), div. (0->0), fcn. (670->6), ass. (0->75)
t231 = sin(qJ(3));
t257 = Icges(4,4) * t231;
t234 = cos(qJ(3));
t256 = Icges(4,4) * t234;
t230 = sin(qJ(4));
t232 = sin(qJ(1));
t255 = t230 * t232;
t235 = cos(qJ(1));
t254 = t230 * t235;
t233 = cos(qJ(4));
t253 = t232 * t233;
t252 = t232 * t234;
t251 = t233 * t235;
t250 = t234 * t235;
t218 = qJD(1) * (pkin(1) * t235 + qJ(2) * t232);
t249 = qJD(1) * t235 * pkin(5) + t218;
t248 = qJD(3) * t232;
t247 = qJD(3) * t235;
t246 = qJD(4) * t234;
t222 = pkin(1) * t232 - qJ(2) * t235;
t245 = -pkin(5) * t232 - t222;
t244 = pkin(3) * t231 - pkin(6) * t234;
t243 = rSges(4,1) * t231 + rSges(4,2) * t234;
t242 = Icges(4,1) * t231 + t256;
t241 = Icges(4,2) * t234 + t257;
t240 = Icges(4,5) * t231 + Icges(4,6) * t234;
t202 = Icges(4,6) * t235 + t241 * t232;
t205 = Icges(4,5) * t235 + t242 * t232;
t239 = -t202 * t234 - t205 * t231;
t203 = Icges(4,6) * t232 - t241 * t235;
t206 = Icges(4,5) * t232 - t242 * t235;
t238 = t203 * t234 + t206 * t231;
t220 = -Icges(4,2) * t231 + t256;
t221 = Icges(4,1) * t234 - t257;
t237 = t220 * t234 + t221 * t231;
t229 = qJD(2) * t232;
t227 = qJD(4) * t231 + qJD(1);
t226 = pkin(3) * t234 + pkin(6) * t231;
t225 = rSges(2,1) * t235 - rSges(2,2) * t232;
t224 = rSges(4,1) * t234 - rSges(4,2) * t231;
t223 = rSges(2,1) * t232 + rSges(2,2) * t235;
t219 = Icges(4,5) * t234 - Icges(4,6) * t231;
t217 = -t232 * t246 + t247;
t216 = t235 * t246 + t248;
t215 = t244 * t235;
t214 = t244 * t232;
t213 = -t231 * t251 + t255;
t212 = t231 * t254 + t253;
t211 = t231 * t253 + t254;
t210 = -t231 * t255 + t251;
t209 = rSges(4,3) * t232 - t243 * t235;
t208 = rSges(5,3) * t231 + (rSges(5,1) * t233 - rSges(5,2) * t230) * t234;
t207 = rSges(4,3) * t235 + t243 * t232;
t204 = Icges(5,5) * t231 + (Icges(5,1) * t233 - Icges(5,4) * t230) * t234;
t201 = Icges(5,6) * t231 + (Icges(5,4) * t233 - Icges(5,2) * t230) * t234;
t200 = Icges(4,3) * t232 - t240 * t235;
t199 = Icges(4,3) * t235 + t240 * t232;
t198 = Icges(5,3) * t231 + (Icges(5,5) * t233 - Icges(5,6) * t230) * t234;
t197 = t218 - qJD(2) * t235 + qJD(1) * (-rSges(3,2) * t235 + rSges(3,3) * t232);
t196 = t229 + (rSges(3,2) * t232 + rSges(3,3) * t235 - t222) * qJD(1);
t195 = rSges(5,1) * t213 + rSges(5,2) * t212 + rSges(5,3) * t250;
t194 = rSges(5,1) * t211 + rSges(5,2) * t210 - rSges(5,3) * t252;
t193 = Icges(5,1) * t213 + Icges(5,4) * t212 + Icges(5,5) * t250;
t192 = Icges(5,1) * t211 + Icges(5,4) * t210 - Icges(5,5) * t252;
t191 = Icges(5,4) * t213 + Icges(5,2) * t212 + Icges(5,6) * t250;
t190 = Icges(5,4) * t211 + Icges(5,2) * t210 - Icges(5,6) * t252;
t189 = Icges(5,5) * t213 + Icges(5,6) * t212 + Icges(5,3) * t250;
t188 = Icges(5,5) * t211 + Icges(5,6) * t210 - Icges(5,3) * t252;
t187 = (-t207 * t232 + t209 * t235) * qJD(3);
t186 = qJD(1) * t207 + (-qJD(3) * t224 - qJD(2)) * t235 + t249;
t185 = t224 * t248 + t229 + (-t209 + t245) * qJD(1);
t184 = qJD(1) * t214 + t194 * t227 - t208 * t217 + (-qJD(3) * t226 - qJD(2)) * t235 + t249;
t183 = t226 * t248 - t195 * t227 + t208 * t216 + t229 + (t215 + t245) * qJD(1);
t182 = -t194 * t216 + t195 * t217 + (-t214 * t232 - t215 * t235) * qJD(3);
t1 = m(3) * (t196 ^ 2 + t197 ^ 2) / 0.2e1 + m(4) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + ((t235 * t219 + t237 * t232) * qJD(1) + (t235 ^ 2 * t199 + (t238 * t232 + (t200 - t239) * t235) * t232) * qJD(3)) * t247 / 0.2e1 + ((t232 * t219 - t237 * t235) * qJD(1) + (t232 ^ 2 * t200 + (t239 * t235 + (t199 - t238) * t232) * t235) * qJD(3)) * t248 / 0.2e1 + qJD(1) * ((-t231 * t220 + t234 * t221) * qJD(1) + ((-t202 * t231 + t205 * t234) * t235 + (-t203 * t231 + t206 * t234) * t232) * qJD(3)) / 0.2e1 + m(5) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + t217 * ((-t188 * t252 + t210 * t190 + t211 * t192) * t217 + (-t189 * t252 + t191 * t210 + t193 * t211) * t216 + (-t198 * t252 + t201 * t210 + t204 * t211) * t227) / 0.2e1 + t216 * ((t188 * t250 + t190 * t212 + t192 * t213) * t217 + (t189 * t250 + t212 * t191 + t213 * t193) * t216 + (t198 * t250 + t201 * t212 + t204 * t213) * t227) / 0.2e1 + t227 * ((t188 * t217 + t189 * t216 + t198 * t227) * t231 + ((-t190 * t230 + t192 * t233) * t217 + (-t191 * t230 + t193 * t233) * t216 + (-t201 * t230 + t204 * t233) * t227) * t234) / 0.2e1 + (m(2) * (t223 ^ 2 + t225 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
