% Calculate kinetic energy for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR4_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:13
% EndTime: 2019-12-31 16:50:14
% DurationCPUTime: 0.57s
% Computational Cost: add. (569->138), mult. (695->236), div. (0->0), fcn. (662->8), ass. (0->75)
t223 = sin(qJ(1));
t247 = pkin(1) * t223;
t222 = sin(qJ(3));
t246 = Icges(4,4) * t222;
t225 = cos(qJ(3));
t245 = Icges(4,4) * t225;
t220 = qJ(1) + pkin(7);
t218 = sin(t220);
t244 = t218 * t222;
t219 = cos(t220);
t243 = t219 * t222;
t221 = sin(qJ(4));
t242 = t221 * t225;
t224 = cos(qJ(4));
t241 = t224 * t225;
t226 = cos(qJ(1));
t217 = qJD(1) * t226 * pkin(1);
t240 = qJD(1) * (pkin(2) * t219 + pkin(5) * t218) + t217;
t239 = qJD(3) * t218;
t238 = qJD(3) * t219;
t237 = qJD(4) * t222;
t236 = -pkin(2) * t218 + pkin(5) * t219 - t247;
t235 = pkin(3) * t225 + pkin(6) * t222;
t234 = rSges(4,1) * t225 - rSges(4,2) * t222;
t233 = Icges(4,1) * t225 - t246;
t232 = -Icges(4,2) * t222 + t245;
t231 = Icges(4,5) * t225 - Icges(4,6) * t222;
t187 = -Icges(4,6) * t219 + t232 * t218;
t189 = -Icges(4,5) * t219 + t233 * t218;
t230 = t187 * t222 - t189 * t225;
t188 = Icges(4,6) * t218 + t232 * t219;
t190 = Icges(4,5) * t218 + t233 * t219;
t229 = -t188 * t222 + t190 * t225;
t210 = Icges(4,2) * t225 + t246;
t211 = Icges(4,1) * t222 + t245;
t228 = -t210 * t222 + t211 * t225;
t216 = -qJD(4) * t225 + qJD(1);
t215 = pkin(3) * t222 - pkin(6) * t225;
t214 = rSges(2,1) * t226 - rSges(2,2) * t223;
t213 = rSges(2,1) * t223 + rSges(2,2) * t226;
t212 = rSges(4,1) * t222 + rSges(4,2) * t225;
t209 = Icges(4,5) * t222 + Icges(4,6) * t225;
t206 = t218 * t237 - t238;
t205 = t219 * t237 + t239;
t204 = t235 * t219;
t203 = t235 * t218;
t202 = -rSges(5,3) * t225 + (rSges(5,1) * t224 - rSges(5,2) * t221) * t222;
t201 = -Icges(5,5) * t225 + (Icges(5,1) * t224 - Icges(5,4) * t221) * t222;
t200 = -Icges(5,6) * t225 + (Icges(5,4) * t224 - Icges(5,2) * t221) * t222;
t199 = -Icges(5,3) * t225 + (Icges(5,5) * t224 - Icges(5,6) * t221) * t222;
t198 = t218 * t221 + t219 * t241;
t197 = t218 * t224 - t219 * t242;
t196 = t218 * t241 - t219 * t221;
t195 = -t218 * t242 - t219 * t224;
t194 = t217 + qJD(1) * (rSges(3,1) * t219 - rSges(3,2) * t218);
t193 = (-rSges(3,1) * t218 - rSges(3,2) * t219 - t247) * qJD(1);
t192 = rSges(4,3) * t218 + t234 * t219;
t191 = -rSges(4,3) * t219 + t234 * t218;
t186 = Icges(4,3) * t218 + t231 * t219;
t185 = -Icges(4,3) * t219 + t231 * t218;
t184 = rSges(5,1) * t198 + rSges(5,2) * t197 + rSges(5,3) * t243;
t183 = rSges(5,1) * t196 + rSges(5,2) * t195 + rSges(5,3) * t244;
t182 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t243;
t181 = Icges(5,1) * t196 + Icges(5,4) * t195 + Icges(5,5) * t244;
t180 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t243;
t179 = Icges(5,4) * t196 + Icges(5,2) * t195 + Icges(5,6) * t244;
t178 = Icges(5,5) * t198 + Icges(5,6) * t197 + Icges(5,3) * t243;
t177 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t244;
t176 = qJD(1) * t192 - t212 * t239 + t240;
t175 = -t212 * t238 + (-t191 + t236) * qJD(1);
t174 = qJD(2) + (t191 * t218 + t192 * t219) * qJD(3);
t173 = qJD(1) * t204 + t184 * t216 - t202 * t205 - t215 * t239 + t240;
t172 = -t215 * t238 - t183 * t216 + t202 * t206 + (-t203 + t236) * qJD(1);
t171 = t183 * t205 - t184 * t206 + qJD(2) + (t203 * t218 + t204 * t219) * qJD(3);
t1 = m(3) * (qJD(2) ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(4) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + ((t218 * t209 + t228 * t219) * qJD(1) + (t218 ^ 2 * t186 + (t230 * t219 + (-t185 + t229) * t218) * t219) * qJD(3)) * t239 / 0.2e1 - ((-t219 * t209 + t228 * t218) * qJD(1) + (t219 ^ 2 * t185 + (t229 * t218 + (-t186 + t230) * t219) * t218) * qJD(3)) * t238 / 0.2e1 + qJD(1) * ((t225 * t210 + t222 * t211) * qJD(1) + ((t225 * t188 + t190 * t222) * t218 - (t187 * t225 + t189 * t222) * t219) * qJD(3)) / 0.2e1 + m(5) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + t205 * ((t178 * t243 + t180 * t197 + t182 * t198) * t205 + (t177 * t243 + t179 * t197 + t181 * t198) * t206 + (t197 * t200 + t198 * t201 + t199 * t243) * t216) / 0.2e1 + t206 * ((t178 * t244 + t180 * t195 + t182 * t196) * t205 + (t177 * t244 + t179 * t195 + t181 * t196) * t206 + (t195 * t200 + t196 * t201 + t199 * t244) * t216) / 0.2e1 + t216 * ((-t177 * t206 - t178 * t205 - t199 * t216) * t225 + ((-t180 * t221 + t182 * t224) * t205 + (-t179 * t221 + t181 * t224) * t206 + (-t200 * t221 + t201 * t224) * t216) * t222) / 0.2e1 + (m(2) * (t213 ^ 2 + t214 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
