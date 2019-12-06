% Calculate kinetic energy for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:42
% EndTime: 2019-12-05 16:50:44
% DurationCPUTime: 2.76s
% Computational Cost: add. (1291->278), mult. (1976->411), div. (0->0), fcn. (1936->8), ass. (0->135)
t279 = Icges(5,1) + Icges(6,1);
t278 = -Icges(5,4) + Icges(6,5);
t277 = Icges(6,4) + Icges(5,5);
t276 = Icges(5,2) + Icges(6,3);
t275 = -Icges(6,6) + Icges(5,6);
t274 = -Icges(5,3) - Icges(6,2);
t273 = rSges(6,1) + pkin(4);
t272 = rSges(6,3) + qJ(5);
t211 = qJ(3) + qJ(4);
t208 = sin(t211);
t209 = cos(t211);
t213 = cos(pkin(8));
t212 = sin(pkin(8));
t217 = cos(qJ(2));
t250 = t212 * t217;
t159 = t208 * t250 + t209 * t213;
t160 = -t208 * t213 + t209 * t250;
t215 = sin(qJ(2));
t251 = t212 * t215;
t269 = t276 * t159 + t278 * t160 - t275 * t251;
t247 = t213 * t217;
t161 = t208 * t247 - t212 * t209;
t162 = t208 * t212 + t209 * t247;
t248 = t213 * t215;
t268 = t276 * t161 + t278 * t162 - t275 * t248;
t267 = -t275 * t159 + t277 * t160 - t274 * t251;
t266 = -t275 * t161 + t277 * t162 - t274 * t248;
t265 = t278 * t159 + t279 * t160 + t277 * t251;
t264 = t278 * t161 + t279 * t162 + t277 * t248;
t263 = t275 * t217 + (t276 * t208 + t278 * t209) * t215;
t262 = t274 * t217 + (-t275 * t208 + t277 * t209) * t215;
t261 = -t277 * t217 + (t278 * t208 + t279 * t209) * t215;
t216 = cos(qJ(3));
t257 = pkin(3) * t216;
t255 = Icges(2,4) * t212;
t254 = Icges(3,4) * t215;
t253 = Icges(3,4) * t217;
t214 = sin(qJ(3));
t252 = t212 * t214;
t249 = t213 * t214;
t246 = t214 * t217;
t245 = t216 * t217;
t244 = rSges(6,2) * t251 + t272 * t159 + t273 * t160;
t243 = rSges(6,2) * t248 + t272 * t161 + t273 * t162;
t242 = -rSges(6,2) * t217 + (t272 * t208 + t273 * t209) * t215;
t241 = qJD(3) * t215;
t240 = qJD(4) * t215;
t239 = V_base(5) * qJ(1) + V_base(1);
t235 = qJD(1) + V_base(3);
t200 = qJD(2) * t212 + V_base(4);
t171 = t213 * t241 + t200;
t234 = pkin(2) * t217 + pkin(6) * t215;
t199 = -qJD(2) * t213 + V_base(5);
t233 = rSges(3,1) * t217 - rSges(3,2) * t215;
t232 = Icges(3,1) * t217 - t254;
t231 = -Icges(3,2) * t215 + t253;
t230 = Icges(3,5) * t217 - Icges(3,6) * t215;
t170 = t212 * t241 + t199;
t193 = pkin(1) * t213 + pkin(5) * t212;
t229 = -V_base(4) * qJ(1) + V_base(6) * t193 + V_base(2);
t192 = pkin(1) * t212 - pkin(5) * t213;
t228 = V_base(4) * t192 - t193 * V_base(5) + t235;
t227 = pkin(7) * t215 + t217 * t257;
t176 = t234 * t212;
t198 = t215 * pkin(2) - t217 * pkin(6);
t226 = t199 * t198 + (-t176 - t192) * V_base(6) + t239;
t225 = (-Icges(3,3) * t213 + t212 * t230) * t199 + (Icges(3,3) * t212 + t213 * t230) * t200 + (Icges(3,5) * t215 + Icges(3,6) * t217) * V_base(6);
t177 = t234 * t213;
t224 = V_base(6) * t177 - t198 * t200 + t229;
t223 = t200 * t176 - t177 * t199 + t228;
t133 = -pkin(3) * t249 + t212 * t227;
t140 = -pkin(7) * t217 + t215 * t257;
t202 = -qJD(3) * t217 + V_base(6);
t222 = -t133 * t202 + t170 * t140 + t226;
t134 = pkin(3) * t252 + t213 * t227;
t221 = t202 * t134 - t140 * t171 + t224;
t220 = t171 * t133 - t134 * t170 + t223;
t153 = -Icges(3,6) * t213 + t212 * t231;
t154 = Icges(3,6) * t212 + t213 * t231;
t155 = -Icges(3,5) * t213 + t212 * t232;
t156 = Icges(3,5) * t212 + t213 * t232;
t195 = Icges(3,2) * t217 + t254;
t196 = Icges(3,1) * t215 + t253;
t219 = (-t154 * t215 + t156 * t217) * t200 + (-t153 * t215 + t155 * t217) * t199 + (-t195 * t215 + t196 * t217) * V_base(6);
t207 = Icges(2,4) * t213;
t197 = rSges(3,1) * t215 + rSges(3,2) * t217;
t191 = rSges(2,1) * t213 - rSges(2,2) * t212;
t190 = rSges(2,1) * t212 + rSges(2,2) * t213;
t189 = Icges(2,1) * t213 - t255;
t188 = Icges(2,1) * t212 + t207;
t187 = -Icges(2,2) * t212 + t207;
t186 = Icges(2,2) * t213 + t255;
t183 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t182 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t181 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t180 = V_base(6) + (-qJD(3) - qJD(4)) * t217;
t175 = t213 * t245 + t252;
t174 = t212 * t216 - t213 * t246;
t173 = t212 * t245 - t249;
t172 = -t212 * t246 - t213 * t216;
t166 = -rSges(4,3) * t217 + (rSges(4,1) * t216 - rSges(4,2) * t214) * t215;
t165 = -Icges(4,5) * t217 + (Icges(4,1) * t216 - Icges(4,4) * t214) * t215;
t164 = -Icges(4,6) * t217 + (Icges(4,4) * t216 - Icges(4,2) * t214) * t215;
t163 = -Icges(4,3) * t217 + (Icges(4,5) * t216 - Icges(4,6) * t214) * t215;
t158 = rSges(3,3) * t212 + t213 * t233;
t157 = -rSges(3,3) * t213 + t212 * t233;
t150 = -rSges(5,3) * t217 + (rSges(5,1) * t209 - rSges(5,2) * t208) * t215;
t142 = t213 * t240 + t171;
t141 = t212 * t240 + t170;
t138 = V_base(5) * rSges(2,3) - t190 * V_base(6) + t239;
t137 = t191 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t136 = t190 * V_base(4) - t191 * V_base(5) + t235;
t132 = rSges(4,1) * t175 + rSges(4,2) * t174 + rSges(4,3) * t248;
t131 = rSges(4,1) * t173 + rSges(4,2) * t172 + rSges(4,3) * t251;
t128 = Icges(4,1) * t175 + Icges(4,4) * t174 + Icges(4,5) * t248;
t127 = Icges(4,1) * t173 + Icges(4,4) * t172 + Icges(4,5) * t251;
t126 = Icges(4,4) * t175 + Icges(4,2) * t174 + Icges(4,6) * t248;
t125 = Icges(4,4) * t173 + Icges(4,2) * t172 + Icges(4,6) * t251;
t124 = Icges(4,5) * t175 + Icges(4,6) * t174 + Icges(4,3) * t248;
t123 = Icges(4,5) * t173 + Icges(4,6) * t172 + Icges(4,3) * t251;
t122 = rSges(5,1) * t162 - rSges(5,2) * t161 + rSges(5,3) * t248;
t120 = rSges(5,1) * t160 - rSges(5,2) * t159 + rSges(5,3) * t251;
t104 = t197 * t199 + (-t157 - t192) * V_base(6) + t239;
t103 = t158 * V_base(6) - t197 * t200 + t229;
t102 = t157 * t200 - t158 * t199 + t228;
t101 = -t131 * t202 + t166 * t170 + t226;
t100 = t132 * t202 - t166 * t171 + t224;
t99 = t131 * t171 - t132 * t170 + t223;
t98 = -t120 * t180 + t141 * t150 + t222;
t97 = t122 * t180 - t142 * t150 + t221;
t96 = t120 * t142 - t122 * t141 + t220;
t95 = qJD(5) * t161 + t141 * t242 - t180 * t244 + t222;
t94 = qJD(5) * t159 - t142 * t242 + t180 * t243 + t221;
t93 = qJD(5) * t208 * t215 - t141 * t243 + t142 * t244 + t220;
t1 = m(1) * (t181 ^ 2 + t182 ^ 2 + t183 ^ 2) / 0.2e1 + m(2) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(3) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t200 * (t212 * t225 + t213 * t219) / 0.2e1 + t199 * (t212 * t219 - t213 * t225) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + t171 * ((t124 * t248 + t174 * t126 + t175 * t128) * t171 + (t123 * t248 + t125 * t174 + t127 * t175) * t170 + (t163 * t248 + t164 * t174 + t165 * t175) * t202) / 0.2e1 + t170 * ((t124 * t251 + t126 * t172 + t128 * t173) * t171 + (t123 * t251 + t172 * t125 + t173 * t127) * t170 + (t163 * t251 + t164 * t172 + t165 * t173) * t202) / 0.2e1 + t202 * ((-t123 * t170 - t124 * t171 - t163 * t202) * t217 + ((-t126 * t214 + t128 * t216) * t171 + (-t125 * t214 + t127 * t216) * t170 + (-t164 * t214 + t165 * t216) * t202) * t215) / 0.2e1 + m(5) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + ((t159 * t263 + t160 * t261 + t251 * t262) * t180 + (t159 * t268 + t160 * t264 + t251 * t266) * t142 + (t269 * t159 + t265 * t160 + t267 * t251) * t141) * t141 / 0.2e1 + ((t161 * t263 + t162 * t261 + t248 * t262) * t180 + (t268 * t161 + t264 * t162 + t266 * t248) * t142 + (t161 * t269 + t162 * t265 + t248 * t267) * t141) * t142 / 0.2e1 + ((-t141 * t267 - t142 * t266 - t180 * t262) * t217 + ((t208 * t263 + t209 * t261) * t180 + (t208 * t268 + t209 * t264) * t142 + (t208 * t269 + t209 * t265) * t141) * t215) * t180 / 0.2e1 + ((-t186 * t212 + t188 * t213 + Icges(1,4)) * V_base(5) + (-t187 * t212 + t189 * t213 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t186 * t213 + t188 * t212 + Icges(1,2)) * V_base(5) + (t187 * t213 + t189 * t212 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t154 * t217 + t156 * t215) * t200 + (t153 * t217 + t155 * t215) * t199 + (t195 * t217 + t196 * t215 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t213 - Icges(2,6) * t212 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t212 + Icges(2,6) * t213 + Icges(1,6));
T = t1;
