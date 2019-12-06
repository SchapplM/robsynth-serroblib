% Calculate kinetic energy for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:25
% EndTime: 2019-12-05 18:05:28
% DurationCPUTime: 3.22s
% Computational Cost: add. (1307->287), mult. (1945->410), div. (0->0), fcn. (1895->8), ass. (0->138)
t277 = Icges(5,1) + Icges(6,1);
t276 = Icges(5,4) + Icges(6,4);
t275 = -Icges(6,5) - Icges(5,5);
t274 = Icges(5,2) + Icges(6,2);
t273 = -Icges(6,6) - Icges(5,6);
t272 = -Icges(6,3) - Icges(5,3);
t202 = qJ(3) + qJ(4);
t197 = sin(t202);
t198 = cos(t202);
t208 = cos(qJ(1));
t204 = cos(pkin(8));
t206 = sin(qJ(1));
t243 = t204 * t206;
t152 = t197 * t243 + t198 * t208;
t153 = t197 * t208 - t198 * t243;
t203 = sin(pkin(8));
t245 = t203 * t206;
t271 = -t273 * t152 - t275 * t153 + t272 * t245;
t242 = t204 * t208;
t154 = -t197 * t242 + t198 * t206;
t155 = t197 * t206 + t198 * t242;
t244 = t203 * t208;
t270 = -t273 * t154 - t275 * t155 - t272 * t244;
t269 = t274 * t152 + t276 * t153 + t273 * t245;
t268 = t274 * t154 + t276 * t155 - t273 * t244;
t267 = t276 * t152 + t277 * t153 + t275 * t245;
t266 = t276 * t154 + t277 * t155 - t275 * t244;
t265 = t272 * t204 + (t273 * t197 - t275 * t198) * t203;
t264 = t273 * t204 + (-t274 * t197 + t276 * t198) * t203;
t263 = t275 * t204 + (-t276 * t197 + t277 * t198) * t203;
t246 = Icges(3,4) * t204;
t221 = -Icges(3,2) * t203 + t246;
t145 = Icges(3,6) * t208 - t206 * t221;
t247 = Icges(3,4) * t203;
t222 = Icges(3,1) * t204 - t247;
t147 = Icges(3,5) * t208 - t206 * t222;
t248 = Icges(2,4) * t208;
t262 = -Icges(2,1) * t206 - t145 * t203 + t147 * t204 - t248;
t146 = Icges(3,6) * t206 + t208 * t221;
t148 = Icges(3,5) * t206 + t208 * t222;
t249 = Icges(2,4) * t206;
t261 = Icges(2,1) * t208 - t146 * t203 + t148 * t204 - t249;
t236 = pkin(4) * t198;
t260 = qJ(5) * t203 + t204 * t236;
t207 = cos(qJ(3));
t253 = t207 * pkin(3);
t259 = pkin(7) * t203 + t204 * t253;
t227 = pkin(4) * t197;
t251 = rSges(6,1) * t153 + rSges(6,2) * t152 - rSges(6,3) * t245 - t206 * t260 + t227 * t208;
t250 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t244 + t227 * t206 + t208 * t260;
t205 = sin(qJ(3));
t241 = t205 * t206;
t240 = t205 * t208;
t239 = t206 * t207;
t238 = t207 * t208;
t237 = (-qJ(5) - rSges(6,3)) * t204 + (rSges(6,1) * t198 - rSges(6,2) * t197 + t236) * t203;
t234 = qJD(3) * t203;
t233 = qJD(5) * t203;
t232 = -qJD(3) - qJD(4);
t190 = pkin(1) * t208 + qJ(2) * t206;
t231 = V_base(5) * t190 + V_base(1);
t230 = V_base(6) * pkin(5) + V_base(2);
t174 = t208 * t234 + V_base(6);
t194 = V_base(4) + qJD(1);
t188 = -pkin(1) * t206 + qJ(2) * t208;
t226 = qJD(2) * t206 + t194 * t188 + V_base(3);
t225 = qJD(2) * t208 + t230;
t224 = pkin(2) * t204 + pkin(6) * t203;
t223 = rSges(3,1) * t204 - rSges(3,2) * t203;
t220 = Icges(3,5) * t204 - Icges(3,6) * t203;
t177 = Icges(3,2) * t204 + t247;
t178 = Icges(3,1) * t203 + t246;
t217 = t177 * t203 - t178 * t204;
t164 = t224 * t208;
t181 = t203 * pkin(2) - t204 * pkin(6);
t216 = V_base(6) * t181 + (-t164 - t190) * t194 + t225;
t163 = t224 * t206;
t215 = V_base(5) * t164 + (t163 - t188) * V_base(6) + t231;
t214 = -t194 * t163 + (-pkin(5) - t181) * V_base(5) + t226;
t213 = (Icges(3,3) * t208 - t206 * t220) * V_base(5) + (Icges(3,3) * t206 + t208 * t220) * V_base(6) + (Icges(3,5) * t203 + Icges(3,6) * t204) * t194;
t121 = pkin(3) * t240 - t206 * t259;
t122 = pkin(3) * t241 + t208 * t259;
t175 = -t206 * t234 + V_base(5);
t212 = -t121 * t174 + t175 * t122 + t215;
t130 = -pkin(7) * t204 + t203 * t253;
t180 = -qJD(3) * t204 + t194;
t211 = -t122 * t180 + t174 * t130 + t216;
t210 = t180 * t121 - t130 * t175 + t214;
t191 = rSges(2,1) * t208 - rSges(2,2) * t206;
t189 = -rSges(2,1) * t206 - rSges(2,2) * t208;
t185 = -Icges(2,2) * t206 + t248;
t184 = -Icges(2,2) * t208 - t249;
t183 = Icges(2,5) * t208 - Icges(2,6) * t206;
t182 = -Icges(2,5) * t206 - Icges(2,6) * t208;
t179 = rSges(3,1) * t203 + rSges(3,2) * t204;
t172 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t171 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t170 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t165 = t204 * t232 + t194;
t162 = t204 * t238 + t241;
t161 = -t204 * t240 + t239;
t160 = -t204 * t239 + t240;
t159 = t204 * t241 + t238;
t157 = t232 * t245 + V_base(5);
t156 = qJD(4) * t244 + t174;
t151 = rSges(3,3) * t206 + t208 * t223;
t150 = rSges(3,3) * t208 - t206 * t223;
t149 = -rSges(4,3) * t204 + (rSges(4,1) * t207 - rSges(4,2) * t205) * t203;
t141 = -Icges(4,5) * t204 + (Icges(4,1) * t207 - Icges(4,4) * t205) * t203;
t140 = -Icges(4,6) * t204 + (Icges(4,4) * t207 - Icges(4,2) * t205) * t203;
t139 = -Icges(4,3) * t204 + (Icges(4,5) * t207 - Icges(4,6) * t205) * t203;
t138 = -rSges(5,3) * t204 + (rSges(5,1) * t198 - rSges(5,2) * t197) * t203;
t129 = V_base(6) * rSges(2,3) - t191 * t194 + t230;
t128 = t189 * t194 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t127 = -t189 * V_base(6) + t191 * V_base(5) + V_base(1);
t124 = rSges(4,1) * t162 + rSges(4,2) * t161 + rSges(4,3) * t244;
t123 = rSges(4,1) * t160 + rSges(4,2) * t159 - rSges(4,3) * t245;
t120 = Icges(4,1) * t162 + Icges(4,4) * t161 + Icges(4,5) * t244;
t119 = Icges(4,1) * t160 + Icges(4,4) * t159 - Icges(4,5) * t245;
t118 = Icges(4,4) * t162 + Icges(4,2) * t161 + Icges(4,6) * t244;
t117 = Icges(4,4) * t160 + Icges(4,2) * t159 - Icges(4,6) * t245;
t116 = Icges(4,5) * t162 + Icges(4,6) * t161 + Icges(4,3) * t244;
t115 = Icges(4,5) * t160 + Icges(4,6) * t159 - Icges(4,3) * t245;
t114 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t244;
t112 = rSges(5,1) * t153 + rSges(5,2) * t152 - rSges(5,3) * t245;
t96 = t179 * V_base(6) + (-t151 - t190) * t194 + t225;
t95 = t150 * t194 + (-pkin(5) - t179) * V_base(5) + t226;
t92 = t151 * V_base(5) + (-t150 - t188) * V_base(6) + t231;
t91 = -t124 * t180 + t149 * t174 + t216;
t90 = t123 * t180 - t149 * t175 + t214;
t89 = -t123 * t174 + t124 * t175 + t215;
t88 = -t114 * t165 + t138 * t156 + t211;
t87 = t112 * t165 - t138 * t157 + t210;
t86 = -t112 * t156 + t114 * t157 + t212;
t85 = t156 * t237 - t165 * t250 - t206 * t233 + t211;
t84 = -t157 * t237 + t165 * t251 + t208 * t233 + t210;
t83 = -qJD(5) * t204 - t156 * t251 + t157 * t250 + t212;
t1 = m(1) * (t170 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(2) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t180 * ((-t115 * t175 - t116 * t174 - t139 * t180) * t204 + ((-t140 * t205 + t141 * t207) * t180 + (-t117 * t205 + t119 * t207) * t175 + (-t118 * t205 + t120 * t207) * t174) * t203) / 0.2e1 + t175 * ((-t139 * t245 + t140 * t159 + t141 * t160) * t180 + (-t115 * t245 + t117 * t159 + t119 * t160) * t175 + (-t116 * t245 + t118 * t159 + t120 * t160) * t174) / 0.2e1 + t174 * ((t139 * t244 + t140 * t161 + t141 * t162) * t180 + (t115 * t244 + t117 * t161 + t119 * t162) * t175 + (t116 * t244 + t118 * t161 + t120 * t162) * t174) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t154 * t264 + t155 * t263 + t265 * t244) * t165 + (t269 * t154 + t267 * t155 + t244 * t271) * t157 + (t268 * t154 + t266 * t155 + t270 * t244) * t156) * t156 / 0.2e1 + ((t152 * t264 + t153 * t263 - t245 * t265) * t165 + (t269 * t152 + t267 * t153 - t271 * t245) * t157 + (t152 * t268 + t153 * t266 - t245 * t270) * t156) * t157 / 0.2e1 + ((-t270 * t156 - t157 * t271 - t265 * t165) * t204 + ((-t197 * t264 + t198 * t263) * t165 + (-t197 * t269 + t198 * t267) * t157 + (-t197 * t268 + t198 * t266) * t156) * t203) * t165 / 0.2e1 + ((t146 * t204 + t148 * t203 + t183) * V_base(6) + (t145 * t204 + t147 * t203 + t182) * V_base(5) + (t177 * t204 + t178 * t203 + Icges(2,3)) * t194) * t194 / 0.2e1 + (t213 * t208 + (t217 * t206 + t182) * t194 + (-t185 * t208 - t206 * t261 + Icges(1,6)) * V_base(6) + (-t184 * t208 - t262 * t206 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t213 * t206 + (-t217 * t208 + t183) * t194 + (-t185 * t206 + t261 * t208 + Icges(1,3)) * V_base(6) + (-t184 * t206 + t208 * t262 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
