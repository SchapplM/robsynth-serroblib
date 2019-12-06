% Calculate kinetic energy for
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:30
% EndTime: 2019-12-05 15:37:32
% DurationCPUTime: 2.67s
% Computational Cost: add. (1228->280), mult. (1871->398), div. (0->0), fcn. (1831->8), ass. (0->131)
t273 = Icges(5,1) + Icges(6,1);
t272 = -Icges(5,4) + Icges(6,5);
t271 = Icges(6,4) + Icges(5,5);
t270 = Icges(5,2) + Icges(6,3);
t269 = -Icges(6,6) + Icges(5,6);
t268 = -Icges(5,3) - Icges(6,2);
t267 = rSges(6,1) + pkin(4);
t266 = rSges(6,3) + qJ(5);
t204 = pkin(8) + qJ(4);
t200 = sin(t204);
t201 = cos(t204);
t208 = cos(pkin(7));
t206 = sin(pkin(7));
t211 = cos(qJ(2));
t243 = t206 * t211;
t151 = t200 * t243 + t201 * t208;
t152 = -t200 * t208 + t201 * t243;
t210 = sin(qJ(2));
t244 = t206 * t210;
t263 = t151 * t270 + t152 * t272 - t244 * t269;
t240 = t208 * t211;
t153 = t200 * t240 - t201 * t206;
t154 = t200 * t206 + t201 * t240;
t241 = t208 * t210;
t262 = t153 * t270 + t154 * t272 - t241 * t269;
t261 = -t151 * t269 + t152 * t271 - t244 * t268;
t260 = -t153 * t269 + t154 * t271 - t241 * t268;
t259 = t151 * t272 + t152 * t273 + t244 * t271;
t258 = t153 * t272 + t154 * t273 + t241 * t271;
t257 = t269 * t211 + (t200 * t270 + t201 * t272) * t210;
t256 = t268 * t211 + (-t200 * t269 + t201 * t271) * t210;
t255 = -t271 * t211 + (t200 * t272 + t201 * t273) * t210;
t207 = cos(pkin(8));
t249 = pkin(3) * t207;
t248 = Icges(2,4) * t206;
t247 = Icges(3,4) * t210;
t246 = Icges(3,4) * t211;
t205 = sin(pkin(8));
t245 = t206 * t205;
t242 = t208 * t205;
t238 = rSges(6,2) * t244 + t151 * t266 + t152 * t267;
t237 = rSges(6,2) * t241 + t153 * t266 + t154 * t267;
t236 = -rSges(6,2) * t211 + (t200 * t266 + t201 * t267) * t210;
t225 = pkin(2) * t211 + qJ(3) * t210;
t170 = t225 * t206;
t185 = pkin(1) * t206 - pkin(5) * t208;
t235 = -t170 - t185;
t234 = qJD(3) * t210;
t233 = qJD(4) * t210;
t232 = V_base(5) * qJ(1) + V_base(1);
t228 = qJD(1) + V_base(3);
t194 = qJD(2) * t206 + V_base(4);
t190 = pkin(2) * t210 - qJ(3) * t211;
t193 = -qJD(2) * t208 + V_base(5);
t227 = t190 * t193 + t208 * t234 + t232;
t226 = rSges(3,1) * t211 - rSges(3,2) * t210;
t224 = Icges(3,1) * t211 - t247;
t223 = -Icges(3,2) * t210 + t246;
t222 = Icges(3,5) * t211 - Icges(3,6) * t210;
t186 = pkin(1) * t208 + pkin(5) * t206;
t221 = -V_base(4) * qJ(1) + t186 * V_base(6) + V_base(2);
t220 = t185 * V_base(4) - t186 * V_base(5) + t228;
t219 = pkin(6) * t210 + t211 * t249;
t171 = t225 * t208;
t218 = t171 * V_base(6) + t206 * t234 + t221;
t217 = (-Icges(3,3) * t208 + t206 * t222) * t193 + (Icges(3,3) * t206 + t208 * t222) * t194 + (Icges(3,5) * t210 + Icges(3,6) * t211) * V_base(6);
t216 = -qJD(3) * t211 + t170 * t194 + t220;
t129 = -pkin(3) * t242 + t206 * t219;
t133 = -pkin(6) * t211 + t210 * t249;
t215 = t193 * t133 + (-t129 + t235) * V_base(6) + t227;
t130 = pkin(3) * t245 + t208 * t219;
t214 = V_base(6) * t130 + (-t133 - t190) * t194 + t218;
t213 = t194 * t129 + (-t130 - t171) * t193 + t216;
t147 = -Icges(3,6) * t208 + t206 * t223;
t148 = Icges(3,6) * t206 + t208 * t223;
t149 = -Icges(3,5) * t208 + t206 * t224;
t150 = Icges(3,5) * t206 + t208 * t224;
t188 = Icges(3,2) * t211 + t247;
t189 = Icges(3,1) * t210 + t246;
t212 = (-t148 * t210 + t150 * t211) * t194 + (-t147 * t210 + t149 * t211) * t193 + (-t188 * t210 + t189 * t211) * V_base(6);
t202 = Icges(2,4) * t208;
t195 = -qJD(4) * t211 + V_base(6);
t191 = rSges(3,1) * t210 + rSges(3,2) * t211;
t184 = rSges(2,1) * t208 - rSges(2,2) * t206;
t183 = rSges(2,1) * t206 + rSges(2,2) * t208;
t182 = Icges(2,1) * t208 - t248;
t181 = Icges(2,1) * t206 + t202;
t180 = -Icges(2,2) * t206 + t202;
t179 = Icges(2,2) * t208 + t248;
t176 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t175 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t174 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t169 = t208 * t233 + t194;
t168 = t206 * t233 + t193;
t167 = t207 * t240 + t245;
t166 = -t205 * t240 + t206 * t207;
t165 = t207 * t243 - t242;
t164 = -t205 * t243 - t207 * t208;
t160 = -rSges(4,3) * t211 + (rSges(4,1) * t207 - rSges(4,2) * t205) * t210;
t159 = rSges(3,3) * t206 + t208 * t226;
t158 = -rSges(3,3) * t208 + t206 * t226;
t157 = -Icges(4,5) * t211 + (Icges(4,1) * t207 - Icges(4,4) * t205) * t210;
t156 = -Icges(4,6) * t211 + (Icges(4,4) * t207 - Icges(4,2) * t205) * t210;
t155 = -Icges(4,3) * t211 + (Icges(4,5) * t207 - Icges(4,6) * t205) * t210;
t144 = -rSges(5,3) * t211 + (rSges(5,1) * t201 - rSges(5,2) * t200) * t210;
t135 = V_base(5) * rSges(2,3) - t183 * V_base(6) + t232;
t134 = t184 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t132 = t183 * V_base(4) - t184 * V_base(5) + t228;
t128 = rSges(4,1) * t167 + rSges(4,2) * t166 + rSges(4,3) * t241;
t127 = rSges(4,1) * t165 + rSges(4,2) * t164 + rSges(4,3) * t244;
t124 = Icges(4,1) * t167 + Icges(4,4) * t166 + Icges(4,5) * t241;
t123 = Icges(4,1) * t165 + Icges(4,4) * t164 + Icges(4,5) * t244;
t122 = Icges(4,4) * t167 + Icges(4,2) * t166 + Icges(4,6) * t241;
t121 = Icges(4,4) * t165 + Icges(4,2) * t164 + Icges(4,6) * t244;
t120 = Icges(4,5) * t167 + Icges(4,6) * t166 + Icges(4,3) * t241;
t119 = Icges(4,5) * t165 + Icges(4,6) * t164 + Icges(4,3) * t244;
t117 = rSges(5,1) * t154 - rSges(5,2) * t153 + rSges(5,3) * t241;
t115 = rSges(5,1) * t152 - rSges(5,2) * t151 + rSges(5,3) * t244;
t100 = t191 * t193 + (-t158 - t185) * V_base(6) + t232;
t99 = t159 * V_base(6) - t191 * t194 + t221;
t98 = t158 * t194 - t159 * t193 + t220;
t97 = t160 * t193 + (-t127 + t235) * V_base(6) + t227;
t96 = t128 * V_base(6) + (-t160 - t190) * t194 + t218;
t95 = t194 * t127 + (-t128 - t171) * t193 + t216;
t94 = -t115 * t195 + t144 * t168 + t215;
t93 = t117 * t195 - t144 * t169 + t214;
t92 = t115 * t169 - t117 * t168 + t213;
t91 = qJD(5) * t153 + t168 * t236 - t195 * t238 + t215;
t90 = qJD(5) * t151 - t169 * t236 + t195 * t237 + t214;
t89 = qJD(5) * t200 * t210 - t168 * t237 + t169 * t238 + t213;
t1 = m(1) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(2) * (t132 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + ((t151 * t257 + t152 * t255 + t244 * t256) * t195 + (t151 * t262 + t152 * t258 + t244 * t260) * t169 + (t263 * t151 + t259 * t152 + t261 * t244) * t168) * t168 / 0.2e1 + ((t153 * t257 + t154 * t255 + t241 * t256) * t195 + (t262 * t153 + t258 * t154 + t260 * t241) * t169 + (t153 * t263 + t154 * t259 + t241 * t261) * t168) * t169 / 0.2e1 + (t206 * t212 - t208 * t217 + (t120 * t244 + t122 * t164 + t124 * t165) * t194 + (t119 * t244 + t121 * t164 + t123 * t165) * t193 + (t155 * t244 + t156 * t164 + t157 * t165) * V_base(6)) * t193 / 0.2e1 + (t206 * t217 + t208 * t212 + (t120 * t241 + t122 * t166 + t124 * t167) * t194 + (t119 * t241 + t121 * t166 + t123 * t167) * t193 + (t155 * t241 + t156 * t166 + t157 * t167) * V_base(6)) * t194 / 0.2e1 + ((-t168 * t261 - t169 * t260 - t195 * t256) * t211 + ((t200 * t257 + t201 * t255) * t195 + (t200 * t262 + t201 * t258) * t169 + (t200 * t263 + t201 * t259) * t168) * t210) * t195 / 0.2e1 + ((-t179 * t206 + t181 * t208 + Icges(1,4)) * V_base(5) + (-t180 * t206 + t182 * t208 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t179 * t208 + t181 * t206 + Icges(1,2)) * V_base(5) + (t180 * t208 + t182 * t206 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t148 * t211 + t210 * t150) * t194 + (t147 * t211 + t210 * t149) * t193 + (-t119 * t193 - t120 * t194) * t211 + ((-t122 * t205 + t124 * t207) * t194 + (-t121 * t205 + t123 * t207) * t193) * t210 + (Icges(1,3) + Icges(2,3) + (t188 - t155) * t211 + (-t156 * t205 + t157 * t207 + t189) * t210) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t208 - Icges(2,6) * t206 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t206 + Icges(2,6) * t208 + Icges(1,6));
T = t1;
