% Calculate kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:43:59
% EndTime: 2019-12-05 17:44:02
% DurationCPUTime: 3.30s
% Computational Cost: add. (1399->337), mult. (1865->476), div. (0->0), fcn. (1815->10), ass. (0->156)
t222 = sin(qJ(1));
t223 = cos(qJ(1));
t218 = sin(pkin(8));
t220 = cos(pkin(8));
t264 = Icges(3,4) * t220;
t235 = -Icges(3,2) * t218 + t264;
t156 = Icges(3,6) * t223 - t222 * t235;
t265 = Icges(3,4) * t218;
t236 = Icges(3,1) * t220 - t265;
t158 = Icges(3,5) * t223 - t222 * t236;
t266 = Icges(2,4) * t223;
t275 = -Icges(2,1) * t222 - t156 * t218 + t158 * t220 - t266;
t157 = Icges(3,6) * t222 + t223 * t235;
t159 = Icges(3,5) * t222 + t223 * t236;
t267 = Icges(2,4) * t222;
t274 = Icges(2,1) * t223 - t157 * t218 + t159 * t220 - t267;
t216 = pkin(9) + qJ(4);
t208 = cos(t216);
t250 = pkin(4) * t208;
t273 = pkin(7) * t218 + t220 * t250;
t219 = cos(pkin(9));
t269 = t219 * pkin(3);
t272 = pkin(6) * t218 + t220 * t269;
t189 = pkin(2) * t218 - qJ(3) * t220;
t268 = -pkin(5) - t189;
t217 = sin(pkin(9));
t263 = t217 * t223;
t262 = t218 * t222;
t261 = t218 * t223;
t260 = t220 * t223;
t209 = qJ(5) + t216;
t204 = sin(t209);
t259 = t222 * t204;
t205 = cos(t209);
t258 = t222 * t205;
t207 = sin(t216);
t257 = t222 * t207;
t256 = t222 * t208;
t255 = t222 * t217;
t254 = t222 * t219;
t237 = pkin(2) * t220 + qJ(3) * t218;
t173 = t237 * t222;
t198 = -t222 * pkin(1) + qJ(2) * t223;
t252 = t173 - t198;
t174 = t237 * t223;
t200 = pkin(1) * t223 + t222 * qJ(2);
t251 = -t174 - t200;
t248 = qJD(3) * t218;
t247 = qJD(4) * t218;
t246 = -qJD(4) - qJD(5);
t245 = V_base(5) * t200 + V_base(1);
t244 = V_base(6) * pkin(5) + V_base(2);
t184 = t223 * t247 + V_base(6);
t210 = V_base(4) + qJD(1);
t241 = pkin(4) * t207;
t240 = qJD(2) * t222 + t210 * t198 + V_base(3);
t239 = qJD(2) * t223 + t244;
t238 = rSges(3,1) * t220 - rSges(3,2) * t218;
t234 = Icges(3,5) * t220 - Icges(3,6) * t218;
t187 = Icges(3,2) * t220 + t265;
t188 = Icges(3,1) * t218 + t264;
t231 = t187 * t218 - t188 * t220;
t230 = -t210 * t173 + t223 * t248 + t240;
t229 = -qJD(3) * t220 + V_base(5) * t174 + t245;
t228 = V_base(6) * t189 - t222 * t248 + t239;
t227 = (Icges(3,3) * t223 - t222 * t234) * V_base(5) + (Icges(3,3) * t222 + t223 * t234) * V_base(6) + (Icges(3,5) * t218 + Icges(3,6) * t220) * t210;
t127 = pkin(3) * t263 - t222 * t272;
t138 = -pkin(6) * t220 + t218 * t269;
t226 = t210 * t127 + (-t138 + t268) * V_base(5) + t230;
t128 = pkin(3) * t255 + t223 * t272;
t225 = V_base(5) * t128 + (-t127 + t252) * V_base(6) + t229;
t224 = V_base(6) * t138 + (-t128 + t251) * t210 + t228;
t201 = rSges(2,1) * t223 - t222 * rSges(2,2);
t199 = -t222 * rSges(2,1) - rSges(2,2) * t223;
t195 = -Icges(2,2) * t222 + t266;
t194 = -Icges(2,2) * t223 - t267;
t193 = Icges(2,5) * t223 - Icges(2,6) * t222;
t192 = -Icges(2,5) * t222 - Icges(2,6) * t223;
t191 = -qJD(4) * t220 + t210;
t190 = rSges(3,1) * t218 + rSges(3,2) * t220;
t185 = -t222 * t247 + V_base(5);
t183 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t182 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t181 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t175 = t220 * t246 + t210;
t172 = t219 * t260 + t255;
t171 = -t217 * t260 + t254;
t170 = -t220 * t254 + t263;
t169 = t219 * t223 + t220 * t255;
t167 = t246 * t262 + V_base(5);
t166 = qJD(5) * t261 + t184;
t165 = t222 * rSges(3,3) + t223 * t238;
t164 = rSges(3,3) * t223 - t222 * t238;
t163 = t208 * t260 + t257;
t162 = -t207 * t260 + t256;
t161 = t207 * t223 - t220 * t256;
t160 = t208 * t223 + t220 * t257;
t152 = -rSges(4,3) * t220 + (rSges(4,1) * t219 - rSges(4,2) * t217) * t218;
t151 = -Icges(4,5) * t220 + (Icges(4,1) * t219 - Icges(4,4) * t217) * t218;
t150 = -Icges(4,6) * t220 + (Icges(4,4) * t219 - Icges(4,2) * t217) * t218;
t149 = -Icges(4,3) * t220 + (Icges(4,5) * t219 - Icges(4,6) * t217) * t218;
t148 = t205 * t260 + t259;
t147 = -t204 * t260 + t258;
t146 = t204 * t223 - t220 * t258;
t145 = t205 * t223 + t220 * t259;
t144 = -rSges(5,3) * t220 + (rSges(5,1) * t208 - rSges(5,2) * t207) * t218;
t143 = -Icges(5,5) * t220 + (Icges(5,1) * t208 - Icges(5,4) * t207) * t218;
t142 = -Icges(5,6) * t220 + (Icges(5,4) * t208 - Icges(5,2) * t207) * t218;
t141 = -Icges(5,3) * t220 + (Icges(5,5) * t208 - Icges(5,6) * t207) * t218;
t140 = V_base(6) * rSges(2,3) - t201 * t210 + t244;
t139 = t199 * t210 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t137 = -rSges(6,3) * t220 + (rSges(6,1) * t205 - rSges(6,2) * t204) * t218;
t136 = -Icges(6,5) * t220 + (Icges(6,1) * t205 - Icges(6,4) * t204) * t218;
t135 = -Icges(6,6) * t220 + (Icges(6,4) * t205 - Icges(6,2) * t204) * t218;
t134 = -Icges(6,3) * t220 + (Icges(6,5) * t205 - Icges(6,6) * t204) * t218;
t133 = -t199 * V_base(6) + t201 * V_base(5) + V_base(1);
t131 = -pkin(7) * t220 + t218 * t250;
t130 = t172 * rSges(4,1) + t171 * rSges(4,2) + rSges(4,3) * t261;
t129 = rSges(4,1) * t170 + rSges(4,2) * t169 - rSges(4,3) * t262;
t126 = Icges(4,1) * t172 + Icges(4,4) * t171 + Icges(4,5) * t261;
t125 = Icges(4,1) * t170 + Icges(4,4) * t169 - Icges(4,5) * t262;
t124 = Icges(4,4) * t172 + Icges(4,2) * t171 + Icges(4,6) * t261;
t123 = Icges(4,4) * t170 + Icges(4,2) * t169 - Icges(4,6) * t262;
t122 = Icges(4,5) * t172 + Icges(4,6) * t171 + Icges(4,3) * t261;
t121 = Icges(4,5) * t170 + Icges(4,6) * t169 - Icges(4,3) * t262;
t118 = t163 * rSges(5,1) + t162 * rSges(5,2) + rSges(5,3) * t261;
t117 = rSges(5,1) * t161 + rSges(5,2) * t160 - rSges(5,3) * t262;
t116 = Icges(5,1) * t163 + Icges(5,4) * t162 + Icges(5,5) * t261;
t115 = Icges(5,1) * t161 + Icges(5,4) * t160 - Icges(5,5) * t262;
t114 = Icges(5,4) * t163 + Icges(5,2) * t162 + Icges(5,6) * t261;
t113 = Icges(5,4) * t161 + Icges(5,2) * t160 - Icges(5,6) * t262;
t112 = Icges(5,5) * t163 + Icges(5,6) * t162 + Icges(5,3) * t261;
t111 = Icges(5,5) * t161 + Icges(5,6) * t160 - Icges(5,3) * t262;
t110 = t148 * rSges(6,1) + t147 * rSges(6,2) + rSges(6,3) * t261;
t109 = rSges(6,1) * t146 + rSges(6,2) * t145 - rSges(6,3) * t262;
t108 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t261;
t107 = Icges(6,1) * t146 + Icges(6,4) * t145 - Icges(6,5) * t262;
t106 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t261;
t105 = Icges(6,4) * t146 + Icges(6,2) * t145 - Icges(6,6) * t262;
t104 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t261;
t103 = Icges(6,5) * t146 + Icges(6,6) * t145 - Icges(6,3) * t262;
t102 = t190 * V_base(6) + (-t165 - t200) * t210 + t239;
t101 = t164 * t210 + (-pkin(5) - t190) * V_base(5) + t240;
t100 = t241 * t222 + t223 * t273;
t99 = -t222 * t273 + t241 * t223;
t98 = t165 * V_base(5) + (-t164 - t198) * V_base(6) + t245;
t97 = t152 * V_base(6) + (-t130 + t251) * t210 + t228;
t96 = t129 * t210 + (-t152 + t268) * V_base(5) + t230;
t95 = t130 * V_base(5) + (-t129 + t252) * V_base(6) + t229;
t94 = -t118 * t191 + t144 * t184 + t224;
t93 = t117 * t191 - t144 * t185 + t226;
t92 = -t117 * t184 + t118 * t185 + t225;
t91 = -t100 * t191 - t110 * t175 + t131 * t184 + t137 * t166 + t224;
t90 = t109 * t175 - t131 * t185 - t137 * t167 + t191 * t99 + t226;
t89 = t100 * t185 - t109 * t166 + t110 * t167 - t184 * t99 + t225;
t1 = m(1) * (t181 ^ 2 + t182 ^ 2 + t183 ^ 2) / 0.2e1 + m(2) * (t133 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(3) * (t101 ^ 2 + t102 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t191 * ((-t111 * t185 - t112 * t184 - t141 * t191) * t220 + ((-t142 * t207 + t143 * t208) * t191 + (-t113 * t207 + t115 * t208) * t185 + (-t114 * t207 + t116 * t208) * t184) * t218) / 0.2e1 + t185 * ((-t141 * t262 + t142 * t160 + t143 * t161) * t191 + (-t111 * t262 + t113 * t160 + t115 * t161) * t185 + (-t112 * t262 + t114 * t160 + t116 * t161) * t184) / 0.2e1 + t184 * ((t141 * t261 + t162 * t142 + t163 * t143) * t191 + (t111 * t261 + t162 * t113 + t163 * t115) * t185 + (t112 * t261 + t162 * t114 + t163 * t116) * t184) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t175 * ((-t103 * t167 - t104 * t166 - t134 * t175) * t220 + ((-t135 * t204 + t136 * t205) * t175 + (-t105 * t204 + t107 * t205) * t167 + (-t106 * t204 + t108 * t205) * t166) * t218) / 0.2e1 + t167 * ((-t134 * t262 + t135 * t145 + t136 * t146) * t175 + (-t103 * t262 + t105 * t145 + t107 * t146) * t167 + (-t104 * t262 + t106 * t145 + t108 * t146) * t166) / 0.2e1 + t166 * ((t134 * t261 + t147 * t135 + t148 * t136) * t175 + (t103 * t261 + t147 * t105 + t148 * t107) * t167 + (t104 * t261 + t147 * t106 + t148 * t108) * t166) / 0.2e1 + ((t193 + (t157 - t122) * t220 + (-t124 * t217 + t126 * t219 + t159) * t218) * V_base(6) + (t192 + (t156 - t121) * t220 + (-t123 * t217 + t125 * t219 + t158) * t218) * V_base(5) + (Icges(2,3) + (t187 - t149) * t220 + (-t150 * t217 + t151 * t219 + t188) * t218) * t210) * t210 / 0.2e1 + (t227 * t223 + (-t149 * t262 + t150 * t169 + t151 * t170 + t231 * t222 + t192) * t210 + (-t122 * t262 + t124 * t169 + t126 * t170 - t195 * t223 - t222 * t274 + Icges(1,6)) * V_base(6) + (-t121 * t262 + t123 * t169 + t125 * t170 - t194 * t223 - t275 * t222 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t227 * t222 + (t149 * t261 + t171 * t150 + t172 * t151 - t231 * t223 + t193) * t210 + (t122 * t261 + t171 * t124 + t172 * t126 - t222 * t195 + t274 * t223 + Icges(1,3)) * V_base(6) + (t121 * t261 + t171 * t123 + t172 * t125 - t222 * t194 + t223 * t275 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
