% Calculate kinetic energy for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:46
% EndTime: 2019-12-05 17:55:50
% DurationCPUTime: 3.63s
% Computational Cost: add. (1426->329), mult. (1910->466), div. (0->0), fcn. (1860->10), ass. (0->152)
t280 = -Icges(5,3) - Icges(4,3);
t216 = qJ(3) + pkin(9);
t208 = cos(t216);
t218 = cos(pkin(8));
t223 = cos(qJ(1));
t207 = sin(t216);
t221 = sin(qJ(1));
t255 = t221 * t207;
t160 = t208 * t223 + t218 * t255;
t254 = t221 * t208;
t161 = t207 * t223 - t218 * t254;
t222 = cos(qJ(3));
t251 = t222 * t223;
t220 = sin(qJ(3));
t253 = t221 * t220;
t169 = t218 * t253 + t251;
t252 = t221 * t222;
t258 = t220 * t223;
t170 = -t218 * t252 + t258;
t217 = sin(pkin(8));
t261 = t217 * t221;
t279 = Icges(4,5) * t170 + Icges(5,5) * t161 + Icges(4,6) * t169 + Icges(5,6) * t160 + t280 * t261;
t259 = t218 * t223;
t162 = -t207 * t259 + t254;
t163 = t208 * t259 + t255;
t171 = -t218 * t258 + t252;
t172 = t218 * t251 + t253;
t260 = t217 * t223;
t278 = Icges(4,5) * t172 + Icges(5,5) * t163 + Icges(4,6) * t171 + Icges(5,6) * t162 - t280 * t260;
t277 = t280 * t218 + (Icges(4,5) * t222 + Icges(5,5) * t208 - Icges(4,6) * t220 - Icges(5,6) * t207) * t217;
t262 = Icges(3,4) * t218;
t235 = -Icges(3,2) * t217 + t262;
t155 = Icges(3,6) * t223 - t221 * t235;
t263 = Icges(3,4) * t217;
t236 = Icges(3,1) * t218 - t263;
t157 = Icges(3,5) * t223 - t221 * t236;
t264 = Icges(2,4) * t223;
t276 = -Icges(2,1) * t221 - t155 * t217 + t157 * t218 - t264;
t156 = Icges(3,6) * t221 + t223 * t235;
t158 = Icges(3,5) * t221 + t223 * t236;
t265 = Icges(2,4) * t221;
t275 = Icges(2,1) * t223 - t156 * t217 + t158 * t218 - t265;
t250 = pkin(4) * t208;
t274 = pkin(7) * t217 + t218 * t250;
t267 = t222 * pkin(3);
t273 = qJ(4) * t217 + t218 * t267;
t209 = qJ(5) + t216;
t204 = sin(t209);
t257 = t221 * t204;
t205 = cos(t209);
t256 = t221 * t205;
t248 = qJD(3) * t217;
t247 = qJD(4) * t217;
t246 = -qJD(3) - qJD(5);
t200 = pkin(1) * t223 + t221 * qJ(2);
t245 = V_base(5) * t200 + V_base(1);
t244 = V_base(6) * pkin(5) + V_base(2);
t184 = t223 * t248 + V_base(6);
t210 = V_base(4) + qJD(1);
t241 = pkin(4) * t207;
t198 = -t221 * pkin(1) + qJ(2) * t223;
t240 = qJD(2) * t221 + t210 * t198 + V_base(3);
t239 = qJD(2) * t223 + t244;
t238 = pkin(2) * t218 + pkin(6) * t217;
t237 = rSges(3,1) * t218 - rSges(3,2) * t217;
t234 = Icges(3,5) * t218 - Icges(3,6) * t217;
t187 = Icges(3,2) * t218 + t263;
t188 = Icges(3,1) * t217 + t262;
t231 = t187 * t217 - t188 * t218;
t174 = t238 * t223;
t191 = pkin(2) * t217 - pkin(6) * t218;
t230 = V_base(6) * t191 + (-t174 - t200) * t210 + t239;
t173 = t238 * t221;
t229 = V_base(5) * t174 + (t173 - t198) * V_base(6) + t245;
t228 = -t210 * t173 + (-pkin(5) - t191) * V_base(5) + t240;
t227 = (Icges(3,3) * t223 - t221 * t234) * V_base(5) + (Icges(3,3) * t221 + t223 * t234) * V_base(6) + (Icges(3,5) * t217 + Icges(3,6) * t218) * t210;
t127 = pkin(3) * t258 - t221 * t273;
t190 = -qJD(3) * t218 + t210;
t226 = t190 * t127 + t223 * t247 + t228;
t128 = pkin(3) * t253 + t223 * t273;
t185 = -t221 * t248 + V_base(5);
t225 = -qJD(4) * t218 + t185 * t128 + t229;
t140 = -qJ(4) * t218 + t217 * t267;
t224 = t184 * t140 - t221 * t247 + t230;
t201 = rSges(2,1) * t223 - t221 * rSges(2,2);
t199 = -t221 * rSges(2,1) - rSges(2,2) * t223;
t195 = -Icges(2,2) * t221 + t264;
t194 = -Icges(2,2) * t223 - t265;
t193 = Icges(2,5) * t223 - Icges(2,6) * t221;
t192 = -Icges(2,5) * t221 - Icges(2,6) * t223;
t189 = rSges(3,1) * t217 + rSges(3,2) * t218;
t182 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t181 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t180 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t175 = t218 * t246 + t210;
t167 = t246 * t261 + V_base(5);
t166 = qJD(5) * t260 + t184;
t165 = t221 * rSges(3,3) + t223 * t237;
t164 = rSges(3,3) * t223 - t221 * t237;
t159 = -rSges(4,3) * t218 + (rSges(4,1) * t222 - rSges(4,2) * t220) * t217;
t151 = -Icges(4,5) * t218 + (Icges(4,1) * t222 - Icges(4,4) * t220) * t217;
t150 = -Icges(4,6) * t218 + (Icges(4,4) * t222 - Icges(4,2) * t220) * t217;
t148 = t205 * t259 + t257;
t147 = -t204 * t259 + t256;
t146 = t204 * t223 - t218 * t256;
t145 = t205 * t223 + t218 * t257;
t144 = -rSges(5,3) * t218 + (rSges(5,1) * t208 - rSges(5,2) * t207) * t217;
t143 = -Icges(5,5) * t218 + (Icges(5,1) * t208 - Icges(5,4) * t207) * t217;
t142 = -Icges(5,6) * t218 + (Icges(5,4) * t208 - Icges(5,2) * t207) * t217;
t139 = V_base(6) * rSges(2,3) - t201 * t210 + t244;
t138 = t199 * t210 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t137 = -rSges(6,3) * t218 + (rSges(6,1) * t205 - rSges(6,2) * t204) * t217;
t136 = -Icges(6,5) * t218 + (Icges(6,1) * t205 - Icges(6,4) * t204) * t217;
t135 = -Icges(6,6) * t218 + (Icges(6,4) * t205 - Icges(6,2) * t204) * t217;
t134 = -Icges(6,3) * t218 + (Icges(6,5) * t205 - Icges(6,6) * t204) * t217;
t133 = -t199 * V_base(6) + t201 * V_base(5) + V_base(1);
t131 = -pkin(7) * t218 + t217 * t250;
t130 = t172 * rSges(4,1) + t171 * rSges(4,2) + rSges(4,3) * t260;
t129 = rSges(4,1) * t170 + rSges(4,2) * t169 - rSges(4,3) * t261;
t126 = Icges(4,1) * t172 + Icges(4,4) * t171 + Icges(4,5) * t260;
t125 = Icges(4,1) * t170 + Icges(4,4) * t169 - Icges(4,5) * t261;
t124 = Icges(4,4) * t172 + Icges(4,2) * t171 + Icges(4,6) * t260;
t123 = Icges(4,4) * t170 + Icges(4,2) * t169 - Icges(4,6) * t261;
t120 = t163 * rSges(5,1) + t162 * rSges(5,2) + rSges(5,3) * t260;
t119 = rSges(5,1) * t161 + rSges(5,2) * t160 - rSges(5,3) * t261;
t118 = Icges(5,1) * t163 + Icges(5,4) * t162 + Icges(5,5) * t260;
t117 = Icges(5,1) * t161 + Icges(5,4) * t160 - Icges(5,5) * t261;
t116 = Icges(5,4) * t163 + Icges(5,2) * t162 + Icges(5,6) * t260;
t115 = Icges(5,4) * t161 + Icges(5,2) * t160 - Icges(5,6) * t261;
t110 = t148 * rSges(6,1) + t147 * rSges(6,2) + rSges(6,3) * t260;
t109 = rSges(6,1) * t146 + rSges(6,2) * t145 - rSges(6,3) * t261;
t108 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t260;
t107 = Icges(6,1) * t146 + Icges(6,4) * t145 - Icges(6,5) * t261;
t106 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t260;
t105 = Icges(6,4) * t146 + Icges(6,2) * t145 - Icges(6,6) * t261;
t104 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t260;
t103 = Icges(6,5) * t146 + Icges(6,6) * t145 - Icges(6,3) * t261;
t102 = t189 * V_base(6) + (-t165 - t200) * t210 + t239;
t101 = t164 * t210 + (-pkin(5) - t189) * V_base(5) + t240;
t100 = t241 * t221 + t223 * t274;
t99 = -t221 * t274 + t241 * t223;
t98 = t165 * V_base(5) + (-t164 - t198) * V_base(6) + t245;
t97 = -t130 * t190 + t159 * t184 + t230;
t96 = t129 * t190 - t159 * t185 + t228;
t95 = -t129 * t184 + t130 * t185 + t229;
t94 = t144 * t184 + (-t120 - t128) * t190 + t224;
t93 = t119 * t190 + (-t140 - t144) * t185 + t226;
t92 = t120 * t185 + (-t119 - t127) * t184 + t225;
t91 = -t110 * t175 + t131 * t184 + t137 * t166 + (-t100 - t128) * t190 + t224;
t90 = t109 * t175 - t137 * t167 + t190 * t99 + (-t131 - t140) * t185 + t226;
t89 = t100 * t185 - t109 * t166 + t110 * t167 + (-t127 - t99) * t184 + t225;
t1 = m(1) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(2) * (t133 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(3) * (t101 ^ 2 + t102 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t175 * ((-t103 * t167 - t104 * t166 - t134 * t175) * t218 + ((-t135 * t204 + t136 * t205) * t175 + (-t105 * t204 + t107 * t205) * t167 + (-t106 * t204 + t108 * t205) * t166) * t217) / 0.2e1 + t167 * ((-t134 * t261 + t135 * t145 + t136 * t146) * t175 + (-t103 * t261 + t105 * t145 + t107 * t146) * t167 + (-t104 * t261 + t106 * t145 + t108 * t146) * t166) / 0.2e1 + t166 * ((t134 * t260 + t147 * t135 + t148 * t136) * t175 + (t103 * t260 + t147 * t105 + t148 * t107) * t167 + (t104 * t260 + t147 * t106 + t148 * t108) * t166) / 0.2e1 + ((t162 * t142 + t163 * t143 + t171 * t150 + t172 * t151 + t260 * t277) * t190 + (t162 * t115 + t163 * t117 + t171 * t123 + t172 * t125 + t260 * t279) * t185 + (t162 * t116 + t163 * t118 + t171 * t124 + t172 * t126 + t278 * t260) * t184) * t184 / 0.2e1 + ((t142 * t160 + t143 * t161 + t150 * t169 + t151 * t170 - t261 * t277) * t190 + (t115 * t160 + t117 * t161 + t123 * t169 + t125 * t170 - t261 * t279) * t185 + (t116 * t160 + t118 * t161 + t124 * t169 + t126 * t170 - t261 * t278) * t184) * t185 / 0.2e1 + ((-t278 * t184 - t185 * t279 - t277 * t190) * t218 + ((-t142 * t207 + t143 * t208 - t150 * t220 + t151 * t222) * t190 + (-t115 * t207 + t117 * t208 - t123 * t220 + t125 * t222) * t185 + (-t116 * t207 + t118 * t208 - t124 * t220 + t126 * t222) * t184) * t217) * t190 / 0.2e1 + ((t156 * t218 + t158 * t217 + t193) * V_base(6) + (t155 * t218 + t157 * t217 + t192) * V_base(5) + (t187 * t218 + t188 * t217 + Icges(2,3)) * t210) * t210 / 0.2e1 + (t227 * t223 + (t231 * t221 + t192) * t210 + (-t195 * t223 - t221 * t275 + Icges(1,6)) * V_base(6) + (-t194 * t223 - t221 * t276 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t227 * t221 + (-t231 * t223 + t193) * t210 + (-t221 * t195 + t223 * t275 + Icges(1,3)) * V_base(6) + (-t221 * t194 + t223 * t276 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
