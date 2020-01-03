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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:41:19
% EndTime: 2020-01-03 11:41:23
% DurationCPUTime: 3.89s
% Computational Cost: add. (1426->329), mult. (1910->466), div. (0->0), fcn. (1860->10), ass. (0->152)
t280 = -Icges(5,3) - Icges(4,3);
t217 = qJ(3) + pkin(9);
t210 = cos(t217);
t219 = cos(pkin(8));
t224 = cos(qJ(1));
t209 = sin(t217);
t222 = sin(qJ(1));
t256 = t222 * t209;
t162 = -t210 * t224 - t219 * t256;
t255 = t222 * t210;
t163 = -t209 * t224 + t219 * t255;
t223 = cos(qJ(3));
t252 = t223 * t224;
t221 = sin(qJ(3));
t254 = t222 * t221;
t171 = -t219 * t254 - t252;
t253 = t222 * t223;
t259 = t221 * t224;
t172 = t219 * t253 - t259;
t218 = sin(pkin(8));
t262 = t218 * t222;
t279 = Icges(4,5) * t172 + Icges(5,5) * t163 + Icges(4,6) * t171 + Icges(5,6) * t162 - t280 * t262;
t260 = t219 * t224;
t164 = t209 * t260 - t255;
t165 = -t210 * t260 - t256;
t173 = t219 * t259 - t253;
t174 = -t219 * t252 - t254;
t261 = t218 * t224;
t278 = Icges(4,5) * t174 + Icges(5,5) * t165 + Icges(4,6) * t173 + Icges(5,6) * t164 + t280 * t261;
t277 = t280 * t219 + (Icges(4,5) * t223 + Icges(5,5) * t210 - Icges(4,6) * t221 - Icges(5,6) * t209) * t218;
t263 = Icges(3,4) * t219;
t238 = -Icges(3,2) * t218 + t263;
t157 = -Icges(3,6) * t224 + t222 * t238;
t264 = Icges(3,4) * t218;
t239 = Icges(3,1) * t219 - t264;
t159 = -Icges(3,5) * t224 + t222 * t239;
t265 = Icges(2,4) * t224;
t276 = Icges(2,1) * t222 - t157 * t218 + t159 * t219 + t265;
t158 = -Icges(3,6) * t222 - t224 * t238;
t160 = -Icges(3,5) * t222 - t224 * t239;
t213 = Icges(2,4) * t222;
t275 = -Icges(2,1) * t224 - t158 * t218 + t160 * t219 + t213;
t251 = pkin(4) * t210;
t274 = pkin(7) * t218 + t219 * t251;
t267 = t223 * pkin(3);
t273 = qJ(4) * t218 + t219 * t267;
t211 = qJ(5) + t217;
t206 = sin(t211);
t258 = t222 * t206;
t207 = cos(t211);
t257 = t222 * t207;
t249 = qJD(3) * t218;
t248 = qJD(4) * t218;
t247 = -qJD(3) - qJD(5);
t202 = -pkin(1) * t224 - t222 * qJ(2);
t246 = V_base(5) * t202 + V_base(1);
t245 = V_base(6) * pkin(5) + V_base(2);
t187 = t222 * t249 + V_base(5);
t212 = V_base(4) + qJD(1);
t242 = pkin(4) * t209;
t241 = pkin(2) * t219 + pkin(6) * t218;
t240 = rSges(3,1) * t219 - rSges(3,2) * t218;
t237 = Icges(3,5) * t219 - Icges(3,6) * t218;
t189 = Icges(3,2) * t219 + t264;
t190 = Icges(3,1) * t218 + t263;
t234 = t189 * t218 - t190 * t219;
t200 = t222 * pkin(1) - qJ(2) * t224;
t233 = -qJD(2) * t222 + t212 * t200 + V_base(3);
t232 = -qJD(2) * t224 + t245;
t175 = t241 * t222;
t176 = t241 * t224;
t231 = -V_base(5) * t176 + (-t175 - t200) * V_base(6) + t246;
t230 = -(-Icges(3,3) * t224 + t222 * t237) * V_base(5) - (-Icges(3,3) * t222 - t224 * t237) * V_base(6) - (Icges(3,5) * t218 + Icges(3,6) * t219) * t212;
t193 = pkin(2) * t218 - pkin(6) * t219;
t229 = V_base(6) * t193 + (t176 - t202) * t212 + t232;
t228 = t212 * t175 + (-pkin(5) - t193) * V_base(5) + t233;
t142 = -qJ(4) * t219 + t218 * t267;
t186 = -t224 * t249 + V_base(6);
t227 = t186 * t142 + t222 * t248 + t229;
t130 = -pkin(3) * t254 - t273 * t224;
t226 = -qJD(4) * t219 + t187 * t130 + t231;
t129 = -pkin(3) * t259 + t273 * t222;
t192 = -qJD(3) * t219 + t212;
t225 = t192 * t129 - t224 * t248 + t228;
t203 = -rSges(2,1) * t224 + t222 * rSges(2,2);
t201 = t222 * rSges(2,1) + rSges(2,2) * t224;
t197 = Icges(2,2) * t222 - t265;
t196 = Icges(2,2) * t224 + t213;
t195 = -Icges(2,5) * t224 + Icges(2,6) * t222;
t194 = Icges(2,5) * t222 + Icges(2,6) * t224;
t191 = rSges(3,1) * t218 + rSges(3,2) * t219;
t184 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t183 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t182 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t177 = t219 * t247 + t212;
t169 = qJD(5) * t262 + t187;
t168 = t247 * t261 + V_base(6);
t167 = -t222 * rSges(3,3) - t224 * t240;
t166 = -rSges(3,3) * t224 + t222 * t240;
t161 = -rSges(4,3) * t219 + (rSges(4,1) * t223 - rSges(4,2) * t221) * t218;
t153 = -Icges(4,5) * t219 + (Icges(4,1) * t223 - Icges(4,4) * t221) * t218;
t152 = -Icges(4,6) * t219 + (Icges(4,4) * t223 - Icges(4,2) * t221) * t218;
t150 = -t207 * t260 - t258;
t149 = t206 * t260 - t257;
t148 = -t206 * t224 + t219 * t257;
t147 = -t207 * t224 - t219 * t258;
t146 = -rSges(5,3) * t219 + (rSges(5,1) * t210 - rSges(5,2) * t209) * t218;
t145 = -Icges(5,5) * t219 + (Icges(5,1) * t210 - Icges(5,4) * t209) * t218;
t144 = -Icges(5,6) * t219 + (Icges(5,4) * t210 - Icges(5,2) * t209) * t218;
t141 = V_base(6) * rSges(2,3) - t203 * t212 + t245;
t140 = t201 * t212 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t139 = -rSges(6,3) * t219 + (rSges(6,1) * t207 - rSges(6,2) * t206) * t218;
t138 = -Icges(6,5) * t219 + (Icges(6,1) * t207 - Icges(6,4) * t206) * t218;
t137 = -Icges(6,6) * t219 + (Icges(6,4) * t207 - Icges(6,2) * t206) * t218;
t136 = -Icges(6,3) * t219 + (Icges(6,5) * t207 - Icges(6,6) * t206) * t218;
t135 = -t201 * V_base(6) + t203 * V_base(5) + V_base(1);
t133 = -pkin(7) * t219 + t218 * t251;
t132 = t174 * rSges(4,1) + t173 * rSges(4,2) - rSges(4,3) * t261;
t131 = rSges(4,1) * t172 + rSges(4,2) * t171 + rSges(4,3) * t262;
t128 = Icges(4,1) * t174 + Icges(4,4) * t173 - Icges(4,5) * t261;
t127 = Icges(4,1) * t172 + Icges(4,4) * t171 + Icges(4,5) * t262;
t126 = Icges(4,4) * t174 + Icges(4,2) * t173 - Icges(4,6) * t261;
t125 = Icges(4,4) * t172 + Icges(4,2) * t171 + Icges(4,6) * t262;
t122 = t165 * rSges(5,1) + t164 * rSges(5,2) - rSges(5,3) * t261;
t121 = rSges(5,1) * t163 + rSges(5,2) * t162 + rSges(5,3) * t262;
t120 = Icges(5,1) * t165 + Icges(5,4) * t164 - Icges(5,5) * t261;
t119 = Icges(5,1) * t163 + Icges(5,4) * t162 + Icges(5,5) * t262;
t118 = Icges(5,4) * t165 + Icges(5,2) * t164 - Icges(5,6) * t261;
t117 = Icges(5,4) * t163 + Icges(5,2) * t162 + Icges(5,6) * t262;
t112 = t150 * rSges(6,1) + t149 * rSges(6,2) - rSges(6,3) * t261;
t111 = rSges(6,1) * t148 + rSges(6,2) * t147 + rSges(6,3) * t262;
t110 = Icges(6,1) * t150 + Icges(6,4) * t149 - Icges(6,5) * t261;
t109 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t262;
t108 = Icges(6,4) * t150 + Icges(6,2) * t149 - Icges(6,6) * t261;
t107 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t262;
t106 = Icges(6,5) * t150 + Icges(6,6) * t149 - Icges(6,3) * t261;
t105 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t262;
t104 = V_base(6) * t191 + (-t167 - t202) * t212 + t232;
t103 = t166 * t212 + (-pkin(5) - t191) * V_base(5) + t233;
t102 = -t242 * t222 - t274 * t224;
t101 = t274 * t222 - t242 * t224;
t100 = t167 * V_base(5) + (-t166 - t200) * V_base(6) + t246;
t99 = -t192 * t132 + t186 * t161 + t229;
t98 = t131 * t192 - t161 * t187 + t228;
t97 = -t131 * t186 + t132 * t187 + t231;
t96 = t186 * t146 + (-t122 - t130) * t192 + t227;
t95 = t192 * t121 + (-t142 - t146) * t187 + t225;
t94 = t122 * t187 + (-t121 - t129) * t186 + t226;
t93 = -t177 * t112 + t186 * t133 + t168 * t139 + (-t102 - t130) * t192 + t227;
t92 = t192 * t101 + t177 * t111 - t169 * t139 + (-t133 - t142) * t187 + t225;
t91 = t102 * t187 - t111 * t168 + t112 * t169 + (-t101 - t129) * t186 + t226;
t1 = m(1) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(2) * (t135 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t177 * ((-t105 * t169 - t106 * t168 - t136 * t177) * t219 + ((-t137 * t206 + t138 * t207) * t177 + (-t107 * t206 + t109 * t207) * t169 + (-t108 * t206 + t110 * t207) * t168) * t218) / 0.2e1 + t169 * ((t136 * t262 + t137 * t147 + t138 * t148) * t177 + (t105 * t262 + t147 * t107 + t148 * t109) * t169 + (t106 * t262 + t108 * t147 + t110 * t148) * t168) / 0.2e1 + t168 * ((-t136 * t261 + t149 * t137 + t150 * t138) * t177 + (-t105 * t261 + t149 * t107 + t150 * t109) * t169 + (-t106 * t261 + t149 * t108 + t150 * t110) * t168) / 0.2e1 + ((t164 * t144 + t165 * t145 + t173 * t152 + t174 * t153 - t277 * t261) * t192 + (t164 * t117 + t165 * t119 + t173 * t125 + t174 * t127 - t279 * t261) * t187 + (t164 * t118 + t165 * t120 + t173 * t126 + t174 * t128 - t278 * t261) * t186) * t186 / 0.2e1 + ((t144 * t162 + t145 * t163 + t152 * t171 + t153 * t172 + t277 * t262) * t192 + (t162 * t117 + t163 * t119 + t171 * t125 + t172 * t127 + t279 * t262) * t187 + (t118 * t162 + t120 * t163 + t126 * t171 + t128 * t172 + t278 * t262) * t186) * t187 / 0.2e1 + ((-t278 * t186 - t279 * t187 - t277 * t192) * t219 + ((-t144 * t209 + t145 * t210 - t152 * t221 + t153 * t223) * t192 + (-t117 * t209 + t119 * t210 - t125 * t221 + t127 * t223) * t187 + (-t118 * t209 + t120 * t210 - t126 * t221 + t128 * t223) * t186) * t218) * t192 / 0.2e1 + ((t158 * t219 + t160 * t218 + t195) * V_base(6) + (t157 * t219 + t159 * t218 + t194) * V_base(5) + (t219 * t189 + t218 * t190 + Icges(2,3)) * t212) * t212 / 0.2e1 + (t230 * t224 + (-t234 * t222 + t194) * t212 + (t197 * t224 + t275 * t222 + Icges(1,6)) * V_base(6) + (t224 * t196 + t276 * t222 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t230 * t222 + (t234 * t224 + t195) * t212 + (t222 * t197 - t275 * t224 + Icges(1,3)) * V_base(6) + (t222 * t196 - t276 * t224 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
