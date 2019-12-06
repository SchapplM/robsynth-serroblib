% Calculate kinetic energy for
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:42
% EndTime: 2019-12-05 17:11:45
% DurationCPUTime: 2.82s
% Computational Cost: add. (1481->331), mult. (2034->509), div. (0->0), fcn. (1984->10), ass. (0->154)
t221 = cos(qJ(3));
t265 = t221 * pkin(3);
t217 = sin(pkin(9));
t263 = Icges(2,4) * t217;
t220 = sin(qJ(2));
t262 = Icges(3,4) * t220;
t222 = cos(qJ(2));
t261 = Icges(3,4) * t222;
t219 = sin(qJ(3));
t260 = t217 * t219;
t259 = t217 * t220;
t258 = t217 * t222;
t218 = cos(pkin(9));
t257 = t218 * t219;
t256 = t218 * t220;
t255 = t218 * t222;
t254 = t219 * t222;
t253 = t221 * t222;
t216 = qJ(3) + qJ(4);
t211 = cos(t216);
t252 = pkin(4) * t211;
t250 = qJD(3) * t220;
t249 = qJD(4) * t220;
t248 = qJD(5) * t220;
t247 = -qJD(3) - qJD(4);
t246 = V_base(5) * qJ(1) + V_base(1);
t242 = qJD(1) + V_base(3);
t199 = qJD(2) * t217 + V_base(4);
t210 = sin(t216);
t241 = pkin(4) * t210;
t168 = t218 * t250 + t199;
t240 = pkin(2) * t222 + pkin(6) * t220;
t198 = -qJD(2) * t218 + V_base(5);
t239 = rSges(3,1) * t222 - rSges(3,2) * t220;
t139 = t218 * t249 + t168;
t238 = Icges(3,1) * t222 - t262;
t237 = -Icges(3,2) * t220 + t261;
t236 = Icges(3,5) * t222 - Icges(3,6) * t220;
t167 = t217 * t250 + t198;
t192 = pkin(1) * t218 + pkin(5) * t217;
t235 = -V_base(4) * qJ(1) + V_base(6) * t192 + V_base(2);
t138 = t217 * t249 + t167;
t191 = pkin(1) * t217 - pkin(5) * t218;
t234 = V_base(4) * t191 - t192 * V_base(5) + t242;
t233 = pkin(7) * t220 + t222 * t265;
t232 = pkin(8) * t220 + t222 * t252;
t173 = t240 * t217;
t197 = t220 * pkin(2) - t222 * pkin(6);
t231 = t198 * t197 + (-t173 - t191) * V_base(6) + t246;
t230 = (-Icges(3,3) * t218 + t217 * t236) * t198 + (Icges(3,3) * t217 + t218 * t236) * t199 + (Icges(3,5) * t220 + Icges(3,6) * t222) * V_base(6);
t174 = t240 * t218;
t229 = V_base(6) * t174 - t197 * t199 + t235;
t228 = t199 * t173 - t174 * t198 + t234;
t123 = -pkin(3) * t257 + t217 * t233;
t137 = -pkin(7) * t222 + t220 * t265;
t200 = -qJD(3) * t222 + V_base(6);
t227 = -t123 * t200 + t167 * t137 + t231;
t124 = pkin(3) * t260 + t218 * t233;
t226 = t200 * t124 - t137 * t168 + t229;
t225 = t168 * t123 - t124 * t167 + t228;
t150 = -Icges(3,6) * t218 + t217 * t237;
t151 = Icges(3,6) * t217 + t218 * t237;
t152 = -Icges(3,5) * t218 + t217 * t238;
t153 = Icges(3,5) * t217 + t218 * t238;
t194 = Icges(3,2) * t222 + t262;
t195 = Icges(3,1) * t220 + t261;
t224 = (-t151 * t220 + t153 * t222) * t199 + (-t150 * t220 + t152 * t222) * t198 + (-t194 * t220 + t195 * t222) * V_base(6);
t213 = qJ(5) + t216;
t209 = Icges(2,4) * t218;
t206 = cos(t213);
t205 = sin(t213);
t196 = rSges(3,1) * t220 + rSges(3,2) * t222;
t190 = rSges(2,1) * t218 - rSges(2,2) * t217;
t189 = rSges(2,1) * t217 + rSges(2,2) * t218;
t188 = Icges(2,1) * t218 - t263;
t187 = Icges(2,1) * t217 + t209;
t186 = -Icges(2,2) * t217 + t209;
t185 = Icges(2,2) * t218 + t263;
t181 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t180 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t179 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t177 = t222 * t247 + V_base(6);
t172 = t218 * t253 + t260;
t171 = t217 * t221 - t218 * t254;
t170 = t217 * t253 - t257;
t169 = -t217 * t254 - t218 * t221;
t166 = V_base(6) + (-qJD(5) + t247) * t222;
t163 = -rSges(4,3) * t222 + (rSges(4,1) * t221 - rSges(4,2) * t219) * t220;
t162 = -Icges(4,5) * t222 + (Icges(4,1) * t221 - Icges(4,4) * t219) * t220;
t161 = -Icges(4,6) * t222 + (Icges(4,4) * t221 - Icges(4,2) * t219) * t220;
t160 = -Icges(4,3) * t222 + (Icges(4,5) * t221 - Icges(4,6) * t219) * t220;
t159 = t210 * t217 + t211 * t255;
t158 = -t210 * t255 + t211 * t217;
t157 = -t210 * t218 + t211 * t258;
t156 = -t210 * t258 - t211 * t218;
t155 = rSges(3,3) * t217 + t218 * t239;
t154 = -rSges(3,3) * t218 + t217 * t239;
t147 = t205 * t217 + t206 * t255;
t146 = -t205 * t255 + t206 * t217;
t145 = -t205 * t218 + t206 * t258;
t144 = -t205 * t258 - t206 * t218;
t143 = -rSges(5,3) * t222 + (rSges(5,1) * t211 - rSges(5,2) * t210) * t220;
t142 = -Icges(5,5) * t222 + (Icges(5,1) * t211 - Icges(5,4) * t210) * t220;
t141 = -Icges(5,6) * t222 + (Icges(5,4) * t211 - Icges(5,2) * t210) * t220;
t140 = -Icges(5,3) * t222 + (Icges(5,5) * t211 - Icges(5,6) * t210) * t220;
t135 = V_base(5) * rSges(2,3) - t189 * V_base(6) + t246;
t134 = t190 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t133 = -rSges(6,3) * t222 + (rSges(6,1) * t206 - rSges(6,2) * t205) * t220;
t132 = -Icges(6,5) * t222 + (Icges(6,1) * t206 - Icges(6,4) * t205) * t220;
t131 = -Icges(6,6) * t222 + (Icges(6,4) * t206 - Icges(6,2) * t205) * t220;
t130 = -Icges(6,3) * t222 + (Icges(6,5) * t206 - Icges(6,6) * t205) * t220;
t129 = t189 * V_base(4) - t190 * V_base(5) + t242;
t128 = t218 * t248 + t139;
t127 = t217 * t248 + t138;
t126 = -pkin(8) * t222 + t220 * t252;
t122 = rSges(4,1) * t172 + rSges(4,2) * t171 + rSges(4,3) * t256;
t121 = rSges(4,1) * t170 + rSges(4,2) * t169 + rSges(4,3) * t259;
t120 = Icges(4,1) * t172 + Icges(4,4) * t171 + Icges(4,5) * t256;
t119 = Icges(4,1) * t170 + Icges(4,4) * t169 + Icges(4,5) * t259;
t118 = Icges(4,4) * t172 + Icges(4,2) * t171 + Icges(4,6) * t256;
t117 = Icges(4,4) * t170 + Icges(4,2) * t169 + Icges(4,6) * t259;
t116 = Icges(4,5) * t172 + Icges(4,6) * t171 + Icges(4,3) * t256;
t115 = Icges(4,5) * t170 + Icges(4,6) * t169 + Icges(4,3) * t259;
t114 = rSges(5,1) * t159 + rSges(5,2) * t158 + rSges(5,3) * t256;
t113 = rSges(5,1) * t157 + rSges(5,2) * t156 + rSges(5,3) * t259;
t112 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t256;
t111 = Icges(5,1) * t157 + Icges(5,4) * t156 + Icges(5,5) * t259;
t110 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t256;
t109 = Icges(5,4) * t157 + Icges(5,2) * t156 + Icges(5,6) * t259;
t108 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t256;
t107 = Icges(5,5) * t157 + Icges(5,6) * t156 + Icges(5,3) * t259;
t105 = rSges(6,1) * t147 + rSges(6,2) * t146 + rSges(6,3) * t256;
t104 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t259;
t103 = Icges(6,1) * t147 + Icges(6,4) * t146 + Icges(6,5) * t256;
t102 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t259;
t101 = Icges(6,4) * t147 + Icges(6,2) * t146 + Icges(6,6) * t256;
t100 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t259;
t99 = Icges(6,5) * t147 + Icges(6,6) * t146 + Icges(6,3) * t256;
t98 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t259;
t96 = t196 * t198 + (-t154 - t191) * V_base(6) + t246;
t95 = t155 * V_base(6) - t196 * t199 + t235;
t94 = t217 * t241 + t218 * t232;
t93 = t217 * t232 - t218 * t241;
t92 = t154 * t199 - t155 * t198 + t234;
t91 = -t121 * t200 + t163 * t167 + t231;
t90 = t122 * t200 - t163 * t168 + t229;
t89 = t121 * t168 - t122 * t167 + t228;
t88 = -t113 * t177 + t138 * t143 + t227;
t87 = t114 * t177 - t139 * t143 + t226;
t86 = t113 * t139 - t114 * t138 + t225;
t85 = -t104 * t166 + t126 * t138 + t127 * t133 - t177 * t93 + t227;
t84 = t105 * t166 - t126 * t139 - t128 * t133 + t177 * t94 + t226;
t83 = t104 * t128 - t105 * t127 - t138 * t94 + t139 * t93 + t225;
t1 = m(1) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(2) * (t129 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t199 * (t230 * t217 + t224 * t218) / 0.2e1 + t198 * (t224 * t217 - t230 * t218) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t168 * ((t116 * t256 + t171 * t118 + t172 * t120) * t168 + (t115 * t256 + t117 * t171 + t119 * t172) * t167 + (t160 * t256 + t161 * t171 + t162 * t172) * t200) / 0.2e1 + t167 * ((t116 * t259 + t118 * t169 + t120 * t170) * t168 + (t115 * t259 + t169 * t117 + t170 * t119) * t167 + (t160 * t259 + t161 * t169 + t162 * t170) * t200) / 0.2e1 + t200 * ((-t115 * t167 - t116 * t168 - t160 * t200) * t222 + ((-t118 * t219 + t120 * t221) * t168 + (-t117 * t219 + t119 * t221) * t167 + (-t161 * t219 + t162 * t221) * t200) * t220) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + t139 * ((t108 * t256 + t158 * t110 + t159 * t112) * t139 + (t107 * t256 + t109 * t158 + t111 * t159) * t138 + (t140 * t256 + t141 * t158 + t142 * t159) * t177) / 0.2e1 + t138 * ((t108 * t259 + t110 * t156 + t112 * t157) * t139 + (t107 * t259 + t156 * t109 + t157 * t111) * t138 + (t140 * t259 + t141 * t156 + t142 * t157) * t177) / 0.2e1 + t177 * ((-t107 * t138 - t108 * t139 - t140 * t177) * t222 + ((-t110 * t210 + t112 * t211) * t139 + (-t109 * t210 + t111 * t211) * t138 + (-t141 * t210 + t142 * t211) * t177) * t220) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t128 * ((t146 * t101 + t147 * t103 + t99 * t256) * t128 + (t100 * t146 + t102 * t147 + t256 * t98) * t127 + (t130 * t256 + t131 * t146 + t132 * t147) * t166) / 0.2e1 + t127 * ((t101 * t144 + t103 * t145 + t259 * t99) * t128 + (t144 * t100 + t145 * t102 + t98 * t259) * t127 + (t130 * t259 + t131 * t144 + t132 * t145) * t166) / 0.2e1 + t166 * ((-t98 * t127 - t99 * t128 - t130 * t166) * t222 + ((-t101 * t205 + t103 * t206) * t128 + (-t100 * t205 + t102 * t206) * t127 + (-t131 * t205 + t132 * t206) * t166) * t220) / 0.2e1 + ((-t185 * t217 + t187 * t218 + Icges(1,4)) * V_base(5) + (-t186 * t217 + t188 * t218 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t185 * t218 + t187 * t217 + Icges(1,2)) * V_base(5) + (t186 * t218 + t188 * t217 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t151 * t222 + t153 * t220) * t199 + (t150 * t222 + t152 * t220) * t198 + (t222 * t194 + t195 * t220 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t218 - Icges(2,6) * t217 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t217 + Icges(2,6) * t218 + Icges(1,6));
T = t1;
