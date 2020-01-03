% Calculate kinetic energy for
% S5RRRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:02
% EndTime: 2019-12-31 21:16:05
% DurationCPUTime: 2.33s
% Computational Cost: add. (1457->309), mult. (1622->460), div. (0->0), fcn. (1516->10), ass. (0->149)
t204 = sin(qJ(2));
t255 = pkin(2) * t204;
t206 = cos(qJ(2));
t254 = pkin(2) * t206;
t202 = cos(pkin(9));
t253 = pkin(4) * t202;
t205 = sin(qJ(1));
t251 = Icges(2,4) * t205;
t250 = Icges(3,4) * t204;
t249 = Icges(3,4) * t206;
t200 = qJ(2) + qJ(3);
t195 = sin(t200);
t248 = Icges(4,4) * t195;
t196 = cos(t200);
t247 = Icges(4,4) * t196;
t246 = t195 * t205;
t207 = cos(qJ(1));
t245 = t195 * t207;
t244 = t196 * t205;
t243 = t196 * t207;
t201 = sin(pkin(9));
t242 = t201 * t205;
t241 = t201 * t207;
t240 = t202 * t205;
t239 = t202 * t207;
t125 = -pkin(7) * t207 + t205 * t254;
t184 = t205 * pkin(1) - t207 * pkin(6);
t237 = -t125 - t184;
t236 = qJD(4) * t195;
t235 = qJD(5) * t195;
t234 = V_base(5) * pkin(5) + V_base(1);
t227 = pkin(3) * t196 + qJ(4) * t195;
t150 = t227 * t205;
t231 = -t150 + t237;
t187 = qJD(2) * t205 + V_base(4);
t192 = V_base(6) + qJD(1);
t186 = -qJD(2) * t207 + V_base(5);
t230 = t186 * t255 + t234;
t163 = qJD(3) * t205 + t187;
t229 = rSges(3,1) * t206 - rSges(3,2) * t204;
t228 = rSges(4,1) * t196 - rSges(4,2) * t195;
t226 = Icges(3,1) * t206 - t250;
t225 = Icges(4,1) * t196 - t248;
t224 = -Icges(3,2) * t204 + t249;
t223 = -Icges(4,2) * t195 + t247;
t222 = Icges(3,5) * t206 - Icges(3,6) * t204;
t221 = Icges(4,5) * t196 - Icges(4,6) * t195;
t160 = pkin(3) * t195 - qJ(4) * t196;
t162 = V_base(5) + (-qJD(2) - qJD(3)) * t207;
t220 = t162 * t160 + t207 * t236 + t230;
t185 = t207 * pkin(1) + t205 * pkin(6);
t219 = -V_base(4) * pkin(5) + t192 * t185 + V_base(2);
t218 = V_base(4) * t184 - t185 * V_base(5) + V_base(3);
t217 = (-Icges(4,3) * t207 + t205 * t221) * t162 + (Icges(4,3) * t205 + t207 * t221) * t163 + (Icges(4,5) * t195 + Icges(4,6) * t196) * t192;
t216 = (-Icges(3,3) * t207 + t205 * t222) * t186 + (Icges(3,3) * t205 + t207 * t222) * t187 + (Icges(3,5) * t204 + Icges(3,6) * t206) * t192;
t215 = pkin(8) * t195 + t196 * t253;
t126 = pkin(7) * t205 + t207 * t254;
t214 = t187 * t125 - t126 * t186 + t218;
t213 = t192 * t126 - t187 * t255 + t219;
t151 = t227 * t207;
t212 = t192 * t151 + t205 * t236 + t213;
t211 = -qJD(4) * t196 + t163 * t150 + t214;
t130 = -Icges(4,6) * t207 + t205 * t223;
t131 = Icges(4,6) * t205 + t207 * t223;
t132 = -Icges(4,5) * t207 + t205 * t225;
t133 = Icges(4,5) * t205 + t207 * t225;
t158 = Icges(4,2) * t196 + t248;
t159 = Icges(4,1) * t195 + t247;
t210 = (-t131 * t195 + t133 * t196) * t163 + (-t130 * t195 + t132 * t196) * t162 + (-t158 * t195 + t159 * t196) * t192;
t144 = -Icges(3,6) * t207 + t205 * t224;
t145 = Icges(3,6) * t205 + t207 * t224;
t146 = -Icges(3,5) * t207 + t205 * t226;
t147 = Icges(3,5) * t205 + t207 * t226;
t173 = Icges(3,2) * t206 + t250;
t176 = Icges(3,1) * t204 + t249;
t209 = (-t145 * t204 + t147 * t206) * t187 + (-t144 * t204 + t146 * t206) * t186 + (-t173 * t204 + t176 * t206) * t192;
t199 = pkin(9) + qJ(5);
t197 = Icges(2,4) * t207;
t191 = cos(t199);
t190 = sin(t199);
t181 = rSges(2,1) * t207 - rSges(2,2) * t205;
t180 = rSges(2,1) * t205 + rSges(2,2) * t207;
t179 = rSges(3,1) * t204 + rSges(3,2) * t206;
t178 = Icges(2,1) * t207 - t251;
t177 = Icges(2,1) * t205 + t197;
t175 = -Icges(2,2) * t205 + t197;
t174 = Icges(2,2) * t207 + t251;
t169 = -qJD(5) * t196 + t192;
t168 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t167 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t166 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = rSges(4,1) * t195 + rSges(4,2) * t196;
t155 = t196 * t239 + t242;
t154 = -t196 * t241 + t240;
t153 = t196 * t240 - t241;
t152 = -t196 * t242 - t239;
t149 = rSges(3,3) * t205 + t207 * t229;
t148 = -rSges(3,3) * t207 + t205 * t229;
t141 = t190 * t205 + t191 * t243;
t140 = -t190 * t243 + t191 * t205;
t139 = -t190 * t207 + t191 * t244;
t138 = -t190 * t244 - t191 * t207;
t137 = t207 * t235 + t163;
t136 = t205 * t235 + t162;
t135 = rSges(4,3) * t205 + t207 * t228;
t134 = -rSges(4,3) * t207 + t205 * t228;
t124 = -rSges(5,3) * t196 + (rSges(5,1) * t202 - rSges(5,2) * t201) * t195;
t123 = -Icges(5,5) * t196 + (Icges(5,1) * t202 - Icges(5,4) * t201) * t195;
t122 = -Icges(5,6) * t196 + (Icges(5,4) * t202 - Icges(5,2) * t201) * t195;
t121 = -Icges(5,3) * t196 + (Icges(5,5) * t202 - Icges(5,6) * t201) * t195;
t120 = V_base(5) * rSges(2,3) - t180 * t192 + t234;
t119 = t181 * t192 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t117 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t115 = -rSges(6,3) * t196 + (rSges(6,1) * t191 - rSges(6,2) * t190) * t195;
t114 = -Icges(6,5) * t196 + (Icges(6,1) * t191 - Icges(6,4) * t190) * t195;
t113 = -Icges(6,6) * t196 + (Icges(6,4) * t191 - Icges(6,2) * t190) * t195;
t112 = -Icges(6,3) * t196 + (Icges(6,5) * t191 - Icges(6,6) * t190) * t195;
t109 = -pkin(8) * t196 + t195 * t253;
t108 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t245;
t107 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t246;
t106 = Icges(5,1) * t155 + Icges(5,4) * t154 + Icges(5,5) * t245;
t105 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t246;
t104 = Icges(5,4) * t155 + Icges(5,2) * t154 + Icges(5,6) * t245;
t103 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t246;
t102 = Icges(5,5) * t155 + Icges(5,6) * t154 + Icges(5,3) * t245;
t101 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t246;
t100 = pkin(4) * t242 + t207 * t215;
t99 = -pkin(4) * t241 + t205 * t215;
t98 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t245;
t97 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t246;
t96 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t245;
t95 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t246;
t94 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t245;
t93 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t246;
t92 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t245;
t91 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t246;
t90 = t179 * t186 + (-t148 - t184) * t192 + t234;
t89 = t149 * t192 - t179 * t187 + t219;
t88 = t148 * t187 - t149 * t186 + t218;
t87 = t161 * t162 + (-t134 + t237) * t192 + t230;
t86 = t135 * t192 - t161 * t163 + t213;
t85 = t134 * t163 - t135 * t162 + t214;
t84 = t124 * t162 + (-t107 + t231) * t192 + t220;
t83 = t108 * t192 + (-t124 - t160) * t163 + t212;
t82 = t107 * t163 + (-t108 - t151) * t162 + t211;
t81 = t109 * t162 + t115 * t136 - t169 * t97 + (-t99 + t231) * t192 + t220;
t80 = t100 * t192 - t115 * t137 + t169 * t98 + (-t109 - t160) * t163 + t212;
t79 = -t136 * t98 + t137 * t97 + t163 * t99 + (-t100 - t151) * t162 + t211;
t1 = m(1) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(3) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + t187 * (t205 * t216 + t207 * t209) / 0.2e1 + t186 * (t205 * t209 - t207 * t216) / 0.2e1 + m(4) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(5) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + t137 * ((t140 * t94 + t141 * t96 + t92 * t245) * t137 + (t140 * t93 + t141 * t95 + t245 * t91) * t136 + (t112 * t245 + t113 * t140 + t114 * t141) * t169) / 0.2e1 + t136 * ((t138 * t94 + t139 * t96 + t246 * t92) * t137 + (t138 * t93 + t139 * t95 + t91 * t246) * t136 + (t112 * t246 + t113 * t138 + t114 * t139) * t169) / 0.2e1 + t169 * ((-t112 * t169 - t91 * t136 - t92 * t137) * t196 + ((-t190 * t94 + t191 * t96) * t137 + (-t190 * t93 + t191 * t95) * t136 + (-t113 * t190 + t114 * t191) * t169) * t195) / 0.2e1 + (t205 * t210 - t207 * t217 + (t102 * t246 + t104 * t152 + t106 * t153) * t163 + (t101 * t246 + t103 * t152 + t105 * t153) * t162 + (t121 * t246 + t122 * t152 + t123 * t153) * t192) * t162 / 0.2e1 + (t205 * t217 + t210 * t207 + (t102 * t245 + t104 * t154 + t106 * t155) * t163 + (t101 * t245 + t103 * t154 + t105 * t155) * t162 + (t121 * t245 + t122 * t154 + t123 * t155) * t192) * t163 / 0.2e1 + ((-t174 * t205 + t177 * t207 + Icges(1,4)) * V_base(5) + (-t175 * t205 + t178 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t174 * t207 + t177 * t205 + Icges(1,2)) * V_base(5) + (t175 * t207 + t178 * t205 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t145 * t206 + t147 * t204) * t187 + (t144 * t206 + t146 * t204) * t186 + (t131 * t196 + t133 * t195) * t163 + (t130 * t196 + t132 * t195) * t162 + (-t101 * t162 - t102 * t163) * t196 + ((-t104 * t201 + t106 * t202) * t163 + (-t103 * t201 + t105 * t202) * t162) * t195 + (t173 * t206 + t176 * t204 + Icges(2,3) + (t158 - t121) * t196 + (-t122 * t201 + t123 * t202 + t159) * t195) * t192) * t192 / 0.2e1 + t192 * V_base(4) * (Icges(2,5) * t207 - Icges(2,6) * t205) + V_base(5) * t192 * (Icges(2,5) * t205 + Icges(2,6) * t207) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
