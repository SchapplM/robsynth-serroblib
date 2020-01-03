% Calculate kinetic energy for
% S5RRPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:19:58
% EndTime: 2019-12-31 20:20:00
% DurationCPUTime: 2.72s
% Computational Cost: add. (1487->304), mult. (1640->447), div. (0->0), fcn. (1534->10), ass. (0->149)
t274 = Icges(3,3) + Icges(4,3);
t210 = qJ(2) + pkin(9);
t201 = sin(t210);
t202 = cos(t210);
t214 = sin(qJ(2));
t217 = cos(qJ(2));
t273 = Icges(3,5) * t217 + Icges(4,5) * t202 - Icges(3,6) * t214 - Icges(4,6) * t201;
t215 = sin(qJ(1));
t218 = cos(qJ(1));
t258 = Icges(4,4) * t202;
t234 = -Icges(4,2) * t201 + t258;
t141 = -Icges(4,6) * t218 + t215 * t234;
t142 = Icges(4,6) * t215 + t218 * t234;
t259 = Icges(4,4) * t201;
t236 = Icges(4,1) * t202 - t259;
t143 = -Icges(4,5) * t218 + t215 * t236;
t144 = Icges(4,5) * t215 + t218 * t236;
t260 = Icges(3,4) * t217;
t235 = -Icges(3,2) * t214 + t260;
t154 = -Icges(3,6) * t218 + t215 * t235;
t155 = Icges(3,6) * t215 + t218 * t235;
t261 = Icges(3,4) * t214;
t237 = Icges(3,1) * t217 - t261;
t156 = -Icges(3,5) * t218 + t215 * t237;
t157 = Icges(3,5) * t215 + t218 * t237;
t171 = Icges(4,2) * t202 + t259;
t172 = Icges(4,1) * t201 + t258;
t184 = Icges(3,2) * t217 + t261;
t187 = Icges(3,1) * t214 + t260;
t197 = -qJD(2) * t218 + V_base(5);
t198 = qJD(2) * t215 + V_base(4);
t203 = V_base(6) + qJD(1);
t272 = (-t171 * t201 + t172 * t202 - t184 * t214 + t187 * t217) * t203 + (-t142 * t201 + t144 * t202 - t155 * t214 + t157 * t217) * t198 + (-t141 * t201 + t143 * t202 - t154 * t214 + t156 * t217) * t197;
t271 = (Icges(3,5) * t214 + Icges(4,5) * t201 + Icges(3,6) * t217 + Icges(4,6) * t202) * t203 + (t274 * t215 + t273 * t218) * t198 + (t273 * t215 - t274 * t218) * t197;
t267 = pkin(2) * t214;
t266 = pkin(2) * t217;
t216 = cos(qJ(4));
t265 = pkin(4) * t216;
t262 = Icges(2,4) * t215;
t257 = t201 * t215;
t256 = t201 * t218;
t211 = qJ(4) + qJ(5);
t206 = sin(t211);
t255 = t206 * t215;
t254 = t206 * t218;
t207 = cos(t211);
t253 = t207 * t215;
t252 = t207 * t218;
t213 = sin(qJ(4));
t251 = t213 * t215;
t250 = t213 * t218;
t249 = t215 * t216;
t248 = t216 * t218;
t136 = -qJ(3) * t218 + t215 * t266;
t195 = pkin(1) * t215 - pkin(6) * t218;
t247 = -t136 - t195;
t246 = qJD(4) * t201;
t245 = qJD(5) * t201;
t244 = V_base(5) * pkin(5) + V_base(1);
t164 = t218 * t246 + t198;
t241 = qJD(3) * t215 + t197 * t267 + t244;
t240 = pkin(3) * t202 + pkin(7) * t201;
t239 = rSges(3,1) * t217 - rSges(3,2) * t214;
t238 = rSges(4,1) * t202 - rSges(4,2) * t201;
t163 = t215 * t246 + t197;
t196 = pkin(1) * t218 + pkin(6) * t215;
t231 = -V_base(4) * pkin(5) + t203 * t196 + V_base(2);
t230 = V_base(4) * t195 - t196 * V_base(5) + V_base(3);
t229 = t198 * t136 + t230;
t228 = pkin(8) * t201 + t202 * t265;
t137 = qJ(3) * t215 + t218 * t266;
t225 = -qJD(3) * t218 + t203 * t137 + t231;
t160 = t240 * t215;
t174 = t201 * pkin(3) - t202 * pkin(7);
t224 = t197 * t174 + (-t160 + t247) * t203 + t241;
t161 = t240 * t218;
t223 = t198 * t160 + (-t137 - t161) * t197 + t229;
t222 = t203 * t161 + (-t174 - t267) * t198 + t225;
t208 = Icges(2,4) * t218;
t194 = rSges(2,1) * t218 - rSges(2,2) * t215;
t193 = rSges(2,1) * t215 + rSges(2,2) * t218;
t192 = rSges(3,1) * t214 + rSges(3,2) * t217;
t189 = Icges(2,1) * t218 - t262;
t188 = Icges(2,1) * t215 + t208;
t186 = -Icges(2,2) * t215 + t208;
t185 = Icges(2,2) * t218 + t262;
t180 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t179 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t178 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t177 = -qJD(4) * t202 + t203;
t173 = rSges(4,1) * t201 + rSges(4,2) * t202;
t168 = t202 * t248 + t251;
t167 = -t202 * t250 + t249;
t166 = t202 * t249 - t250;
t165 = -t202 * t251 - t248;
t162 = (-qJD(4) - qJD(5)) * t202 + t203;
t159 = rSges(3,3) * t215 + t218 * t239;
t158 = -rSges(3,3) * t218 + t215 * t239;
t151 = t202 * t252 + t255;
t150 = -t202 * t254 + t253;
t149 = t202 * t253 - t254;
t148 = -t202 * t255 - t252;
t146 = rSges(4,3) * t215 + t218 * t238;
t145 = -rSges(4,3) * t218 + t215 * t238;
t135 = -rSges(5,3) * t202 + (rSges(5,1) * t216 - rSges(5,2) * t213) * t201;
t134 = V_base(5) * rSges(2,3) - t193 * t203 + t244;
t133 = t194 * t203 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t132 = -Icges(5,5) * t202 + (Icges(5,1) * t216 - Icges(5,4) * t213) * t201;
t131 = -Icges(5,6) * t202 + (Icges(5,4) * t216 - Icges(5,2) * t213) * t201;
t130 = -Icges(5,3) * t202 + (Icges(5,5) * t216 - Icges(5,6) * t213) * t201;
t129 = t218 * t245 + t164;
t128 = t215 * t245 + t163;
t126 = t193 * V_base(4) - t194 * V_base(5) + V_base(3);
t124 = -rSges(6,3) * t202 + (rSges(6,1) * t207 - rSges(6,2) * t206) * t201;
t123 = -Icges(6,5) * t202 + (Icges(6,1) * t207 - Icges(6,4) * t206) * t201;
t122 = -Icges(6,6) * t202 + (Icges(6,4) * t207 - Icges(6,2) * t206) * t201;
t121 = -Icges(6,3) * t202 + (Icges(6,5) * t207 - Icges(6,6) * t206) * t201;
t119 = -pkin(8) * t202 + t201 * t265;
t118 = rSges(5,1) * t168 + rSges(5,2) * t167 + rSges(5,3) * t256;
t117 = rSges(5,1) * t166 + rSges(5,2) * t165 + rSges(5,3) * t257;
t116 = Icges(5,1) * t168 + Icges(5,4) * t167 + Icges(5,5) * t256;
t115 = Icges(5,1) * t166 + Icges(5,4) * t165 + Icges(5,5) * t257;
t114 = Icges(5,4) * t168 + Icges(5,2) * t167 + Icges(5,6) * t256;
t113 = Icges(5,4) * t166 + Icges(5,2) * t165 + Icges(5,6) * t257;
t112 = Icges(5,5) * t168 + Icges(5,6) * t167 + Icges(5,3) * t256;
t111 = Icges(5,5) * t166 + Icges(5,6) * t165 + Icges(5,3) * t257;
t110 = pkin(4) * t251 + t218 * t228;
t109 = -pkin(4) * t250 + t215 * t228;
t108 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t256;
t107 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t257;
t106 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t256;
t105 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t257;
t104 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t256;
t103 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t257;
t102 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t256;
t101 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t257;
t100 = t192 * t197 + (-t158 - t195) * t203 + t244;
t99 = t159 * t203 - t192 * t198 + t231;
t98 = t158 * t198 - t159 * t197 + t230;
t97 = t173 * t197 + (-t145 + t247) * t203 + t241;
t96 = t146 * t203 + (-t173 - t267) * t198 + t225;
t95 = t145 * t198 + (-t137 - t146) * t197 + t229;
t94 = -t117 * t177 + t135 * t163 + t224;
t93 = t118 * t177 - t135 * t164 + t222;
t92 = t117 * t164 - t118 * t163 + t223;
t91 = -t107 * t162 - t109 * t177 + t119 * t163 + t124 * t128 + t224;
t90 = t108 * t162 + t110 * t177 - t119 * t164 - t124 * t129 + t222;
t89 = t107 * t129 - t108 * t128 + t109 * t164 - t110 * t163 + t223;
t1 = m(1) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + m(2) * (t126 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t164 * ((t112 * t256 + t167 * t114 + t168 * t116) * t164 + (t111 * t256 + t113 * t167 + t115 * t168) * t163 + (t130 * t256 + t131 * t167 + t132 * t168) * t177) / 0.2e1 + t163 * ((t112 * t257 + t114 * t165 + t116 * t166) * t164 + (t111 * t257 + t165 * t113 + t166 * t115) * t163 + (t130 * t257 + t131 * t165 + t132 * t166) * t177) / 0.2e1 + t177 * ((-t111 * t163 - t112 * t164 - t130 * t177) * t202 + ((-t114 * t213 + t116 * t216) * t164 + (-t113 * t213 + t115 * t216) * t163 + (-t131 * t213 + t132 * t216) * t177) * t201) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t129 * ((t102 * t256 + t150 * t104 + t151 * t106) * t129 + (t101 * t256 + t103 * t150 + t105 * t151) * t128 + (t121 * t256 + t122 * t150 + t123 * t151) * t162) / 0.2e1 + t128 * ((t102 * t257 + t104 * t148 + t106 * t149) * t129 + (t101 * t257 + t148 * t103 + t149 * t105) * t128 + (t121 * t257 + t122 * t148 + t123 * t149) * t162) / 0.2e1 + t162 * ((-t101 * t128 - t102 * t129 - t121 * t162) * t202 + ((-t104 * t206 + t106 * t207) * t129 + (-t103 * t206 + t105 * t207) * t128 + (-t122 * t206 + t123 * t207) * t162) * t201) / 0.2e1 + (t272 * t215 - t271 * t218) * t197 / 0.2e1 + (t271 * t215 + t272 * t218) * t198 / 0.2e1 + ((-t185 * t215 + t188 * t218 + Icges(1,4)) * V_base(5) + (-t215 * t186 + t218 * t189 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t185 * t218 + t188 * t215 + Icges(1,2)) * V_base(5) + (t186 * t218 + t189 * t215 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t142 * t202 + t144 * t201 + t155 * t217 + t157 * t214) * t198 + (t141 * t202 + t143 * t201 + t154 * t217 + t156 * t214) * t197 + (t202 * t171 + t201 * t172 + t217 * t184 + t214 * t187 + Icges(2,3)) * t203) * t203 / 0.2e1 + t203 * V_base(4) * (Icges(2,5) * t218 - Icges(2,6) * t215) + V_base(5) * t203 * (Icges(2,5) * t215 + Icges(2,6) * t218) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
