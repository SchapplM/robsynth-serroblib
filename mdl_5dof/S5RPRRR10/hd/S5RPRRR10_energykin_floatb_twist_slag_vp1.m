% Calculate kinetic energy for
% S5RPRRR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:39
% EndTime: 2019-12-31 19:09:42
% DurationCPUTime: 2.40s
% Computational Cost: add. (1465->310), mult. (1618->461), div. (0->0), fcn. (1512->10), ass. (0->154)
t212 = sin(pkin(9));
t268 = pkin(2) * t212;
t213 = cos(pkin(9));
t267 = pkin(2) * t213;
t217 = cos(qJ(4));
t266 = pkin(4) * t217;
t216 = sin(qJ(1));
t264 = Icges(2,4) * t216;
t263 = Icges(3,4) * t212;
t262 = Icges(3,4) * t213;
t210 = pkin(9) + qJ(3);
t201 = sin(t210);
t261 = Icges(4,4) * t201;
t202 = cos(t210);
t260 = Icges(4,4) * t202;
t259 = t201 * t216;
t218 = cos(qJ(1));
t258 = t201 * t218;
t211 = qJ(4) + qJ(5);
t206 = sin(t211);
t257 = t206 * t216;
t256 = t206 * t218;
t207 = cos(t211);
t255 = t207 * t216;
t254 = t207 * t218;
t215 = sin(qJ(4));
t253 = t215 * t216;
t252 = t215 * t218;
t251 = t216 * t217;
t250 = t217 * t218;
t136 = -pkin(6) * t218 + t216 * t267;
t192 = pkin(1) * t216 - qJ(2) * t218;
t248 = -t136 - t192;
t247 = qJD(4) * t201;
t246 = qJD(5) * t201;
t245 = V_base(4) * t192 + V_base(3);
t244 = V_base(5) * pkin(5) + V_base(1);
t197 = qJD(3) * t216 + V_base(4);
t203 = V_base(6) + qJD(1);
t241 = qJD(2) * t216 + t244;
t164 = t218 * t247 + t197;
t240 = V_base(5) * t268 + t241;
t239 = pkin(3) * t202 + pkin(7) * t201;
t196 = -qJD(3) * t218 + V_base(5);
t238 = rSges(3,1) * t213 - rSges(3,2) * t212;
t237 = rSges(4,1) * t202 - rSges(4,2) * t201;
t236 = Icges(3,1) * t213 - t263;
t235 = Icges(4,1) * t202 - t261;
t234 = -Icges(3,2) * t212 + t262;
t233 = -Icges(4,2) * t201 + t260;
t232 = Icges(3,5) * t213 - Icges(3,6) * t212;
t231 = Icges(4,5) * t202 - Icges(4,6) * t201;
t194 = pkin(1) * t218 + qJ(2) * t216;
t230 = -qJD(2) * t218 + t203 * t194 + V_base(2);
t163 = t216 * t247 + t196;
t229 = pkin(8) * t201 + t202 * t266;
t228 = (-Icges(4,3) * t218 + t216 * t231) * t196 + (Icges(4,3) * t216 + t218 * t231) * t197 + (Icges(4,5) * t201 + Icges(4,6) * t202) * t203;
t137 = pkin(6) * t216 + t218 * t267;
t227 = V_base(4) * t136 + (-t137 - t194) * V_base(5) + t245;
t226 = (-Icges(3,3) * t218 + t216 * t232) * V_base(5) + (Icges(3,3) * t216 + t218 * t232) * V_base(4) + (Icges(3,5) * t212 + Icges(3,6) * t213) * t203;
t160 = t239 * t216;
t174 = t201 * pkin(3) - t202 * pkin(7);
t225 = t196 * t174 + (-t160 + t248) * t203 + t240;
t161 = t239 * t218;
t224 = t197 * t160 - t161 * t196 + t227;
t223 = t203 * t137 + (-pkin(5) - t268) * V_base(4) + t230;
t222 = t203 * t161 - t174 * t197 + t223;
t141 = -Icges(4,6) * t218 + t216 * t233;
t142 = Icges(4,6) * t216 + t218 * t233;
t143 = -Icges(4,5) * t218 + t216 * t235;
t144 = Icges(4,5) * t216 + t218 * t235;
t171 = Icges(4,2) * t202 + t261;
t172 = Icges(4,1) * t201 + t260;
t221 = (-t142 * t201 + t144 * t202) * t197 + (-t141 * t201 + t143 * t202) * t196 + (-t171 * t201 + t172 * t202) * t203;
t154 = -Icges(3,6) * t218 + t216 * t234;
t155 = Icges(3,6) * t216 + t218 * t234;
t156 = -Icges(3,5) * t218 + t216 * t236;
t157 = Icges(3,5) * t216 + t218 * t236;
t181 = Icges(3,2) * t213 + t263;
t182 = Icges(3,1) * t212 + t262;
t220 = (-t155 * t212 + t157 * t213) * V_base(4) + (-t154 * t212 + t156 * t213) * V_base(5) + (-t181 * t212 + t182 * t213) * t203;
t208 = Icges(2,4) * t218;
t195 = rSges(2,1) * t218 - rSges(2,2) * t216;
t193 = rSges(2,1) * t216 + rSges(2,2) * t218;
t189 = Icges(2,1) * t218 - t264;
t188 = Icges(2,1) * t216 + t208;
t187 = -Icges(2,2) * t216 + t208;
t186 = Icges(2,2) * t218 + t264;
t185 = Icges(2,5) * t218 - Icges(2,6) * t216;
t184 = Icges(2,5) * t216 + Icges(2,6) * t218;
t183 = rSges(3,1) * t212 + rSges(3,2) * t213;
t179 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t178 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t177 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t176 = -qJD(4) * t202 + t203;
t173 = rSges(4,1) * t201 + rSges(4,2) * t202;
t168 = t202 * t250 + t253;
t167 = -t202 * t252 + t251;
t166 = t202 * t251 - t252;
t165 = -t202 * t253 - t250;
t162 = (-qJD(4) - qJD(5)) * t202 + t203;
t159 = rSges(3,3) * t216 + t218 * t238;
t158 = -rSges(3,3) * t218 + t216 * t238;
t151 = t202 * t254 + t257;
t150 = -t202 * t256 + t255;
t149 = t202 * t255 - t256;
t148 = -t202 * t257 - t254;
t146 = rSges(4,3) * t216 + t218 * t237;
t145 = -rSges(4,3) * t218 + t216 * t237;
t135 = -rSges(5,3) * t202 + (rSges(5,1) * t217 - rSges(5,2) * t215) * t201;
t134 = V_base(5) * rSges(2,3) - t193 * t203 + t244;
t133 = t195 * t203 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t132 = -Icges(5,5) * t202 + (Icges(5,1) * t217 - Icges(5,4) * t215) * t201;
t131 = -Icges(5,6) * t202 + (Icges(5,4) * t217 - Icges(5,2) * t215) * t201;
t130 = -Icges(5,3) * t202 + (Icges(5,5) * t217 - Icges(5,6) * t215) * t201;
t129 = t218 * t246 + t164;
t128 = t216 * t246 + t163;
t126 = t193 * V_base(4) - t195 * V_base(5) + V_base(3);
t124 = -rSges(6,3) * t202 + (rSges(6,1) * t207 - rSges(6,2) * t206) * t201;
t123 = -Icges(6,5) * t202 + (Icges(6,1) * t207 - Icges(6,4) * t206) * t201;
t122 = -Icges(6,6) * t202 + (Icges(6,4) * t207 - Icges(6,2) * t206) * t201;
t121 = -Icges(6,3) * t202 + (Icges(6,5) * t207 - Icges(6,6) * t206) * t201;
t119 = -pkin(8) * t202 + t201 * t266;
t118 = rSges(5,1) * t168 + rSges(5,2) * t167 + rSges(5,3) * t258;
t117 = rSges(5,1) * t166 + rSges(5,2) * t165 + rSges(5,3) * t259;
t116 = Icges(5,1) * t168 + Icges(5,4) * t167 + Icges(5,5) * t258;
t115 = Icges(5,1) * t166 + Icges(5,4) * t165 + Icges(5,5) * t259;
t114 = Icges(5,4) * t168 + Icges(5,2) * t167 + Icges(5,6) * t258;
t113 = Icges(5,4) * t166 + Icges(5,2) * t165 + Icges(5,6) * t259;
t112 = Icges(5,5) * t168 + Icges(5,6) * t167 + Icges(5,3) * t258;
t111 = Icges(5,5) * t166 + Icges(5,6) * t165 + Icges(5,3) * t259;
t110 = pkin(4) * t253 + t218 * t229;
t109 = -pkin(4) * t252 + t216 * t229;
t108 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t258;
t107 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t259;
t106 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t258;
t105 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t259;
t104 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t258;
t103 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t259;
t102 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t258;
t101 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t259;
t100 = t183 * V_base(5) + (-t158 - t192) * t203 + t241;
t99 = t159 * t203 + (-pkin(5) - t183) * V_base(4) + t230;
t98 = t158 * V_base(4) + (-t159 - t194) * V_base(5) + t245;
t97 = t173 * t196 + (-t145 + t248) * t203 + t240;
t96 = t146 * t203 - t173 * t197 + t223;
t95 = t145 * t197 - t146 * t196 + t227;
t94 = -t117 * t176 + t135 * t163 + t225;
t93 = t118 * t176 - t135 * t164 + t222;
t92 = t117 * t164 - t118 * t163 + t224;
t91 = -t107 * t162 - t109 * t176 + t119 * t163 + t124 * t128 + t225;
t90 = t108 * t162 + t110 * t176 - t119 * t164 - t124 * t129 + t222;
t89 = t107 * t129 - t108 * t128 + t109 * t164 - t110 * t163 + t224;
t1 = m(1) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + m(2) * (t126 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + t197 * (t216 * t228 + t221 * t218) / 0.2e1 + t196 * (t216 * t221 - t218 * t228) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t164 * ((t112 * t258 + t167 * t114 + t168 * t116) * t164 + (t111 * t258 + t113 * t167 + t115 * t168) * t163 + (t130 * t258 + t131 * t167 + t132 * t168) * t176) / 0.2e1 + t163 * ((t112 * t259 + t114 * t165 + t116 * t166) * t164 + (t111 * t259 + t165 * t113 + t166 * t115) * t163 + (t130 * t259 + t131 * t165 + t132 * t166) * t176) / 0.2e1 + t176 * ((-t111 * t163 - t112 * t164 - t130 * t176) * t202 + ((-t114 * t215 + t116 * t217) * t164 + (-t113 * t215 + t115 * t217) * t163 + (-t131 * t215 + t132 * t217) * t176) * t201) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t129 * ((t102 * t258 + t150 * t104 + t151 * t106) * t129 + (t101 * t258 + t103 * t150 + t105 * t151) * t128 + (t121 * t258 + t122 * t150 + t123 * t151) * t162) / 0.2e1 + t128 * ((t102 * t259 + t104 * t148 + t106 * t149) * t129 + (t101 * t259 + t148 * t103 + t149 * t105) * t128 + (t121 * t259 + t122 * t148 + t123 * t149) * t162) / 0.2e1 + t162 * ((-t101 * t128 - t102 * t129 - t121 * t162) * t202 + ((-t104 * t206 + t106 * t207) * t129 + (-t103 * t206 + t105 * t207) * t128 + (-t122 * t206 + t123 * t207) * t162) * t201) / 0.2e1 + ((t142 * t202 + t144 * t201) * t197 + (t141 * t202 + t143 * t201) * t196 + (t154 * t213 + t156 * t212 + t184) * V_base(5) + (t155 * t213 + t157 * t212 + t185) * V_base(4) + (t202 * t171 + t201 * t172 + t213 * t181 + t212 * t182 + Icges(2,3)) * t203) * t203 / 0.2e1 + (t185 * t203 + t216 * t226 + t218 * t220 + (-t186 * t216 + t188 * t218 + Icges(1,4)) * V_base(5) + (-t216 * t187 + t218 * t189 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t184 * t203 + t216 * t220 - t218 * t226 + (t218 * t186 + t216 * t188 + Icges(1,2)) * V_base(5) + (t187 * t218 + t189 * t216 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
