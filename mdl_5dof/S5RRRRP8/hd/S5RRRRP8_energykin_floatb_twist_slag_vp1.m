% Calculate kinetic energy for
% S5RRRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:35
% EndTime: 2019-12-31 21:59:37
% DurationCPUTime: 2.85s
% Computational Cost: add. (1371->283), mult. (2009->422), div. (0->0), fcn. (1959->8), ass. (0->135)
t268 = Icges(5,1) + Icges(6,1);
t267 = Icges(5,4) + Icges(6,4);
t266 = -Icges(6,5) - Icges(5,5);
t265 = Icges(5,2) + Icges(6,2);
t264 = -Icges(6,6) - Icges(5,6);
t263 = -Icges(6,3) - Icges(5,3);
t201 = qJ(3) + qJ(4);
t195 = sin(t201);
t196 = cos(t201);
t207 = cos(qJ(1));
t204 = sin(qJ(1));
t206 = cos(qJ(2));
t237 = t204 * t206;
t151 = -t195 * t237 - t196 * t207;
t152 = -t195 * t207 + t196 * t237;
t203 = sin(qJ(2));
t239 = t203 * t204;
t262 = -t264 * t151 - t266 * t152 - t263 * t239;
t236 = t206 * t207;
t153 = -t195 * t236 + t196 * t204;
t154 = t195 * t204 + t196 * t236;
t238 = t203 * t207;
t261 = -t264 * t153 - t266 * t154 - t263 * t238;
t260 = t265 * t151 + t267 * t152 - t264 * t239;
t259 = t265 * t153 + t267 * t154 - t264 * t238;
t258 = t267 * t151 + t268 * t152 - t266 * t239;
t257 = t267 * t153 + t268 * t154 - t266 * t238;
t256 = t263 * t206 + (t264 * t195 - t266 * t196) * t203;
t255 = t264 * t206 + (-t265 * t195 + t267 * t196) * t203;
t254 = t266 * t206 + (-t267 * t195 + t268 * t196) * t203;
t205 = cos(qJ(3));
t248 = t205 * pkin(3);
t234 = pkin(4) * t196;
t215 = qJ(5) * t203 + t206 * t234;
t226 = pkin(4) * t195;
t246 = rSges(6,1) * t152 + rSges(6,2) * t151 + rSges(6,3) * t239 + t204 * t215 - t207 * t226;
t245 = rSges(6,1) * t154 + rSges(6,2) * t153 + rSges(6,3) * t238 + t204 * t226 + t207 * t215;
t244 = Icges(2,4) * t204;
t243 = Icges(3,4) * t203;
t242 = Icges(3,4) * t206;
t202 = sin(qJ(3));
t241 = t202 * t204;
t240 = t202 * t207;
t235 = (-qJ(5) - rSges(6,3)) * t206 + (rSges(6,1) * t196 - rSges(6,2) * t195 + t234) * t203;
t232 = qJD(3) * t203;
t231 = qJD(4) * t203;
t230 = qJD(5) * t203;
t229 = V_base(5) * pkin(5) + V_base(1);
t189 = qJD(2) * t204 + V_base(4);
t193 = V_base(6) + qJD(1);
t157 = t207 * t232 + t189;
t225 = pkin(2) * t206 + pkin(7) * t203;
t188 = -qJD(2) * t207 + V_base(5);
t224 = rSges(3,1) * t206 - rSges(3,2) * t203;
t223 = Icges(3,1) * t206 - t243;
t222 = -Icges(3,2) * t203 + t242;
t221 = Icges(3,5) * t206 - Icges(3,6) * t203;
t156 = t204 * t232 + t188;
t187 = pkin(1) * t207 + pkin(6) * t204;
t220 = -V_base(4) * pkin(5) + t193 * t187 + V_base(2);
t186 = pkin(1) * t204 - pkin(6) * t207;
t219 = V_base(4) * t186 - t187 * V_base(5) + V_base(3);
t218 = pkin(8) * t203 + t206 * t248;
t163 = t225 * t204;
t185 = t203 * pkin(2) - t206 * pkin(7);
t217 = t188 * t185 + (-t163 - t186) * t193 + t229;
t216 = (-Icges(3,3) * t207 + t204 * t221) * t188 + (Icges(3,3) * t204 + t207 * t221) * t189 + (Icges(3,5) * t203 + Icges(3,6) * t206) * t193;
t164 = t225 * t207;
t214 = t193 * t164 - t185 * t189 + t220;
t213 = t189 * t163 - t164 * t188 + t219;
t119 = -pkin(3) * t240 + t204 * t218;
t126 = -pkin(8) * t206 + t203 * t248;
t181 = -qJD(3) * t206 + t193;
t212 = -t119 * t181 + t156 * t126 + t217;
t120 = pkin(3) * t241 + t207 * t218;
t211 = t181 * t120 - t126 * t157 + t214;
t210 = t157 * t119 - t120 * t156 + t213;
t143 = -Icges(3,6) * t207 + t204 * t222;
t144 = Icges(3,6) * t204 + t207 * t222;
t146 = -Icges(3,5) * t207 + t204 * t223;
t147 = Icges(3,5) * t204 + t207 * t223;
t175 = Icges(3,2) * t206 + t243;
t178 = Icges(3,1) * t203 + t242;
t209 = (-t144 * t203 + t147 * t206) * t189 + (-t143 * t203 + t146 * t206) * t188 + (-t175 * t203 + t178 * t206) * t193;
t197 = Icges(2,4) * t207;
t184 = rSges(2,1) * t207 - rSges(2,2) * t204;
t183 = rSges(2,1) * t204 + rSges(2,2) * t207;
t182 = rSges(3,1) * t203 + rSges(3,2) * t206;
t180 = Icges(2,1) * t207 - t244;
t179 = Icges(2,1) * t204 + t197;
t177 = -Icges(2,2) * t204 + t197;
t176 = Icges(2,2) * t207 + t244;
t170 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t169 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t168 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t165 = (-qJD(3) - qJD(4)) * t206 + t193;
t161 = t205 * t236 + t241;
t160 = -t202 * t236 + t204 * t205;
t159 = t205 * t237 - t240;
t158 = -t202 * t237 - t205 * t207;
t150 = rSges(3,3) * t204 + t207 * t224;
t149 = -rSges(3,3) * t207 + t204 * t224;
t148 = -rSges(4,3) * t206 + (rSges(4,1) * t205 - rSges(4,2) * t202) * t203;
t145 = -Icges(4,5) * t206 + (Icges(4,1) * t205 - Icges(4,4) * t202) * t203;
t142 = -Icges(4,6) * t206 + (Icges(4,4) * t205 - Icges(4,2) * t202) * t203;
t139 = -Icges(4,3) * t206 + (Icges(4,5) * t205 - Icges(4,6) * t202) * t203;
t137 = t207 * t231 + t157;
t136 = t204 * t231 + t156;
t135 = -rSges(5,3) * t206 + (rSges(5,1) * t196 - rSges(5,2) * t195) * t203;
t125 = V_base(5) * rSges(2,3) - t183 * t193 + t229;
t124 = t184 * t193 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t123 = t183 * V_base(4) - t184 * V_base(5) + V_base(3);
t118 = rSges(4,1) * t161 + rSges(4,2) * t160 + rSges(4,3) * t238;
t117 = rSges(4,1) * t159 + rSges(4,2) * t158 + rSges(4,3) * t239;
t116 = Icges(4,1) * t161 + Icges(4,4) * t160 + Icges(4,5) * t238;
t115 = Icges(4,1) * t159 + Icges(4,4) * t158 + Icges(4,5) * t239;
t114 = Icges(4,4) * t161 + Icges(4,2) * t160 + Icges(4,6) * t238;
t113 = Icges(4,4) * t159 + Icges(4,2) * t158 + Icges(4,6) * t239;
t112 = Icges(4,5) * t161 + Icges(4,6) * t160 + Icges(4,3) * t238;
t111 = Icges(4,5) * t159 + Icges(4,6) * t158 + Icges(4,3) * t239;
t110 = rSges(5,1) * t154 + rSges(5,2) * t153 + rSges(5,3) * t238;
t108 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t239;
t92 = t182 * t188 + (-t149 - t186) * t193 + t229;
t91 = t150 * t193 - t182 * t189 + t220;
t88 = t149 * t189 - t150 * t188 + t219;
t87 = -t117 * t181 + t148 * t156 + t217;
t86 = t118 * t181 - t148 * t157 + t214;
t85 = t117 * t157 - t118 * t156 + t213;
t84 = -t108 * t165 + t135 * t136 + t212;
t83 = t110 * t165 - t135 * t137 + t211;
t82 = t108 * t137 - t110 * t136 + t210;
t81 = t136 * t235 - t165 * t246 + t207 * t230 + t212;
t80 = -t137 * t235 + t165 * t245 + t204 * t230 + t211;
t79 = -qJD(5) * t206 - t136 * t245 + t137 * t246 + t210;
t1 = m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(2) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(3) * (t88 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + t189 * (t204 * t216 + t207 * t209) / 0.2e1 + t188 * (t204 * t209 - t216 * t207) / 0.2e1 + m(4) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t157 * ((t112 * t238 + t114 * t160 + t116 * t161) * t157 + (t111 * t238 + t113 * t160 + t115 * t161) * t156 + (t139 * t238 + t142 * t160 + t145 * t161) * t181) / 0.2e1 + t156 * ((t112 * t239 + t114 * t158 + t116 * t159) * t157 + (t111 * t239 + t113 * t158 + t115 * t159) * t156 + (t139 * t239 + t142 * t158 + t145 * t159) * t181) / 0.2e1 + t181 * ((-t111 * t156 - t112 * t157 - t139 * t181) * t206 + ((-t114 * t202 + t116 * t205) * t157 + (-t113 * t202 + t115 * t205) * t156 + (-t142 * t202 + t145 * t205) * t181) * t203) / 0.2e1 + m(5) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + ((t151 * t255 + t152 * t254 + t239 * t256) * t165 + (t151 * t259 + t152 * t257 + t239 * t261) * t137 + (t260 * t151 + t258 * t152 + t262 * t239) * t136) * t136 / 0.2e1 + ((t153 * t255 + t154 * t254 + t238 * t256) * t165 + (t259 * t153 + t257 * t154 + t261 * t238) * t137 + (t260 * t153 + t258 * t154 + t238 * t262) * t136) * t137 / 0.2e1 + ((-t136 * t262 - t261 * t137 - t256 * t165) * t206 + ((-t195 * t255 + t196 * t254) * t165 + (-t195 * t259 + t196 * t257) * t137 + (-t195 * t260 + t196 * t258) * t136) * t203) * t165 / 0.2e1 + ((t144 * t206 + t147 * t203) * t189 + (t143 * t206 + t146 * t203) * t188 + (t175 * t206 + t178 * t203 + Icges(2,3)) * t193) * t193 / 0.2e1 + ((-t176 * t204 + t179 * t207 + Icges(1,4)) * V_base(5) + (-t177 * t204 + t180 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t176 * t207 + t179 * t204 + Icges(1,2)) * V_base(5) + (t177 * t207 + t180 * t204 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t193 * (Icges(2,5) * t207 - Icges(2,6) * t204) + V_base(5) * t193 * (Icges(2,5) * t204 + Icges(2,6) * t207) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
