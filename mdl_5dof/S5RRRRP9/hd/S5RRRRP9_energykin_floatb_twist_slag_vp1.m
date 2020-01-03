% Calculate kinetic energy for
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:43
% EndTime: 2019-12-31 22:03:45
% DurationCPUTime: 2.78s
% Computational Cost: add. (1336->278), mult. (1976->414), div. (0->0), fcn. (1936->8), ass. (0->133)
t271 = Icges(5,1) + Icges(6,1);
t270 = -Icges(5,4) + Icges(6,5);
t269 = Icges(6,4) + Icges(5,5);
t268 = Icges(5,2) + Icges(6,3);
t267 = -Icges(6,6) + Icges(5,6);
t266 = -Icges(5,3) - Icges(6,2);
t265 = rSges(6,1) + pkin(4);
t264 = rSges(6,3) + qJ(5);
t208 = qJ(3) + qJ(4);
t204 = sin(t208);
t205 = cos(t208);
t214 = cos(qJ(1));
t211 = sin(qJ(1));
t213 = cos(qJ(2));
t241 = t211 * t213;
t160 = t204 * t241 + t205 * t214;
t161 = -t204 * t214 + t205 * t241;
t210 = sin(qJ(2));
t243 = t210 * t211;
t263 = t268 * t160 + t270 * t161 - t267 * t243;
t240 = t213 * t214;
t162 = t204 * t240 - t211 * t205;
t163 = t204 * t211 + t205 * t240;
t242 = t210 * t214;
t262 = t268 * t162 + t270 * t163 - t267 * t242;
t261 = -t267 * t160 + t269 * t161 - t266 * t243;
t260 = -t267 * t162 + t269 * t163 - t266 * t242;
t259 = t270 * t160 + t271 * t161 + t269 * t243;
t258 = t270 * t162 + t271 * t163 + t269 * t242;
t257 = t267 * t213 + (t268 * t204 + t270 * t205) * t210;
t256 = t266 * t213 + (-t267 * t204 + t269 * t205) * t210;
t255 = -t269 * t213 + (t270 * t204 + t271 * t205) * t210;
t212 = cos(qJ(3));
t250 = pkin(3) * t212;
t248 = Icges(2,4) * t211;
t247 = Icges(3,4) * t210;
t246 = Icges(3,4) * t213;
t209 = sin(qJ(3));
t245 = t209 * t211;
t244 = t209 * t214;
t239 = rSges(6,2) * t243 + t264 * t160 + t161 * t265;
t238 = rSges(6,2) * t242 + t264 * t162 + t163 * t265;
t237 = -rSges(6,2) * t213 + (t264 * t204 + t205 * t265) * t210;
t236 = qJD(3) * t210;
t235 = qJD(4) * t210;
t234 = V_base(5) * pkin(5) + V_base(1);
t197 = qJD(2) * t211 + V_base(4);
t202 = V_base(6) + qJD(1);
t167 = t214 * t236 + t197;
t231 = pkin(2) * t213 + pkin(7) * t210;
t196 = -qJD(2) * t214 + V_base(5);
t230 = rSges(3,1) * t213 - rSges(3,2) * t210;
t229 = Icges(3,1) * t213 - t247;
t228 = -Icges(3,2) * t210 + t246;
t227 = Icges(3,5) * t213 - Icges(3,6) * t210;
t166 = t211 * t236 + t196;
t195 = pkin(1) * t214 + pkin(6) * t211;
t226 = -V_base(4) * pkin(5) + t202 * t195 + V_base(2);
t194 = pkin(1) * t211 - pkin(6) * t214;
t225 = V_base(4) * t194 - t195 * V_base(5) + V_base(3);
t224 = pkin(8) * t210 + t213 * t250;
t173 = t231 * t211;
t193 = t210 * pkin(2) - t213 * pkin(7);
t223 = t196 * t193 + (-t173 - t194) * t202 + t234;
t222 = (-Icges(3,3) * t214 + t211 * t227) * t196 + (Icges(3,3) * t211 + t214 * t227) * t197 + (Icges(3,5) * t210 + Icges(3,6) * t213) * t202;
t174 = t231 * t214;
t221 = t202 * t174 - t193 * t197 + t226;
t220 = t197 * t173 - t174 * t196 + t225;
t129 = -pkin(3) * t244 + t211 * t224;
t135 = -pkin(8) * t213 + t210 * t250;
t189 = -qJD(3) * t213 + t202;
t219 = -t129 * t189 + t166 * t135 + t223;
t130 = pkin(3) * t245 + t214 * t224;
t218 = t189 * t130 - t135 * t167 + t221;
t217 = t167 * t129 - t130 * t166 + t220;
t152 = -Icges(3,6) * t214 + t211 * t228;
t153 = Icges(3,6) * t211 + t214 * t228;
t155 = -Icges(3,5) * t214 + t211 * t229;
t156 = Icges(3,5) * t211 + t214 * t229;
t183 = Icges(3,2) * t213 + t247;
t186 = Icges(3,1) * t210 + t246;
t216 = (-t153 * t210 + t156 * t213) * t197 + (-t152 * t210 + t155 * t213) * t196 + (-t183 * t210 + t186 * t213) * t202;
t206 = Icges(2,4) * t214;
t192 = rSges(2,1) * t214 - rSges(2,2) * t211;
t191 = rSges(2,1) * t211 + rSges(2,2) * t214;
t190 = rSges(3,1) * t210 + rSges(3,2) * t213;
t188 = Icges(2,1) * t214 - t248;
t187 = Icges(2,1) * t211 + t206;
t185 = -Icges(2,2) * t211 + t206;
t184 = Icges(2,2) * t214 + t248;
t179 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t178 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t177 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t175 = (-qJD(3) - qJD(4)) * t213 + t202;
t171 = t212 * t240 + t245;
t170 = -t209 * t240 + t211 * t212;
t169 = t212 * t241 - t244;
t168 = -t209 * t241 - t212 * t214;
t159 = rSges(3,3) * t211 + t214 * t230;
t158 = -rSges(3,3) * t214 + t211 * t230;
t157 = -rSges(4,3) * t213 + (rSges(4,1) * t212 - rSges(4,2) * t209) * t210;
t154 = -Icges(4,5) * t213 + (Icges(4,1) * t212 - Icges(4,4) * t209) * t210;
t151 = -Icges(4,6) * t213 + (Icges(4,4) * t212 - Icges(4,2) * t209) * t210;
t148 = -Icges(4,3) * t213 + (Icges(4,5) * t212 - Icges(4,6) * t209) * t210;
t146 = t214 * t235 + t167;
t145 = t211 * t235 + t166;
t144 = -rSges(5,3) * t213 + (rSges(5,1) * t205 - rSges(5,2) * t204) * t210;
t134 = V_base(5) * rSges(2,3) - t191 * t202 + t234;
t133 = t192 * t202 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t132 = t191 * V_base(4) - t192 * V_base(5) + V_base(3);
t128 = rSges(4,1) * t171 + rSges(4,2) * t170 + rSges(4,3) * t242;
t127 = rSges(4,1) * t169 + rSges(4,2) * t168 + rSges(4,3) * t243;
t126 = Icges(4,1) * t171 + Icges(4,4) * t170 + Icges(4,5) * t242;
t125 = Icges(4,1) * t169 + Icges(4,4) * t168 + Icges(4,5) * t243;
t124 = Icges(4,4) * t171 + Icges(4,2) * t170 + Icges(4,6) * t242;
t123 = Icges(4,4) * t169 + Icges(4,2) * t168 + Icges(4,6) * t243;
t122 = Icges(4,5) * t171 + Icges(4,6) * t170 + Icges(4,3) * t242;
t121 = Icges(4,5) * t169 + Icges(4,6) * t168 + Icges(4,3) * t243;
t118 = rSges(5,1) * t163 - rSges(5,2) * t162 + rSges(5,3) * t242;
t116 = rSges(5,1) * t161 - rSges(5,2) * t160 + rSges(5,3) * t243;
t100 = t190 * t196 + (-t158 - t194) * t202 + t234;
t99 = t159 * t202 - t190 * t197 + t226;
t98 = t158 * t197 - t159 * t196 + t225;
t97 = -t127 * t189 + t157 * t166 + t223;
t96 = t128 * t189 - t157 * t167 + t221;
t95 = t127 * t167 - t128 * t166 + t220;
t94 = -t116 * t175 + t144 * t145 + t219;
t93 = t118 * t175 - t144 * t146 + t218;
t92 = t116 * t146 - t118 * t145 + t217;
t91 = qJD(5) * t162 + t145 * t237 - t175 * t239 + t219;
t90 = qJD(5) * t160 - t146 * t237 + t175 * t238 + t218;
t89 = qJD(5) * t204 * t210 - t145 * t238 + t146 * t239 + t217;
t1 = m(1) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + m(2) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t197 * (t222 * t211 + t216 * t214) / 0.2e1 + t196 * (t216 * t211 - t222 * t214) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + t167 * ((t122 * t242 + t170 * t124 + t171 * t126) * t167 + (t121 * t242 + t123 * t170 + t125 * t171) * t166 + (t148 * t242 + t151 * t170 + t154 * t171) * t189) / 0.2e1 + t166 * ((t122 * t243 + t124 * t168 + t126 * t169) * t167 + (t121 * t243 + t168 * t123 + t169 * t125) * t166 + (t148 * t243 + t151 * t168 + t154 * t169) * t189) / 0.2e1 + t189 * ((-t121 * t166 - t122 * t167 - t148 * t189) * t213 + ((-t124 * t209 + t126 * t212) * t167 + (-t123 * t209 + t125 * t212) * t166 + (-t151 * t209 + t154 * t212) * t189) * t210) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + ((t160 * t257 + t161 * t255 + t243 * t256) * t175 + (t160 * t262 + t161 * t258 + t243 * t260) * t146 + (t263 * t160 + t259 * t161 + t261 * t243) * t145) * t145 / 0.2e1 + ((t162 * t257 + t163 * t255 + t242 * t256) * t175 + (t262 * t162 + t258 * t163 + t260 * t242) * t146 + (t162 * t263 + t259 * t163 + t261 * t242) * t145) * t146 / 0.2e1 + ((-t145 * t261 - t146 * t260 - t175 * t256) * t213 + ((t204 * t257 + t205 * t255) * t175 + (t204 * t262 + t205 * t258) * t146 + (t204 * t263 + t259 * t205) * t145) * t210) * t175 / 0.2e1 + ((t153 * t213 + t156 * t210) * t197 + (t152 * t213 + t155 * t210) * t196 + (t183 * t213 + t186 * t210 + Icges(2,3)) * t202) * t202 / 0.2e1 + ((-t184 * t211 + t187 * t214 + Icges(1,4)) * V_base(5) + (-t185 * t211 + t188 * t214 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t184 * t214 + t187 * t211 + Icges(1,2)) * V_base(5) + (t185 * t214 + t188 * t211 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t202 * (Icges(2,5) * t214 - Icges(2,6) * t211) + V_base(5) * t202 * (Icges(2,5) * t211 + Icges(2,6) * t214) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
