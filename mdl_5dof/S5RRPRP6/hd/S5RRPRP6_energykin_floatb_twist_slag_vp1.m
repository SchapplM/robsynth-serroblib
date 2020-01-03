% Calculate kinetic energy for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:42
% EndTime: 2019-12-31 19:56:44
% DurationCPUTime: 2.20s
% Computational Cost: add. (1317->256), mult. (1615->357), div. (0->0), fcn. (1509->8), ass. (0->129)
t273 = Icges(5,1) + Icges(6,1);
t272 = Icges(5,4) + Icges(6,4);
t271 = -Icges(6,5) - Icges(5,5);
t270 = Icges(5,2) + Icges(6,2);
t269 = -Icges(6,6) - Icges(5,6);
t268 = Icges(3,3) + Icges(4,3);
t267 = -Icges(6,3) - Icges(5,3);
t193 = qJ(2) + pkin(8);
t186 = sin(t193);
t187 = cos(t193);
t197 = sin(qJ(2));
t200 = cos(qJ(2));
t266 = Icges(3,5) * t200 + Icges(4,5) * t187 - Icges(3,6) * t197 - Icges(4,6) * t186;
t199 = cos(qJ(4));
t201 = cos(qJ(1));
t231 = t199 * t201;
t196 = sin(qJ(4));
t198 = sin(qJ(1));
t233 = t198 * t196;
t152 = -t187 * t233 - t231;
t232 = t198 * t199;
t234 = t196 * t201;
t153 = t187 * t232 - t234;
t236 = t186 * t198;
t265 = -t269 * t152 - t271 * t153 - t267 * t236;
t154 = -t187 * t234 + t232;
t155 = t187 * t231 + t233;
t235 = t186 * t201;
t264 = -t269 * t154 - t271 * t155 - t267 * t235;
t263 = t270 * t152 + t272 * t153 - t269 * t236;
t262 = t270 * t154 + t272 * t155 - t269 * t235;
t261 = t272 * t152 + t273 * t153 - t271 * t236;
t260 = t272 * t154 + t273 * t155 - t271 * t235;
t259 = t267 * t187 + (t269 * t196 - t271 * t199) * t186;
t258 = t269 * t187 + (-t270 * t196 + t272 * t199) * t186;
t257 = t271 * t187 + (-t272 * t196 + t273 * t199) * t186;
t237 = Icges(4,4) * t187;
t216 = -Icges(4,2) * t186 + t237;
t133 = -Icges(4,6) * t201 + t198 * t216;
t134 = Icges(4,6) * t198 + t201 * t216;
t238 = Icges(4,4) * t186;
t218 = Icges(4,1) * t187 - t238;
t135 = -Icges(4,5) * t201 + t198 * t218;
t136 = Icges(4,5) * t198 + t201 * t218;
t239 = Icges(3,4) * t200;
t217 = -Icges(3,2) * t197 + t239;
t142 = -Icges(3,6) * t201 + t198 * t217;
t143 = Icges(3,6) * t198 + t201 * t217;
t240 = Icges(3,4) * t197;
t219 = Icges(3,1) * t200 - t240;
t144 = -Icges(3,5) * t201 + t198 * t219;
t145 = Icges(3,5) * t198 + t201 * t219;
t158 = Icges(4,2) * t187 + t238;
t159 = Icges(4,1) * t186 + t237;
t171 = Icges(3,2) * t200 + t240;
t174 = Icges(3,1) * t197 + t239;
t182 = -qJD(2) * t201 + V_base(5);
t183 = qJD(2) * t198 + V_base(4);
t188 = V_base(6) + qJD(1);
t256 = (-t158 * t186 + t159 * t187 - t171 * t197 + t174 * t200) * t188 + (-t134 * t186 + t136 * t187 - t143 * t197 + t145 * t200) * t183 + (-t133 * t186 + t135 * t187 - t142 * t197 + t144 * t200) * t182;
t255 = (Icges(3,5) * t197 + Icges(4,5) * t186 + Icges(3,6) * t200 + Icges(4,6) * t187) * t188 + (t268 * t198 + t266 * t201) * t183 + (t266 * t198 - t268 * t201) * t182;
t248 = pkin(2) * t197;
t247 = pkin(2) * t200;
t246 = pkin(4) * t199;
t210 = qJ(5) * t186 + t187 * t246;
t243 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t236 - pkin(4) * t234 + t198 * t210;
t242 = t155 * rSges(6,1) + t154 * rSges(6,2) + rSges(6,3) * t235 + pkin(4) * t233 + t201 * t210;
t241 = Icges(2,4) * t198;
t230 = (-qJ(5) - rSges(6,3)) * t187 + (rSges(6,1) * t199 - rSges(6,2) * t196 + t246) * t186;
t128 = -qJ(3) * t201 + t198 * t247;
t180 = t198 * pkin(1) - pkin(6) * t201;
t229 = -t128 - t180;
t228 = qJD(4) * t186;
t227 = qJD(5) * t186;
t226 = V_base(5) * pkin(5) + V_base(1);
t223 = qJD(3) * t198 + t182 * t248 + t226;
t222 = pkin(3) * t187 + pkin(7) * t186;
t221 = rSges(3,1) * t200 - rSges(3,2) * t197;
t220 = rSges(4,1) * t187 - rSges(4,2) * t186;
t181 = pkin(1) * t201 + t198 * pkin(6);
t213 = -V_base(4) * pkin(5) + t188 * t181 + V_base(2);
t212 = V_base(4) * t180 - t181 * V_base(5) + V_base(3);
t211 = t183 * t128 + t212;
t129 = qJ(3) * t198 + t201 * t247;
t207 = -qJD(3) * t201 + t188 * t129 + t213;
t148 = t222 * t198;
t161 = pkin(3) * t186 - pkin(7) * t187;
t206 = t182 * t161 + (-t148 + t229) * t188 + t223;
t149 = t222 * t201;
t205 = t183 * t148 + (-t129 - t149) * t182 + t211;
t204 = t188 * t149 + (-t161 - t248) * t183 + t207;
t191 = Icges(2,4) * t201;
t179 = rSges(2,1) * t201 - t198 * rSges(2,2);
t178 = t198 * rSges(2,1) + rSges(2,2) * t201;
t177 = rSges(3,1) * t197 + rSges(3,2) * t200;
t176 = Icges(2,1) * t201 - t241;
t175 = Icges(2,1) * t198 + t191;
t173 = -Icges(2,2) * t198 + t191;
t172 = Icges(2,2) * t201 + t241;
t167 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t166 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t165 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t164 = -qJD(4) * t187 + t188;
t160 = rSges(4,1) * t186 + rSges(4,2) * t187;
t151 = t201 * t228 + t183;
t150 = t198 * t228 + t182;
t147 = t198 * rSges(3,3) + t201 * t221;
t146 = -rSges(3,3) * t201 + t198 * t221;
t138 = t198 * rSges(4,3) + t201 * t220;
t137 = -rSges(4,3) * t201 + t198 * t220;
t127 = -rSges(5,3) * t187 + (rSges(5,1) * t199 - rSges(5,2) * t196) * t186;
t125 = V_base(5) * rSges(2,3) - t178 * t188 + t226;
t124 = t179 * t188 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t116 = t178 * V_base(4) - t179 * V_base(5) + V_base(3);
t112 = t155 * rSges(5,1) + t154 * rSges(5,2) + rSges(5,3) * t235;
t110 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t236;
t94 = t177 * t182 + (-t146 - t180) * t188 + t226;
t93 = t147 * t188 - t177 * t183 + t213;
t92 = t146 * t183 - t147 * t182 + t212;
t91 = t160 * t182 + (-t137 + t229) * t188 + t223;
t90 = t188 * t138 + (-t160 - t248) * t183 + t207;
t89 = t137 * t183 + (-t129 - t138) * t182 + t211;
t88 = -t110 * t164 + t127 * t150 + t206;
t87 = t164 * t112 - t151 * t127 + t204;
t86 = t110 * t151 - t112 * t150 + t205;
t85 = t150 * t230 - t164 * t243 + t201 * t227 + t206;
t84 = -t151 * t230 + t164 * t242 + t198 * t227 + t204;
t83 = -qJD(5) * t187 - t150 * t242 + t151 * t243 + t205;
t1 = m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t152 * t258 + t153 * t257 + t236 * t259) * t164 + (t152 * t262 + t153 * t260 + t236 * t264) * t151 + (t263 * t152 + t261 * t153 + t265 * t236) * t150) * t150 / 0.2e1 + ((t154 * t258 + t155 * t257 + t235 * t259) * t164 + (t262 * t154 + t260 * t155 + t264 * t235) * t151 + (t263 * t154 + t261 * t155 + t235 * t265) * t150) * t151 / 0.2e1 + ((-t150 * t265 - t264 * t151 - t259 * t164) * t187 + ((-t196 * t258 + t199 * t257) * t164 + (-t196 * t262 + t199 * t260) * t151 + (-t196 * t263 + t199 * t261) * t150) * t186) * t164 / 0.2e1 + (t256 * t198 - t255 * t201) * t182 / 0.2e1 + (t255 * t198 + t256 * t201) * t183 / 0.2e1 + ((-t198 * t172 + t175 * t201 + Icges(1,4)) * V_base(5) + (-t198 * t173 + t176 * t201 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t172 * t201 + t198 * t175 + Icges(1,2)) * V_base(5) + (t173 * t201 + t198 * t176 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t134 * t187 + t136 * t186 + t143 * t200 + t145 * t197) * t183 + (t133 * t187 + t135 * t186 + t142 * t200 + t144 * t197) * t182 + (t187 * t158 + t159 * t186 + t171 * t200 + t174 * t197 + Icges(2,3)) * t188) * t188 / 0.2e1 + V_base(4) * t188 * (Icges(2,5) * t201 - Icges(2,6) * t198) + V_base(5) * t188 * (Icges(2,5) * t198 + Icges(2,6) * t201) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
