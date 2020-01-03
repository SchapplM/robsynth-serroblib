% Calculate kinetic energy for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:45
% EndTime: 2019-12-31 21:52:47
% DurationCPUTime: 2.21s
% Computational Cost: add. (1359->259), mult. (1657->379), div. (0->0), fcn. (1551->8), ass. (0->133)
t269 = Icges(5,1) + Icges(6,1);
t268 = Icges(5,4) + Icges(6,4);
t267 = -Icges(6,5) - Icges(5,5);
t266 = Icges(5,2) + Icges(6,2);
t265 = -Icges(6,6) - Icges(5,6);
t264 = -Icges(6,3) - Icges(5,3);
t195 = qJ(2) + qJ(3);
t192 = cos(t195);
t200 = cos(qJ(4));
t202 = cos(qJ(1));
t233 = t200 * t202;
t197 = sin(qJ(4));
t199 = sin(qJ(1));
t236 = t197 * t199;
t152 = -t192 * t236 - t233;
t234 = t199 * t200;
t235 = t197 * t202;
t153 = t192 * t234 - t235;
t191 = sin(t195);
t238 = t191 * t199;
t263 = -t265 * t152 - t267 * t153 - t264 * t238;
t154 = -t192 * t235 + t234;
t155 = t192 * t233 + t236;
t237 = t191 * t202;
t262 = -t265 * t154 - t267 * t155 - t264 * t237;
t261 = t266 * t152 + t268 * t153 - t265 * t238;
t260 = t266 * t154 + t268 * t155 - t265 * t237;
t259 = t268 * t152 + t269 * t153 - t267 * t238;
t258 = t268 * t154 + t269 * t155 - t267 * t237;
t257 = t264 * t192 + (t265 * t197 - t267 * t200) * t191;
t256 = t265 * t192 + (-t266 * t197 + t268 * t200) * t191;
t255 = t267 * t192 + (-t268 * t197 + t269 * t200) * t191;
t198 = sin(qJ(2));
t250 = pkin(2) * t198;
t201 = cos(qJ(2));
t249 = pkin(2) * t201;
t248 = pkin(4) * t200;
t213 = qJ(5) * t191 + t192 * t248;
t245 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t238 - pkin(4) * t235 + t199 * t213;
t244 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t237 + pkin(4) * t236 + t202 * t213;
t243 = Icges(2,4) * t199;
t242 = Icges(3,4) * t198;
t241 = Icges(3,4) * t201;
t240 = Icges(4,4) * t191;
t239 = Icges(4,4) * t192;
t232 = (-qJ(5) - rSges(6,3)) * t192 + (rSges(6,1) * t200 - rSges(6,2) * t197 + t248) * t191;
t129 = -pkin(7) * t202 + t199 * t249;
t182 = t199 * pkin(1) - t202 * pkin(6);
t231 = -t129 - t182;
t230 = qJD(4) * t191;
t229 = qJD(5) * t191;
t228 = V_base(5) * pkin(5) + V_base(1);
t185 = qJD(2) * t199 + V_base(4);
t188 = V_base(6) + qJD(1);
t184 = -qJD(2) * t202 + V_base(5);
t225 = t184 * t250 + t228;
t163 = qJD(3) * t199 + t185;
t224 = pkin(3) * t192 + pkin(8) * t191;
t223 = rSges(3,1) * t201 - rSges(3,2) * t198;
t222 = rSges(4,1) * t192 - rSges(4,2) * t191;
t221 = Icges(3,1) * t201 - t242;
t220 = Icges(4,1) * t192 - t240;
t219 = -Icges(3,2) * t198 + t241;
t218 = -Icges(4,2) * t191 + t239;
t217 = Icges(3,5) * t201 - Icges(3,6) * t198;
t216 = Icges(4,5) * t192 - Icges(4,6) * t191;
t183 = t202 * pkin(1) + t199 * pkin(6);
t215 = -V_base(4) * pkin(5) + t188 * t183 + V_base(2);
t214 = V_base(4) * t182 - t183 * V_base(5) + V_base(3);
t162 = V_base(5) + (-qJD(2) - qJD(3)) * t202;
t212 = (-Icges(4,3) * t202 + t199 * t216) * t162 + (Icges(4,3) * t199 + t202 * t216) * t163 + (Icges(4,5) * t191 + Icges(4,6) * t192) * t188;
t211 = (-Icges(3,3) * t202 + t199 * t217) * t184 + (Icges(3,3) * t199 + t202 * t217) * t185 + (Icges(3,5) * t198 + Icges(3,6) * t201) * t188;
t150 = t224 * t199;
t161 = pkin(3) * t191 - pkin(8) * t192;
t210 = t162 * t161 + (-t150 + t231) * t188 + t225;
t130 = pkin(7) * t199 + t202 * t249;
t209 = t185 * t129 - t130 * t184 + t214;
t208 = t188 * t130 - t185 * t250 + t215;
t151 = t224 * t202;
t207 = t163 * t150 - t151 * t162 + t209;
t206 = t188 * t151 - t161 * t163 + t208;
t134 = -Icges(4,6) * t202 + t199 * t218;
t135 = Icges(4,6) * t199 + t202 * t218;
t136 = -Icges(4,5) * t202 + t199 * t220;
t137 = Icges(4,5) * t199 + t202 * t220;
t158 = Icges(4,2) * t192 + t240;
t159 = Icges(4,1) * t191 + t239;
t205 = (-t135 * t191 + t137 * t192) * t163 + (-t134 * t191 + t136 * t192) * t162 + (-t158 * t191 + t159 * t192) * t188;
t144 = -Icges(3,6) * t202 + t199 * t219;
t145 = Icges(3,6) * t199 + t202 * t219;
t146 = -Icges(3,5) * t202 + t199 * t221;
t147 = Icges(3,5) * t199 + t202 * t221;
t173 = Icges(3,2) * t201 + t242;
t176 = Icges(3,1) * t198 + t241;
t204 = (-t145 * t198 + t147 * t201) * t185 + (-t144 * t198 + t146 * t201) * t184 + (-t173 * t198 + t176 * t201) * t188;
t193 = Icges(2,4) * t202;
t181 = rSges(2,1) * t202 - rSges(2,2) * t199;
t180 = rSges(2,1) * t199 + rSges(2,2) * t202;
t179 = rSges(3,1) * t198 + rSges(3,2) * t201;
t178 = Icges(2,1) * t202 - t243;
t177 = Icges(2,1) * t199 + t193;
t175 = -Icges(2,2) * t199 + t193;
t174 = Icges(2,2) * t202 + t243;
t169 = -qJD(4) * t192 + t188;
t168 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t167 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t166 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t160 = rSges(4,1) * t191 + rSges(4,2) * t192;
t149 = rSges(3,3) * t199 + t202 * t223;
t148 = -rSges(3,3) * t202 + t199 * t223;
t141 = t202 * t230 + t163;
t140 = t199 * t230 + t162;
t139 = rSges(4,3) * t199 + t202 * t222;
t138 = -rSges(4,3) * t202 + t199 * t222;
t128 = -rSges(5,3) * t192 + (rSges(5,1) * t200 - rSges(5,2) * t197) * t191;
t120 = V_base(5) * rSges(2,3) - t180 * t188 + t228;
t119 = t181 * t188 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t117 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t112 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t237;
t110 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t238;
t94 = t179 * t184 + (-t148 - t182) * t188 + t228;
t93 = t149 * t188 - t179 * t185 + t215;
t92 = t148 * t185 - t149 * t184 + t214;
t91 = t160 * t162 + (-t138 + t231) * t188 + t225;
t90 = t139 * t188 - t160 * t163 + t208;
t89 = t138 * t163 - t139 * t162 + t209;
t88 = -t110 * t169 + t128 * t140 + t210;
t87 = t112 * t169 - t128 * t141 + t206;
t86 = t110 * t141 - t112 * t140 + t207;
t85 = t140 * t232 - t169 * t245 + t202 * t229 + t210;
t84 = -t141 * t232 + t169 * t244 + t199 * t229 + t206;
t83 = -qJD(5) * t192 - t140 * t244 + t141 * t245 + t207;
t1 = m(1) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t185 * (t211 * t199 + t204 * t202) / 0.2e1 + t184 * (t204 * t199 - t211 * t202) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t163 * (t212 * t199 + t205 * t202) / 0.2e1 + t162 * (t205 * t199 - t212 * t202) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t152 * t256 + t153 * t255 + t238 * t257) * t169 + (t152 * t260 + t153 * t258 + t238 * t262) * t141 + (t261 * t152 + t259 * t153 + t263 * t238) * t140) * t140 / 0.2e1 + ((t154 * t256 + t155 * t255 + t237 * t257) * t169 + (t260 * t154 + t258 * t155 + t262 * t237) * t141 + (t261 * t154 + t259 * t155 + t237 * t263) * t140) * t141 / 0.2e1 + ((-t140 * t263 - t262 * t141 - t257 * t169) * t192 + ((-t197 * t256 + t200 * t255) * t169 + (-t197 * t260 + t200 * t258) * t141 + (-t197 * t261 + t200 * t259) * t140) * t191) * t169 / 0.2e1 + ((-t174 * t199 + t177 * t202 + Icges(1,4)) * V_base(5) + (-t175 * t199 + t178 * t202 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t174 * t202 + t177 * t199 + Icges(1,2)) * V_base(5) + (t175 * t202 + t178 * t199 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t145 * t201 + t147 * t198) * t185 + (t144 * t201 + t146 * t198) * t184 + (t135 * t192 + t137 * t191) * t163 + (t134 * t192 + t136 * t191) * t162 + (t158 * t192 + t191 * t159 + t173 * t201 + t176 * t198 + Icges(2,3)) * t188) * t188 / 0.2e1 + V_base(4) * t188 * (Icges(2,5) * t202 - Icges(2,6) * t199) + V_base(5) * t188 * (Icges(2,5) * t199 + Icges(2,6) * t202) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
