% Calculate kinetic energy for
% S5RPRRP10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:51
% EndTime: 2019-12-31 18:50:53
% DurationCPUTime: 2.34s
% Computational Cost: add. (1295->262), mult. (1593->371), div. (0->0), fcn. (1487->8), ass. (0->134)
t268 = Icges(5,1) + Icges(6,1);
t267 = Icges(5,4) + Icges(6,4);
t266 = -Icges(6,5) - Icges(5,5);
t265 = Icges(5,2) + Icges(6,2);
t264 = -Icges(6,6) - Icges(5,6);
t263 = -Icges(6,3) - Icges(5,3);
t193 = pkin(8) + qJ(3);
t187 = cos(t193);
t200 = cos(qJ(4));
t201 = cos(qJ(1));
t233 = t200 * t201;
t198 = sin(qJ(4));
t199 = sin(qJ(1));
t235 = t199 * t198;
t152 = -t187 * t235 - t233;
t234 = t199 * t200;
t236 = t198 * t201;
t153 = t187 * t234 - t236;
t186 = sin(t193);
t238 = t186 * t199;
t262 = -t264 * t152 - t266 * t153 - t263 * t238;
t154 = -t187 * t236 + t234;
t155 = t187 * t233 + t235;
t237 = t186 * t201;
t261 = -t264 * t154 - t266 * t155 - t263 * t237;
t260 = t265 * t152 + t267 * t153 - t264 * t238;
t259 = t265 * t154 + t267 * t155 - t264 * t237;
t258 = t267 * t152 + t268 * t153 - t266 * t238;
t257 = t267 * t154 + t268 * t155 - t266 * t237;
t256 = t263 * t187 + (t264 * t198 - t266 * t200) * t186;
t255 = t264 * t187 + (-t265 * t198 + t267 * t200) * t186;
t254 = t266 * t187 + (-t267 * t198 + t268 * t200) * t186;
t194 = sin(pkin(8));
t249 = pkin(2) * t194;
t195 = cos(pkin(8));
t248 = pkin(2) * t195;
t247 = pkin(4) * t200;
t211 = qJ(5) * t186 + t187 * t247;
t245 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t238 - pkin(4) * t236 + t199 * t211;
t244 = t155 * rSges(6,1) + t154 * rSges(6,2) + rSges(6,3) * t237 + pkin(4) * t235 + t201 * t211;
t243 = Icges(2,4) * t199;
t242 = Icges(3,4) * t194;
t241 = Icges(3,4) * t195;
t240 = Icges(4,4) * t186;
t239 = Icges(4,4) * t187;
t231 = (-qJ(5) - rSges(6,3)) * t187 + (rSges(6,1) * t200 - rSges(6,2) * t198 + t247) * t186;
t128 = -pkin(6) * t201 + t199 * t248;
t177 = t199 * pkin(1) - qJ(2) * t201;
t230 = -t128 - t177;
t229 = qJD(4) * t186;
t228 = qJD(5) * t186;
t227 = V_base(4) * t177 + V_base(3);
t226 = V_base(5) * pkin(5) + V_base(1);
t182 = qJD(3) * t199 + V_base(4);
t188 = V_base(6) + qJD(1);
t223 = qJD(2) * t199 + t226;
t222 = V_base(5) * t249 + t223;
t221 = pkin(3) * t187 + pkin(7) * t186;
t181 = -qJD(3) * t201 + V_base(5);
t220 = rSges(3,1) * t195 - rSges(3,2) * t194;
t219 = rSges(4,1) * t187 - rSges(4,2) * t186;
t218 = Icges(3,1) * t195 - t242;
t217 = Icges(4,1) * t187 - t240;
t216 = -Icges(3,2) * t194 + t241;
t215 = -Icges(4,2) * t186 + t239;
t214 = Icges(3,5) * t195 - Icges(3,6) * t194;
t213 = Icges(4,5) * t187 - Icges(4,6) * t186;
t179 = pkin(1) * t201 + t199 * qJ(2);
t212 = -qJD(2) * t201 + t188 * t179 + V_base(2);
t210 = (-Icges(4,3) * t201 + t199 * t213) * t181 + (Icges(4,3) * t199 + t201 * t213) * t182 + (Icges(4,5) * t186 + Icges(4,6) * t187) * t188;
t129 = pkin(6) * t199 + t201 * t248;
t209 = V_base(4) * t128 + (-t129 - t179) * V_base(5) + t227;
t208 = (-Icges(3,3) * t201 + t199 * t214) * V_base(5) + (Icges(3,3) * t199 + t201 * t214) * V_base(4) + (Icges(3,5) * t194 + Icges(3,6) * t195) * t188;
t148 = t221 * t199;
t161 = pkin(3) * t186 - pkin(7) * t187;
t207 = t181 * t161 + (-t148 + t230) * t188 + t222;
t149 = t221 * t201;
t206 = t182 * t148 - t149 * t181 + t209;
t205 = t188 * t129 + (-pkin(5) - t249) * V_base(4) + t212;
t204 = t188 * t149 - t182 * t161 + t205;
t133 = -Icges(4,6) * t201 + t199 * t215;
t134 = Icges(4,6) * t199 + t201 * t215;
t135 = -Icges(4,5) * t201 + t199 * t217;
t136 = Icges(4,5) * t199 + t201 * t217;
t158 = Icges(4,2) * t187 + t240;
t159 = Icges(4,1) * t186 + t239;
t203 = (-t134 * t186 + t136 * t187) * t182 + (-t133 * t186 + t135 * t187) * t181 + (-t158 * t186 + t159 * t187) * t188;
t142 = -Icges(3,6) * t201 + t199 * t216;
t143 = Icges(3,6) * t199 + t201 * t216;
t144 = -Icges(3,5) * t201 + t199 * t218;
t145 = Icges(3,5) * t199 + t201 * t218;
t168 = Icges(3,2) * t195 + t242;
t169 = Icges(3,1) * t194 + t241;
t202 = (-t143 * t194 + t145 * t195) * V_base(4) + (-t142 * t194 + t144 * t195) * V_base(5) + (-t168 * t194 + t169 * t195) * t188;
t191 = Icges(2,4) * t201;
t180 = rSges(2,1) * t201 - t199 * rSges(2,2);
t178 = t199 * rSges(2,1) + rSges(2,2) * t201;
t176 = Icges(2,1) * t201 - t243;
t175 = Icges(2,1) * t199 + t191;
t174 = -Icges(2,2) * t199 + t191;
t173 = Icges(2,2) * t201 + t243;
t172 = Icges(2,5) * t201 - Icges(2,6) * t199;
t171 = Icges(2,5) * t199 + Icges(2,6) * t201;
t170 = rSges(3,1) * t194 + rSges(3,2) * t195;
t166 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t165 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t164 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t163 = -qJD(4) * t187 + t188;
t160 = rSges(4,1) * t186 + rSges(4,2) * t187;
t151 = t201 * t229 + t182;
t150 = t199 * t229 + t181;
t147 = t199 * rSges(3,3) + t201 * t220;
t146 = -rSges(3,3) * t201 + t199 * t220;
t138 = t199 * rSges(4,3) + t201 * t219;
t137 = -rSges(4,3) * t201 + t199 * t219;
t127 = -rSges(5,3) * t187 + (rSges(5,1) * t200 - rSges(5,2) * t198) * t186;
t125 = V_base(5) * rSges(2,3) - t178 * t188 + t226;
t124 = t180 * t188 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t116 = t178 * V_base(4) - t180 * V_base(5) + V_base(3);
t112 = t155 * rSges(5,1) + t154 * rSges(5,2) + rSges(5,3) * t237;
t110 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t238;
t94 = t170 * V_base(5) + (-t146 - t177) * t188 + t223;
t93 = t188 * t147 + (-pkin(5) - t170) * V_base(4) + t212;
t92 = t146 * V_base(4) + (-t147 - t179) * V_base(5) + t227;
t91 = t160 * t181 + (-t137 + t230) * t188 + t222;
t90 = t188 * t138 - t182 * t160 + t205;
t89 = t137 * t182 - t138 * t181 + t209;
t88 = -t110 * t163 + t127 * t150 + t207;
t87 = t163 * t112 - t151 * t127 + t204;
t86 = t110 * t151 - t112 * t150 + t206;
t85 = t150 * t231 - t163 * t245 + t201 * t228 + t207;
t84 = -t151 * t231 + t163 * t244 + t199 * t228 + t204;
t83 = -qJD(5) * t187 - t150 * t244 + t151 * t245 + t206;
t1 = m(1) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t182 * (t199 * t210 + t201 * t203) / 0.2e1 + t181 * (t199 * t203 - t201 * t210) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t255 * t152 + t254 * t153 + t256 * t238) * t163 + (t259 * t152 + t257 * t153 + t261 * t238) * t151 + (t260 * t152 + t258 * t153 + t262 * t238) * t150) * t150 / 0.2e1 + ((t255 * t154 + t254 * t155 + t256 * t237) * t163 + (t259 * t154 + t257 * t155 + t261 * t237) * t151 + (t260 * t154 + t258 * t155 + t262 * t237) * t150) * t151 / 0.2e1 + ((-t262 * t150 - t261 * t151 - t256 * t163) * t187 + ((-t255 * t198 + t254 * t200) * t163 + (-t259 * t198 + t257 * t200) * t151 + (-t260 * t198 + t258 * t200) * t150) * t186) * t163 / 0.2e1 + ((t134 * t187 + t136 * t186) * t182 + (t133 * t187 + t135 * t186) * t181 + (t142 * t195 + t144 * t194 + t171) * V_base(5) + (t143 * t195 + t145 * t194 + t172) * V_base(4) + (t158 * t187 + t159 * t186 + t168 * t195 + t169 * t194 + Icges(2,3)) * t188) * t188 / 0.2e1 + (t172 * t188 + t199 * t208 + t201 * t202 + (-t199 * t173 + t175 * t201 + Icges(1,4)) * V_base(5) + (-t199 * t174 + t176 * t201 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t171 * t188 + t199 * t202 - t201 * t208 + (t173 * t201 + t199 * t175 + Icges(1,2)) * V_base(5) + (t174 * t201 + t199 * t176 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
