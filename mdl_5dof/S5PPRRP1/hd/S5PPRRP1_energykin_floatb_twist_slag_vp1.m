% Calculate kinetic energy for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:35
% EndTime: 2019-12-05 15:06:37
% DurationCPUTime: 1.95s
% Computational Cost: add. (1250->262), mult. (1593->366), div. (0->0), fcn. (1487->8), ass. (0->134)
t270 = Icges(5,1) + Icges(6,1);
t269 = Icges(5,4) + Icges(6,4);
t268 = -Icges(6,5) - Icges(5,5);
t267 = Icges(5,2) + Icges(6,2);
t266 = -Icges(6,6) - Icges(5,6);
t265 = -Icges(6,3) - Icges(5,3);
t192 = pkin(8) + qJ(3);
t189 = cos(t192);
t196 = cos(pkin(7));
t200 = cos(qJ(4));
t234 = t196 * t200;
t194 = sin(pkin(7));
t199 = sin(qJ(4));
t237 = t194 * t199;
t152 = -t189 * t237 - t234;
t235 = t196 * t199;
t236 = t194 * t200;
t153 = t189 * t236 - t235;
t188 = sin(t192);
t239 = t188 * t194;
t264 = -t266 * t152 - t268 * t153 - t265 * t239;
t154 = -t189 * t235 + t236;
t155 = t189 * t234 + t237;
t238 = t188 * t196;
t263 = -t266 * t154 - t268 * t155 - t265 * t238;
t262 = Icges(2,5) * t196 - Icges(2,6) * t194 + Icges(1,5);
t261 = Icges(2,5) * t194 + Icges(2,6) * t196 + Icges(1,6);
t260 = t267 * t152 + t269 * t153 - t266 * t239;
t259 = t267 * t154 + t269 * t155 - t266 * t238;
t258 = t269 * t152 + t270 * t153 - t268 * t239;
t257 = t269 * t154 + t270 * t155 - t268 * t238;
t256 = t265 * t189 + (t266 * t199 - t268 * t200) * t188;
t255 = t266 * t189 + (-t267 * t199 + t269 * t200) * t188;
t254 = t268 * t189 + (-t269 * t199 + t270 * t200) * t188;
t193 = sin(pkin(8));
t250 = pkin(2) * t193;
t195 = cos(pkin(8));
t249 = pkin(2) * t195;
t248 = pkin(4) * t200;
t210 = qJ(5) * t188 + t189 * t248;
t246 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t239 - pkin(4) * t235 + t194 * t210;
t245 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t238 + pkin(4) * t237 + t196 * t210;
t244 = Icges(2,4) * t194;
t243 = Icges(3,4) * t193;
t242 = Icges(3,4) * t195;
t241 = Icges(4,4) * t188;
t240 = Icges(4,4) * t189;
t232 = (-qJ(5) - rSges(6,3)) * t189 + (rSges(6,1) * t200 - rSges(6,2) * t199 + t248) * t188;
t118 = -pkin(5) * t196 + t194 * t249;
t177 = pkin(1) * t194 - qJ(2) * t196;
t231 = -t118 - t177;
t230 = qJD(4) * t188;
t229 = qJD(5) * t188;
t228 = V_base(5) * qJ(1) + V_base(1);
t224 = qJD(1) + V_base(3);
t182 = qJD(3) * t194 + V_base(4);
t223 = qJD(2) * t194 + t228;
t222 = V_base(4) * t177 + t224;
t221 = V_base(5) * t250 + t223;
t220 = pkin(3) * t189 + pkin(6) * t188;
t181 = -qJD(3) * t196 + V_base(5);
t219 = rSges(3,1) * t195 - rSges(3,2) * t193;
t218 = rSges(4,1) * t189 - rSges(4,2) * t188;
t217 = Icges(3,1) * t195 - t243;
t216 = Icges(4,1) * t189 - t241;
t215 = -Icges(3,2) * t193 + t242;
t214 = -Icges(4,2) * t188 + t240;
t213 = Icges(3,5) * t195 - Icges(3,6) * t193;
t212 = Icges(4,5) * t189 - Icges(4,6) * t188;
t179 = pkin(1) * t196 + qJ(2) * t194;
t211 = -qJD(2) * t196 + V_base(6) * t179 + V_base(2);
t209 = (-Icges(4,3) * t196 + t194 * t212) * t181 + (Icges(4,3) * t194 + t196 * t212) * t182 + (Icges(4,5) * t188 + Icges(4,6) * t189) * V_base(6);
t119 = pkin(5) * t194 + t196 * t249;
t208 = V_base(4) * t118 + (-t119 - t179) * V_base(5) + t222;
t207 = (-Icges(3,3) * t196 + t194 * t213) * V_base(5) + (Icges(3,3) * t194 + t196 * t213) * V_base(4) + (Icges(3,5) * t193 + Icges(3,6) * t195) * V_base(6);
t148 = t220 * t194;
t160 = pkin(3) * t188 - pkin(6) * t189;
t206 = t181 * t160 + (-t148 + t231) * V_base(6) + t221;
t205 = V_base(6) * t119 + (-qJ(1) - t250) * V_base(4) + t211;
t149 = t220 * t196;
t204 = t182 * t148 - t149 * t181 + t208;
t203 = V_base(6) * t149 - t160 * t182 + t205;
t130 = -Icges(4,6) * t196 + t194 * t214;
t131 = Icges(4,6) * t194 + t196 * t214;
t132 = -Icges(4,5) * t196 + t194 * t216;
t133 = Icges(4,5) * t194 + t196 * t216;
t157 = Icges(4,2) * t189 + t241;
t158 = Icges(4,1) * t188 + t240;
t202 = (-t131 * t188 + t133 * t189) * t182 + (-t130 * t188 + t132 * t189) * t181 + (-t157 * t188 + t158 * t189) * V_base(6);
t142 = -Icges(3,6) * t196 + t194 * t215;
t143 = Icges(3,6) * t194 + t196 * t215;
t144 = -Icges(3,5) * t196 + t194 * t217;
t145 = Icges(3,5) * t194 + t196 * t217;
t169 = Icges(3,2) * t195 + t243;
t172 = Icges(3,1) * t193 + t242;
t201 = (-t143 * t193 + t145 * t195) * V_base(4) + (-t142 * t193 + t144 * t195) * V_base(5) + (-t169 * t193 + t172 * t195) * V_base(6);
t190 = Icges(2,4) * t196;
t180 = rSges(2,1) * t196 - rSges(2,2) * t194;
t178 = rSges(2,1) * t194 + rSges(2,2) * t196;
t176 = rSges(3,1) * t193 + rSges(3,2) * t195;
t175 = -qJD(4) * t189 + V_base(6);
t174 = Icges(2,1) * t196 - t244;
t173 = Icges(2,1) * t194 + t190;
t171 = -Icges(2,2) * t194 + t190;
t170 = Icges(2,2) * t196 + t244;
t165 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t164 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t163 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t159 = rSges(4,1) * t188 + rSges(4,2) * t189;
t151 = t196 * t230 + t182;
t150 = t194 * t230 + t181;
t147 = rSges(3,3) * t194 + t196 * t219;
t146 = -rSges(3,3) * t196 + t194 * t219;
t137 = rSges(4,3) * t194 + t196 * t218;
t136 = -rSges(4,3) * t196 + t194 * t218;
t135 = V_base(5) * rSges(2,3) - t178 * V_base(6) + t228;
t134 = t180 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t127 = -t189 * rSges(5,3) + (rSges(5,1) * t200 - rSges(5,2) * t199) * t188;
t114 = t178 * V_base(4) - t180 * V_base(5) + t224;
t112 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t238;
t110 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t239;
t94 = t176 * V_base(5) + (-t146 - t177) * V_base(6) + t223;
t93 = t147 * V_base(6) + (-qJ(1) - t176) * V_base(4) + t211;
t92 = t146 * V_base(4) + (-t147 - t179) * V_base(5) + t222;
t91 = t159 * t181 + (-t136 + t231) * V_base(6) + t221;
t90 = t137 * V_base(6) - t159 * t182 + t205;
t89 = t136 * t182 - t137 * t181 + t208;
t88 = -t110 * t175 + t127 * t150 + t206;
t87 = t112 * t175 - t127 * t151 + t203;
t86 = t110 * t151 - t112 * t150 + t204;
t85 = t150 * t232 - t175 * t246 + t196 * t229 + t206;
t84 = -t151 * t232 + t175 * t245 + t194 * t229 + t203;
t83 = -qJD(5) * t189 - t150 * t245 + t151 * t246 + t204;
t1 = m(1) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(2) * (t114 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t182 * (t209 * t194 + t202 * t196) / 0.2e1 + t181 * (t202 * t194 - t209 * t196) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t255 * t152 + t254 * t153 + t256 * t239) * t175 + (t259 * t152 + t257 * t153 + t263 * t239) * t151 + (t260 * t152 + t258 * t153 + t264 * t239) * t150) * t150 / 0.2e1 + ((t255 * t154 + t254 * t155 + t256 * t238) * t175 + (t259 * t154 + t257 * t155 + t263 * t238) * t151 + (t260 * t154 + t258 * t155 + t264 * t238) * t150) * t151 / 0.2e1 + ((-t264 * t150 - t263 * t151 - t256 * t175) * t189 + ((-t255 * t199 + t254 * t200) * t175 + (-t259 * t199 + t257 * t200) * t151 + (-t260 * t199 + t258 * t200) * t150) * t188) * t175 / 0.2e1 + (t194 * t207 + t196 * t201 + t262 * V_base(6) + (-t170 * t194 + t173 * t196 + Icges(1,4)) * V_base(5) + (-t171 * t194 + t174 * t196 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t194 * t201 - t196 * t207 + t261 * V_base(6) + (t170 * t196 + t173 * t194 + Icges(1,2)) * V_base(5) + (t171 * t196 + t174 * t194 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t131 * t189 + t133 * t188) * t182 + (t130 * t189 + t132 * t188) * t181 + (t142 * t195 + t144 * t193 + t261) * V_base(5) + (t143 * t195 + t145 * t193 + t262) * V_base(4) + (t157 * t189 + t158 * t188 + t169 * t195 + t172 * t193 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1;
T = t1;
