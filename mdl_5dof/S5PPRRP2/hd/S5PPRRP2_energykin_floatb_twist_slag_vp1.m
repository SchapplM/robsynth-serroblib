% Calculate kinetic energy for
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:28
% EndTime: 2019-12-05 15:08:30
% DurationCPUTime: 1.99s
% Computational Cost: add. (1220->257), mult. (1584->359), div. (0->0), fcn. (1488->8), ass. (0->133)
t269 = Icges(5,1) + Icges(6,1);
t268 = -Icges(5,4) + Icges(6,5);
t267 = Icges(6,4) + Icges(5,5);
t266 = Icges(5,2) + Icges(6,3);
t265 = -Icges(6,6) + Icges(5,6);
t264 = -Icges(5,3) - Icges(6,2);
t263 = rSges(6,1) + pkin(4);
t262 = rSges(6,3) + qJ(5);
t194 = pkin(8) + qJ(3);
t191 = cos(t194);
t198 = cos(pkin(7));
t201 = cos(qJ(4));
t235 = t198 * t201;
t196 = sin(pkin(7));
t200 = sin(qJ(4));
t238 = t196 * t200;
t154 = t191 * t238 + t235;
t236 = t198 * t200;
t237 = t196 * t201;
t155 = t191 * t237 - t236;
t190 = sin(t194);
t240 = t190 * t196;
t261 = -t265 * t154 + t267 * t155 - t264 * t240;
t260 = t266 * t154 + t268 * t155 - t265 * t240;
t156 = t191 * t236 - t237;
t157 = t191 * t235 + t238;
t239 = t190 * t198;
t259 = t266 * t156 + t268 * t157 - t265 * t239;
t258 = Icges(2,5) * t198 - Icges(2,6) * t196 + Icges(1,5);
t257 = Icges(2,5) * t196 + Icges(2,6) * t198 + Icges(1,6);
t256 = -t265 * t156 + t267 * t157 - t264 * t239;
t255 = t268 * t154 + t269 * t155 + t267 * t240;
t254 = t268 * t156 + t269 * t157 + t267 * t239;
t253 = t265 * t191 + (t266 * t200 + t268 * t201) * t190;
t252 = t264 * t191 + (-t265 * t200 + t267 * t201) * t190;
t251 = -t267 * t191 + (t268 * t200 + t269 * t201) * t190;
t195 = sin(pkin(8));
t247 = pkin(2) * t195;
t197 = cos(pkin(8));
t246 = pkin(2) * t197;
t245 = Icges(2,4) * t196;
t244 = Icges(3,4) * t195;
t243 = Icges(3,4) * t197;
t242 = Icges(4,4) * t190;
t241 = Icges(4,4) * t191;
t233 = rSges(6,2) * t240 + t262 * t154 + t263 * t155;
t232 = rSges(6,2) * t239 + t262 * t156 + t263 * t157;
t119 = -pkin(5) * t198 + t196 * t246;
t179 = pkin(1) * t196 - qJ(2) * t198;
t231 = -t119 - t179;
t230 = -t191 * rSges(6,2) + (t262 * t200 + t263 * t201) * t190;
t229 = qJD(4) * t190;
t228 = V_base(5) * qJ(1) + V_base(1);
t224 = qJD(1) + V_base(3);
t184 = qJD(3) * t196 + V_base(4);
t223 = qJD(2) * t196 + t228;
t222 = V_base(4) * t179 + t224;
t221 = V_base(5) * t247 + t223;
t220 = pkin(3) * t191 + pkin(6) * t190;
t183 = -qJD(3) * t198 + V_base(5);
t219 = rSges(3,1) * t197 - rSges(3,2) * t195;
t218 = rSges(4,1) * t191 - rSges(4,2) * t190;
t217 = Icges(3,1) * t197 - t244;
t216 = Icges(4,1) * t191 - t242;
t215 = -Icges(3,2) * t195 + t243;
t214 = -Icges(4,2) * t190 + t241;
t213 = Icges(3,5) * t197 - Icges(3,6) * t195;
t212 = Icges(4,5) * t191 - Icges(4,6) * t190;
t181 = pkin(1) * t198 + qJ(2) * t196;
t211 = -qJD(2) * t198 + V_base(6) * t181 + V_base(2);
t210 = (-Icges(4,3) * t198 + t196 * t212) * t183 + (Icges(4,3) * t196 + t198 * t212) * t184 + (Icges(4,5) * t190 + Icges(4,6) * t191) * V_base(6);
t120 = pkin(5) * t196 + t198 * t246;
t209 = V_base(4) * t119 + (-t120 - t181) * V_base(5) + t222;
t208 = (-Icges(3,3) * t198 + t196 * t213) * V_base(5) + (Icges(3,3) * t196 + t198 * t213) * V_base(4) + (Icges(3,5) * t195 + Icges(3,6) * t197) * V_base(6);
t149 = t220 * t196;
t162 = pkin(3) * t190 - pkin(6) * t191;
t207 = t183 * t162 + (-t149 + t231) * V_base(6) + t221;
t206 = V_base(6) * t120 + (-qJ(1) - t247) * V_base(4) + t211;
t150 = t220 * t198;
t205 = t184 * t149 - t150 * t183 + t209;
t204 = V_base(6) * t150 - t162 * t184 + t206;
t131 = -Icges(4,6) * t198 + t196 * t214;
t132 = Icges(4,6) * t196 + t198 * t214;
t133 = -Icges(4,5) * t198 + t196 * t216;
t134 = Icges(4,5) * t196 + t198 * t216;
t159 = Icges(4,2) * t191 + t242;
t160 = Icges(4,1) * t190 + t241;
t203 = (-t132 * t190 + t134 * t191) * t184 + (-t131 * t190 + t133 * t191) * t183 + (-t159 * t190 + t160 * t191) * V_base(6);
t143 = -Icges(3,6) * t198 + t196 * t215;
t144 = Icges(3,6) * t196 + t198 * t215;
t145 = -Icges(3,5) * t198 + t196 * t217;
t146 = Icges(3,5) * t196 + t198 * t217;
t171 = Icges(3,2) * t197 + t244;
t174 = Icges(3,1) * t195 + t243;
t202 = (-t144 * t195 + t146 * t197) * V_base(4) + (-t143 * t195 + t145 * t197) * V_base(5) + (-t171 * t195 + t174 * t197) * V_base(6);
t192 = Icges(2,4) * t198;
t182 = rSges(2,1) * t198 - rSges(2,2) * t196;
t180 = rSges(2,1) * t196 + rSges(2,2) * t198;
t178 = rSges(3,1) * t195 + rSges(3,2) * t197;
t177 = -qJD(4) * t191 + V_base(6);
t176 = Icges(2,1) * t198 - t245;
t175 = Icges(2,1) * t196 + t192;
t173 = -Icges(2,2) * t196 + t192;
t172 = Icges(2,2) * t198 + t245;
t167 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t166 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t165 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = rSges(4,1) * t190 + rSges(4,2) * t191;
t153 = t198 * t229 + t184;
t152 = t196 * t229 + t183;
t148 = rSges(3,3) * t196 + t198 * t219;
t147 = -rSges(3,3) * t198 + t196 * t219;
t138 = rSges(4,3) * t196 + t198 * t218;
t137 = -rSges(4,3) * t198 + t196 * t218;
t136 = V_base(5) * rSges(2,3) - t180 * V_base(6) + t228;
t135 = t182 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t128 = -t191 * rSges(5,3) + (rSges(5,1) * t201 - rSges(5,2) * t200) * t190;
t115 = t180 * V_base(4) - t182 * V_base(5) + t224;
t112 = rSges(5,1) * t157 - rSges(5,2) * t156 + rSges(5,3) * t239;
t110 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t240;
t96 = t178 * V_base(5) + (-t147 - t179) * V_base(6) + t223;
t95 = t148 * V_base(6) + (-qJ(1) - t178) * V_base(4) + t211;
t94 = t147 * V_base(4) + (-t148 - t181) * V_base(5) + t222;
t93 = t161 * t183 + (-t137 + t231) * V_base(6) + t221;
t92 = t138 * V_base(6) - t161 * t184 + t206;
t91 = t137 * t184 - t138 * t183 + t209;
t90 = -t110 * t177 + t128 * t152 + t207;
t89 = t112 * t177 - t128 * t153 + t204;
t88 = t110 * t153 - t112 * t152 + t205;
t87 = qJD(5) * t156 + t152 * t230 - t177 * t233 + t207;
t86 = qJD(5) * t154 - t153 * t230 + t177 * t232 + t204;
t85 = qJD(5) * t190 * t200 - t152 * t232 + t153 * t233 + t205;
t1 = m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(2) * (t115 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t184 * (t210 * t196 + t203 * t198) / 0.2e1 + t183 * (t203 * t196 - t210 * t198) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((t154 * t253 + t155 * t251 + t240 * t252) * t177 + (t154 * t259 + t155 * t254 + t240 * t256) * t153 + (t260 * t154 + t255 * t155 + t261 * t240) * t152) * t152 / 0.2e1 + ((t156 * t253 + t157 * t251 + t239 * t252) * t177 + (t259 * t156 + t254 * t157 + t256 * t239) * t153 + (t260 * t156 + t255 * t157 + t239 * t261) * t152) * t153 / 0.2e1 + ((-t152 * t261 - t256 * t153 - t252 * t177) * t191 + ((t200 * t253 + t201 * t251) * t177 + (t200 * t259 + t201 * t254) * t153 + (t200 * t260 + t201 * t255) * t152) * t190) * t177 / 0.2e1 + (t208 * t196 + t202 * t198 + t258 * V_base(6) + (-t172 * t196 + t175 * t198 + Icges(1,4)) * V_base(5) + (-t173 * t196 + t176 * t198 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t202 * t196 - t208 * t198 + t257 * V_base(6) + (t172 * t198 + t175 * t196 + Icges(1,2)) * V_base(5) + (t173 * t198 + t176 * t196 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t132 * t191 + t134 * t190) * t184 + (t131 * t191 + t133 * t190) * t183 + (t143 * t197 + t145 * t195 + t257) * V_base(5) + (t144 * t197 + t146 * t195 + t258) * V_base(4) + (t159 * t191 + t160 * t190 + t171 * t197 + t174 * t195 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1;
T = t1;
