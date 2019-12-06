% Calculate kinetic energy for
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:32
% EndTime: 2019-12-05 16:45:35
% DurationCPUTime: 2.24s
% Computational Cost: add. (1284->254), mult. (1648->367), div. (0->0), fcn. (1552->8), ass. (0->132)
t270 = Icges(5,1) + Icges(6,1);
t269 = -Icges(5,4) + Icges(6,5);
t268 = Icges(6,4) + Icges(5,5);
t267 = Icges(5,2) + Icges(6,3);
t266 = -Icges(6,6) + Icges(5,6);
t265 = -Icges(5,3) - Icges(6,2);
t264 = rSges(6,1) + pkin(4);
t263 = rSges(6,3) + qJ(5);
t196 = qJ(2) + qJ(3);
t194 = cos(t196);
t198 = cos(pkin(8));
t201 = cos(qJ(4));
t235 = t198 * t201;
t197 = sin(pkin(8));
t199 = sin(qJ(4));
t238 = t197 * t199;
t154 = t194 * t238 + t235;
t236 = t198 * t199;
t237 = t197 * t201;
t155 = t194 * t237 - t236;
t193 = sin(t196);
t240 = t193 * t197;
t262 = -t266 * t154 + t268 * t155 - t265 * t240;
t261 = t267 * t154 + t269 * t155 - t266 * t240;
t156 = t194 * t236 - t237;
t157 = t194 * t235 + t238;
t239 = t193 * t198;
t260 = t267 * t156 + t269 * t157 - t266 * t239;
t257 = -t266 * t156 + t268 * t157 - t265 * t239;
t256 = t269 * t154 + t270 * t155 + t268 * t240;
t255 = t269 * t156 + t270 * t157 + t268 * t239;
t254 = t266 * t194 + (t267 * t199 + t269 * t201) * t193;
t253 = t265 * t194 + (-t266 * t199 + t268 * t201) * t193;
t252 = -t268 * t194 + (t269 * t199 + t270 * t201) * t193;
t200 = sin(qJ(2));
t248 = pkin(2) * t200;
t202 = cos(qJ(2));
t247 = pkin(2) * t202;
t245 = Icges(2,4) * t197;
t244 = Icges(3,4) * t200;
t243 = Icges(3,4) * t202;
t242 = Icges(4,4) * t193;
t241 = Icges(4,4) * t194;
t234 = rSges(6,2) * t240 + t263 * t154 + t264 * t155;
t233 = rSges(6,2) * t239 + t263 * t156 + t264 * t157;
t120 = -pkin(6) * t198 + t197 * t247;
t180 = t197 * pkin(1) - t198 * pkin(5);
t232 = -t120 - t180;
t231 = -rSges(6,2) * t194 + (t263 * t199 + t264 * t201) * t193;
t230 = qJD(4) * t193;
t229 = V_base(5) * qJ(1) + V_base(1);
t225 = qJD(1) + V_base(3);
t187 = qJD(2) * t197 + V_base(4);
t186 = -qJD(2) * t198 + V_base(5);
t224 = t186 * t248 + t229;
t162 = qJD(3) * t197 + t187;
t223 = pkin(3) * t194 + pkin(7) * t193;
t222 = rSges(3,1) * t202 - rSges(3,2) * t200;
t221 = rSges(4,1) * t194 - rSges(4,2) * t193;
t220 = Icges(3,1) * t202 - t244;
t219 = Icges(4,1) * t194 - t242;
t218 = -Icges(3,2) * t200 + t243;
t217 = -Icges(4,2) * t193 + t241;
t216 = Icges(3,5) * t202 - Icges(3,6) * t200;
t215 = Icges(4,5) * t194 - Icges(4,6) * t193;
t181 = t198 * pkin(1) + t197 * pkin(5);
t214 = -V_base(4) * qJ(1) + V_base(6) * t181 + V_base(2);
t161 = V_base(5) + (-qJD(2) - qJD(3)) * t198;
t213 = V_base(4) * t180 - t181 * V_base(5) + t225;
t212 = (-Icges(4,3) * t198 + t197 * t215) * t161 + (Icges(4,3) * t197 + t198 * t215) * t162 + (Icges(4,5) * t193 + Icges(4,6) * t194) * V_base(6);
t211 = (-Icges(3,3) * t198 + t197 * t216) * t186 + (Icges(3,3) * t197 + t198 * t216) * t187 + (Icges(3,5) * t200 + Icges(3,6) * t202) * V_base(6);
t121 = pkin(6) * t197 + t198 * t247;
t210 = V_base(6) * t121 - t187 * t248 + t214;
t151 = t223 * t197;
t166 = pkin(3) * t193 - pkin(7) * t194;
t209 = t161 * t166 + (-t151 + t232) * V_base(6) + t224;
t208 = t187 * t120 - t121 * t186 + t213;
t152 = t223 * t198;
t207 = V_base(6) * t152 - t162 * t166 + t210;
t206 = t162 * t151 - t152 * t161 + t208;
t134 = -Icges(4,6) * t198 + t197 * t217;
t135 = Icges(4,6) * t197 + t198 * t217;
t136 = -Icges(4,5) * t198 + t197 * t219;
t137 = Icges(4,5) * t197 + t198 * t219;
t159 = Icges(4,2) * t194 + t242;
t160 = Icges(4,1) * t193 + t241;
t205 = (-t135 * t193 + t137 * t194) * t162 + (-t134 * t193 + t136 * t194) * t161 + (-t159 * t193 + t160 * t194) * V_base(6);
t145 = -Icges(3,6) * t198 + t197 * t218;
t146 = Icges(3,6) * t197 + t198 * t218;
t147 = -Icges(3,5) * t198 + t197 * t220;
t148 = Icges(3,5) * t197 + t198 * t220;
t183 = Icges(3,2) * t202 + t244;
t184 = Icges(3,1) * t200 + t243;
t204 = (-t146 * t200 + t148 * t202) * t187 + (-t145 * t200 + t147 * t202) * t186 + (-t183 * t200 + t184 * t202) * V_base(6);
t192 = Icges(2,4) * t198;
t185 = rSges(3,1) * t200 + rSges(3,2) * t202;
t179 = -qJD(4) * t194 + V_base(6);
t178 = rSges(2,1) * t198 - rSges(2,2) * t197;
t177 = rSges(2,1) * t197 + rSges(2,2) * t198;
t176 = Icges(2,1) * t198 - t245;
t175 = Icges(2,1) * t197 + t192;
t174 = -Icges(2,2) * t197 + t192;
t173 = Icges(2,2) * t198 + t245;
t170 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t169 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t168 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t163 = rSges(4,1) * t193 + rSges(4,2) * t194;
t150 = rSges(3,3) * t197 + t198 * t222;
t149 = -rSges(3,3) * t198 + t197 * t222;
t141 = t198 * t230 + t162;
t140 = t197 * t230 + t161;
t139 = rSges(4,3) * t197 + t198 * t221;
t138 = -rSges(4,3) * t198 + t197 * t221;
t131 = -rSges(5,3) * t194 + (rSges(5,1) * t201 - rSges(5,2) * t199) * t193;
t129 = V_base(5) * rSges(2,3) - t177 * V_base(6) + t229;
t128 = t178 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t117 = t177 * V_base(4) - t178 * V_base(5) + t225;
t112 = rSges(5,1) * t157 - rSges(5,2) * t156 + rSges(5,3) * t239;
t110 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t240;
t96 = t185 * t186 + (-t149 - t180) * V_base(6) + t229;
t95 = t150 * V_base(6) - t185 * t187 + t214;
t94 = t149 * t187 - t150 * t186 + t213;
t93 = t161 * t163 + (-t138 + t232) * V_base(6) + t224;
t92 = t139 * V_base(6) - t162 * t163 + t210;
t91 = t138 * t162 - t139 * t161 + t208;
t90 = -t110 * t179 + t131 * t140 + t209;
t89 = t112 * t179 - t131 * t141 + t207;
t88 = t110 * t141 - t112 * t140 + t206;
t87 = qJD(5) * t156 + t140 * t231 - t179 * t234 + t209;
t86 = qJD(5) * t154 - t141 * t231 + t179 * t233 + t207;
t85 = qJD(5) * t193 * t199 - t140 * t233 + t141 * t234 + t206;
t1 = m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t187 * (t211 * t197 + t204 * t198) / 0.2e1 + t186 * (t204 * t197 - t211 * t198) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t162 * (t212 * t197 + t205 * t198) / 0.2e1 + t161 * (t197 * t205 - t198 * t212) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((t154 * t254 + t155 * t252 + t240 * t253) * t179 + (t154 * t260 + t155 * t255 + t240 * t257) * t141 + (t261 * t154 + t256 * t155 + t262 * t240) * t140) * t140 / 0.2e1 + ((t156 * t254 + t157 * t252 + t239 * t253) * t179 + (t260 * t156 + t255 * t157 + t257 * t239) * t141 + (t261 * t156 + t256 * t157 + t239 * t262) * t140) * t141 / 0.2e1 + ((-t140 * t262 - t257 * t141 - t253 * t179) * t194 + ((t199 * t254 + t201 * t252) * t179 + (t199 * t260 + t201 * t255) * t141 + (t199 * t261 + t201 * t256) * t140) * t193) * t179 / 0.2e1 + ((-t173 * t197 + t175 * t198 + Icges(1,4)) * V_base(5) + (-t174 * t197 + t176 * t198 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t173 * t198 + t175 * t197 + Icges(1,2)) * V_base(5) + (t174 * t198 + t176 * t197 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t146 * t202 + t148 * t200) * t187 + (t145 * t202 + t147 * t200) * t186 + (t135 * t194 + t137 * t193) * t162 + (t134 * t194 + t136 * t193) * t161 + (t159 * t194 + t160 * t193 + t183 * t202 + t184 * t200 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t198 - Icges(2,6) * t197 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t197 + Icges(2,6) * t198 + Icges(1,6));
T = t1;
