% Calculate kinetic energy for
% S5RRPRP7
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:51
% EndTime: 2019-12-31 19:59:53
% DurationCPUTime: 2.25s
% Computational Cost: add. (1287->251), mult. (1606->350), div. (0->0), fcn. (1510->8), ass. (0->128)
t272 = Icges(5,1) + Icges(6,1);
t271 = -Icges(5,4) + Icges(6,5);
t270 = Icges(6,4) + Icges(5,5);
t269 = Icges(5,2) + Icges(6,3);
t268 = -Icges(6,6) + Icges(5,6);
t267 = Icges(3,3) + Icges(4,3);
t266 = -Icges(5,3) - Icges(6,2);
t195 = qJ(2) + pkin(8);
t188 = sin(t195);
t189 = cos(t195);
t198 = sin(qJ(2));
t201 = cos(qJ(2));
t265 = Icges(3,5) * t201 + Icges(4,5) * t189 - Icges(3,6) * t198 - Icges(4,6) * t188;
t264 = rSges(6,1) + pkin(4);
t263 = rSges(6,3) + qJ(5);
t200 = cos(qJ(4));
t202 = cos(qJ(1));
t232 = t200 * t202;
t197 = sin(qJ(4));
t199 = sin(qJ(1));
t234 = t199 * t197;
t154 = t189 * t234 + t232;
t233 = t199 * t200;
t235 = t197 * t202;
t155 = t189 * t233 - t235;
t237 = t188 * t199;
t262 = -t268 * t154 + t270 * t155 - t266 * t237;
t261 = t269 * t154 + t271 * t155 - t268 * t237;
t156 = t189 * t235 - t233;
t157 = t189 * t232 + t234;
t236 = t188 * t202;
t260 = t269 * t156 + t271 * t157 - t268 * t236;
t259 = -t268 * t156 + t270 * t157 - t266 * t236;
t258 = t271 * t154 + t272 * t155 + t270 * t237;
t257 = t271 * t156 + t272 * t157 + t270 * t236;
t256 = t268 * t189 + (t269 * t197 + t271 * t200) * t188;
t255 = t266 * t189 + (-t268 * t197 + t270 * t200) * t188;
t254 = -t270 * t189 + (t271 * t197 + t272 * t200) * t188;
t238 = Icges(4,4) * t189;
t216 = -Icges(4,2) * t188 + t238;
t134 = -Icges(4,6) * t202 + t199 * t216;
t135 = Icges(4,6) * t199 + t202 * t216;
t239 = Icges(4,4) * t188;
t218 = Icges(4,1) * t189 - t239;
t136 = -Icges(4,5) * t202 + t199 * t218;
t137 = Icges(4,5) * t199 + t202 * t218;
t240 = Icges(3,4) * t201;
t217 = -Icges(3,2) * t198 + t240;
t143 = -Icges(3,6) * t202 + t199 * t217;
t144 = Icges(3,6) * t199 + t202 * t217;
t241 = Icges(3,4) * t198;
t219 = Icges(3,1) * t201 - t241;
t145 = -Icges(3,5) * t202 + t199 * t219;
t146 = Icges(3,5) * t199 + t202 * t219;
t160 = Icges(4,2) * t189 + t239;
t161 = Icges(4,1) * t188 + t238;
t173 = Icges(3,2) * t201 + t241;
t176 = Icges(3,1) * t198 + t240;
t184 = -qJD(2) * t202 + V_base(5);
t185 = qJD(2) * t199 + V_base(4);
t190 = V_base(6) + qJD(1);
t253 = (-t160 * t188 + t161 * t189 - t173 * t198 + t176 * t201) * t190 + (-t135 * t188 + t137 * t189 - t144 * t198 + t146 * t201) * t185 + (-t134 * t188 + t136 * t189 - t143 * t198 + t145 * t201) * t184;
t252 = (Icges(3,5) * t198 + Icges(4,5) * t188 + Icges(3,6) * t201 + Icges(4,6) * t189) * t190 + (t267 * t199 + t265 * t202) * t185 + (t265 * t199 - t267 * t202) * t184;
t245 = pkin(2) * t198;
t244 = pkin(2) * t201;
t242 = Icges(2,4) * t199;
t231 = rSges(6,2) * t237 + t263 * t154 + t264 * t155;
t230 = rSges(6,2) * t236 + t263 * t156 + t264 * t157;
t229 = -rSges(6,2) * t189 + (t263 * t197 + t264 * t200) * t188;
t129 = -qJ(3) * t202 + t199 * t244;
t182 = t199 * pkin(1) - pkin(6) * t202;
t228 = -t129 - t182;
t227 = qJD(4) * t188;
t226 = V_base(5) * pkin(5) + V_base(1);
t223 = qJD(3) * t199 + t184 * t245 + t226;
t222 = pkin(3) * t189 + pkin(7) * t188;
t221 = rSges(3,1) * t201 - rSges(3,2) * t198;
t220 = rSges(4,1) * t189 - rSges(4,2) * t188;
t183 = pkin(1) * t202 + t199 * pkin(6);
t213 = -V_base(4) * pkin(5) + t190 * t183 + V_base(2);
t212 = V_base(4) * t182 - t183 * V_base(5) + V_base(3);
t211 = t185 * t129 + t212;
t130 = qJ(3) * t199 + t202 * t244;
t208 = -qJD(3) * t202 + t190 * t130 + t213;
t150 = t222 * t199;
t163 = pkin(3) * t188 - pkin(7) * t189;
t207 = t184 * t163 + (-t150 + t228) * t190 + t223;
t151 = t222 * t202;
t206 = t185 * t150 + (-t130 - t151) * t184 + t211;
t205 = t190 * t151 + (-t163 - t245) * t185 + t208;
t193 = Icges(2,4) * t202;
t181 = rSges(2,1) * t202 - t199 * rSges(2,2);
t180 = t199 * rSges(2,1) + rSges(2,2) * t202;
t179 = rSges(3,1) * t198 + rSges(3,2) * t201;
t178 = Icges(2,1) * t202 - t242;
t177 = Icges(2,1) * t199 + t193;
t175 = -Icges(2,2) * t199 + t193;
t174 = Icges(2,2) * t202 + t242;
t169 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t168 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t167 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t166 = -qJD(4) * t189 + t190;
t162 = rSges(4,1) * t188 + rSges(4,2) * t189;
t153 = t202 * t227 + t185;
t152 = t199 * t227 + t184;
t149 = t199 * rSges(3,3) + t202 * t221;
t148 = -rSges(3,3) * t202 + t199 * t221;
t139 = t199 * rSges(4,3) + t202 * t220;
t138 = -rSges(4,3) * t202 + t199 * t220;
t128 = -rSges(5,3) * t189 + (rSges(5,1) * t200 - rSges(5,2) * t197) * t188;
t126 = V_base(5) * rSges(2,3) - t180 * t190 + t226;
t125 = t181 * t190 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t117 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t112 = t157 * rSges(5,1) - t156 * rSges(5,2) + rSges(5,3) * t236;
t110 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t237;
t96 = t179 * t184 + (-t148 - t182) * t190 + t226;
t95 = t149 * t190 - t179 * t185 + t213;
t94 = t148 * t185 - t149 * t184 + t212;
t93 = t162 * t184 + (-t138 + t228) * t190 + t223;
t92 = t190 * t139 + (-t162 - t245) * t185 + t208;
t91 = t138 * t185 + (-t130 - t139) * t184 + t211;
t90 = -t110 * t166 + t128 * t152 + t207;
t89 = t166 * t112 - t153 * t128 + t205;
t88 = t110 * t153 - t112 * t152 + t206;
t87 = qJD(5) * t156 + t152 * t229 - t166 * t231 + t207;
t86 = qJD(5) * t154 - t153 * t229 + t166 * t230 + t205;
t85 = qJD(5) * t188 * t197 - t152 * t230 + t153 * t231 + t206;
t1 = m(1) * (t167 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((t256 * t154 + t254 * t155 + t255 * t237) * t166 + (t260 * t154 + t257 * t155 + t259 * t237) * t153 + (t261 * t154 + t258 * t155 + t262 * t237) * t152) * t152 / 0.2e1 + ((t256 * t156 + t254 * t157 + t255 * t236) * t166 + (t260 * t156 + t257 * t157 + t259 * t236) * t153 + (t261 * t156 + t258 * t157 + t262 * t236) * t152) * t153 / 0.2e1 + ((-t262 * t152 - t259 * t153 - t255 * t166) * t189 + ((t256 * t197 + t254 * t200) * t166 + (t260 * t197 + t257 * t200) * t153 + (t261 * t197 + t258 * t200) * t152) * t188) * t166 / 0.2e1 + (t253 * t199 - t252 * t202) * t184 / 0.2e1 + (t252 * t199 + t253 * t202) * t185 / 0.2e1 + ((-t199 * t174 + t177 * t202 + Icges(1,4)) * V_base(5) + (-t199 * t175 + t178 * t202 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t174 * t202 + t199 * t177 + Icges(1,2)) * V_base(5) + (t175 * t202 + t199 * t178 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t135 * t189 + t137 * t188 + t144 * t201 + t146 * t198) * t185 + (t134 * t189 + t136 * t188 + t143 * t201 + t145 * t198) * t184 + (t160 * t189 + t161 * t188 + t173 * t201 + t176 * t198 + Icges(2,3)) * t190) * t190 / 0.2e1 + t190 * V_base(4) * (Icges(2,5) * t202 - Icges(2,6) * t199) + V_base(5) * t190 * (Icges(2,5) * t199 + Icges(2,6) * t202) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
