% Calculate kinetic energy for
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:23
% EndTime: 2019-12-31 18:53:25
% DurationCPUTime: 2.36s
% Computational Cost: add. (1265->257), mult. (1584->364), div. (0->0), fcn. (1488->8), ass. (0->133)
t267 = Icges(5,1) + Icges(6,1);
t266 = -Icges(5,4) + Icges(6,5);
t265 = Icges(6,4) + Icges(5,5);
t264 = Icges(5,2) + Icges(6,3);
t263 = -Icges(6,6) + Icges(5,6);
t262 = -Icges(5,3) - Icges(6,2);
t261 = rSges(6,1) + pkin(4);
t260 = rSges(6,3) + qJ(5);
t195 = pkin(8) + qJ(3);
t189 = cos(t195);
t201 = cos(qJ(4));
t202 = cos(qJ(1));
t234 = t201 * t202;
t199 = sin(qJ(4));
t200 = sin(qJ(1));
t236 = t200 * t199;
t154 = t189 * t236 + t234;
t235 = t200 * t201;
t237 = t199 * t202;
t155 = t189 * t235 - t237;
t188 = sin(t195);
t239 = t188 * t200;
t259 = -t263 * t154 + t265 * t155 - t262 * t239;
t258 = t264 * t154 + t266 * t155 - t263 * t239;
t156 = t189 * t237 - t235;
t157 = t189 * t234 + t236;
t238 = t188 * t202;
t257 = t264 * t156 + t266 * t157 - t263 * t238;
t256 = -t263 * t156 + t265 * t157 - t262 * t238;
t255 = t266 * t154 + t267 * t155 + t265 * t239;
t254 = t266 * t156 + t267 * t157 + t265 * t238;
t253 = t263 * t189 + (t264 * t199 + t266 * t201) * t188;
t252 = t262 * t189 + (-t263 * t199 + t265 * t201) * t188;
t251 = -t265 * t189 + (t266 * t199 + t267 * t201) * t188;
t196 = sin(pkin(8));
t246 = pkin(2) * t196;
t197 = cos(pkin(8));
t245 = pkin(2) * t197;
t244 = Icges(2,4) * t200;
t243 = Icges(3,4) * t196;
t242 = Icges(3,4) * t197;
t241 = Icges(4,4) * t188;
t240 = Icges(4,4) * t189;
t232 = rSges(6,2) * t239 + t260 * t154 + t261 * t155;
t231 = rSges(6,2) * t238 + t260 * t156 + t261 * t157;
t230 = -rSges(6,2) * t189 + (t260 * t199 + t261 * t201) * t188;
t129 = -pkin(6) * t202 + t200 * t245;
t179 = t200 * pkin(1) - qJ(2) * t202;
t229 = -t129 - t179;
t228 = qJD(4) * t188;
t227 = V_base(4) * t179 + V_base(3);
t226 = V_base(5) * pkin(5) + V_base(1);
t184 = qJD(3) * t200 + V_base(4);
t190 = V_base(6) + qJD(1);
t223 = qJD(2) * t200 + t226;
t222 = V_base(5) * t246 + t223;
t221 = pkin(3) * t189 + pkin(7) * t188;
t183 = -qJD(3) * t202 + V_base(5);
t220 = rSges(3,1) * t197 - rSges(3,2) * t196;
t219 = rSges(4,1) * t189 - rSges(4,2) * t188;
t218 = Icges(3,1) * t197 - t243;
t217 = Icges(4,1) * t189 - t241;
t216 = -Icges(3,2) * t196 + t242;
t215 = -Icges(4,2) * t188 + t240;
t214 = Icges(3,5) * t197 - Icges(3,6) * t196;
t213 = Icges(4,5) * t189 - Icges(4,6) * t188;
t181 = pkin(1) * t202 + t200 * qJ(2);
t212 = -qJD(2) * t202 + t190 * t181 + V_base(2);
t211 = (-Icges(4,3) * t202 + t200 * t213) * t183 + (Icges(4,3) * t200 + t202 * t213) * t184 + (Icges(4,5) * t188 + Icges(4,6) * t189) * t190;
t130 = pkin(6) * t200 + t202 * t245;
t210 = V_base(4) * t129 + (-t130 - t181) * V_base(5) + t227;
t209 = (-Icges(3,3) * t202 + t200 * t214) * V_base(5) + (Icges(3,3) * t200 + t202 * t214) * V_base(4) + (Icges(3,5) * t196 + Icges(3,6) * t197) * t190;
t150 = t221 * t200;
t163 = pkin(3) * t188 - pkin(7) * t189;
t208 = t183 * t163 + (-t150 + t229) * t190 + t222;
t151 = t221 * t202;
t207 = t184 * t150 - t151 * t183 + t210;
t206 = t190 * t130 + (-pkin(5) - t246) * V_base(4) + t212;
t205 = t190 * t151 - t184 * t163 + t206;
t134 = -Icges(4,6) * t202 + t200 * t215;
t135 = Icges(4,6) * t200 + t202 * t215;
t136 = -Icges(4,5) * t202 + t200 * t217;
t137 = Icges(4,5) * t200 + t202 * t217;
t160 = Icges(4,2) * t189 + t241;
t161 = Icges(4,1) * t188 + t240;
t204 = (-t135 * t188 + t137 * t189) * t184 + (-t134 * t188 + t136 * t189) * t183 + (-t160 * t188 + t161 * t189) * t190;
t143 = -Icges(3,6) * t202 + t200 * t216;
t144 = Icges(3,6) * t200 + t202 * t216;
t145 = -Icges(3,5) * t202 + t200 * t218;
t146 = Icges(3,5) * t200 + t202 * t218;
t170 = Icges(3,2) * t197 + t243;
t171 = Icges(3,1) * t196 + t242;
t203 = (-t144 * t196 + t146 * t197) * V_base(4) + (-t143 * t196 + t145 * t197) * V_base(5) + (-t170 * t196 + t171 * t197) * t190;
t193 = Icges(2,4) * t202;
t182 = rSges(2,1) * t202 - t200 * rSges(2,2);
t180 = t200 * rSges(2,1) + rSges(2,2) * t202;
t178 = Icges(2,1) * t202 - t244;
t177 = Icges(2,1) * t200 + t193;
t176 = -Icges(2,2) * t200 + t193;
t175 = Icges(2,2) * t202 + t244;
t174 = Icges(2,5) * t202 - Icges(2,6) * t200;
t173 = Icges(2,5) * t200 + Icges(2,6) * t202;
t172 = rSges(3,1) * t196 + rSges(3,2) * t197;
t168 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t167 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t166 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t165 = -qJD(4) * t189 + t190;
t162 = rSges(4,1) * t188 + rSges(4,2) * t189;
t153 = t202 * t228 + t184;
t152 = t200 * t228 + t183;
t148 = t200 * rSges(3,3) + t202 * t220;
t147 = -rSges(3,3) * t202 + t200 * t220;
t139 = t200 * rSges(4,3) + t202 * t219;
t138 = -rSges(4,3) * t202 + t200 * t219;
t128 = -rSges(5,3) * t189 + (rSges(5,1) * t201 - rSges(5,2) * t199) * t188;
t126 = V_base(5) * rSges(2,3) - t180 * t190 + t226;
t125 = t182 * t190 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t117 = t180 * V_base(4) - t182 * V_base(5) + V_base(3);
t112 = t157 * rSges(5,1) - t156 * rSges(5,2) + rSges(5,3) * t238;
t110 = rSges(5,1) * t155 - rSges(5,2) * t154 + rSges(5,3) * t239;
t96 = t172 * V_base(5) + (-t147 - t179) * t190 + t223;
t95 = t190 * t148 + (-pkin(5) - t172) * V_base(4) + t212;
t94 = t147 * V_base(4) + (-t148 - t181) * V_base(5) + t227;
t93 = t162 * t183 + (-t138 + t229) * t190 + t222;
t92 = t190 * t139 - t184 * t162 + t206;
t91 = t138 * t184 - t139 * t183 + t210;
t90 = -t110 * t165 + t128 * t152 + t208;
t89 = t165 * t112 - t153 * t128 + t205;
t88 = t110 * t153 - t112 * t152 + t207;
t87 = qJD(5) * t156 + t152 * t230 - t165 * t232 + t208;
t86 = qJD(5) * t154 - t153 * t230 + t165 * t231 + t205;
t85 = qJD(5) * t188 * t199 - t152 * t231 + t153 * t232 + t207;
t1 = m(1) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t184 * (t200 * t211 + t202 * t204) / 0.2e1 + t183 * (t200 * t204 - t202 * t211) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((t154 * t253 + t155 * t251 + t239 * t252) * t165 + (t154 * t257 + t155 * t254 + t239 * t256) * t153 + (t258 * t154 + t255 * t155 + t259 * t239) * t152) * t152 / 0.2e1 + ((t156 * t253 + t157 * t251 + t238 * t252) * t165 + (t257 * t156 + t254 * t157 + t256 * t238) * t153 + (t258 * t156 + t255 * t157 + t238 * t259) * t152) * t153 / 0.2e1 + ((-t152 * t259 - t256 * t153 - t252 * t165) * t189 + ((t199 * t253 + t201 * t251) * t165 + (t199 * t257 + t201 * t254) * t153 + (t199 * t258 + t201 * t255) * t152) * t188) * t165 / 0.2e1 + ((t135 * t189 + t137 * t188) * t184 + (t134 * t189 + t136 * t188) * t183 + (t143 * t197 + t145 * t196 + t173) * V_base(5) + (t144 * t197 + t146 * t196 + t174) * V_base(4) + (t160 * t189 + t161 * t188 + t170 * t197 + t171 * t196 + Icges(2,3)) * t190) * t190 / 0.2e1 + (t174 * t190 + t200 * t209 + t202 * t203 + (-t200 * t175 + t177 * t202 + Icges(1,4)) * V_base(5) + (-t200 * t176 + t178 * t202 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t173 * t190 + t200 * t203 - t202 * t209 + (t175 * t202 + t200 * t177 + Icges(1,2)) * V_base(5) + (t176 * t202 + t200 * t178 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
