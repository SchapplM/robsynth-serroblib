% Calculate kinetic energy for
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:35
% EndTime: 2019-12-05 15:32:37
% DurationCPUTime: 2.36s
% Computational Cost: add. (1272->256), mult. (1615->352), div. (0->0), fcn. (1509->8), ass. (0->129)
t275 = Icges(5,1) + Icges(6,1);
t274 = Icges(5,4) + Icges(6,4);
t273 = -Icges(6,5) - Icges(5,5);
t272 = Icges(5,2) + Icges(6,2);
t271 = -Icges(6,6) - Icges(5,6);
t270 = Icges(3,3) + Icges(4,3);
t269 = -Icges(6,3) - Icges(5,3);
t192 = qJ(2) + pkin(8);
t188 = sin(t192);
t189 = cos(t192);
t198 = sin(qJ(2));
t200 = cos(qJ(2));
t268 = Icges(3,5) * t200 + Icges(4,5) * t189 - Icges(3,6) * t198 - Icges(4,6) * t188;
t194 = cos(pkin(7));
t199 = cos(qJ(4));
t232 = t194 * t199;
t193 = sin(pkin(7));
t197 = sin(qJ(4));
t235 = t193 * t197;
t152 = -t189 * t235 - t232;
t233 = t194 * t197;
t234 = t193 * t199;
t153 = t189 * t234 - t233;
t237 = t188 * t193;
t267 = -t271 * t152 - t273 * t153 - t269 * t237;
t154 = -t189 * t233 + t234;
t155 = t189 * t232 + t235;
t236 = t188 * t194;
t266 = -t271 * t154 - t273 * t155 - t269 * t236;
t263 = t272 * t152 + t274 * t153 - t271 * t237;
t262 = t272 * t154 + t274 * t155 - t271 * t236;
t261 = t274 * t152 + t275 * t153 - t273 * t237;
t260 = t274 * t154 + t275 * t155 - t273 * t236;
t259 = t269 * t189 + (t271 * t197 - t273 * t199) * t188;
t258 = t271 * t189 + (-t272 * t197 + t274 * t199) * t188;
t257 = t273 * t189 + (-t274 * t197 + t275 * t199) * t188;
t238 = Icges(4,4) * t189;
t215 = -Icges(4,2) * t188 + t238;
t130 = -Icges(4,6) * t194 + t193 * t215;
t131 = Icges(4,6) * t193 + t194 * t215;
t239 = Icges(4,4) * t188;
t217 = Icges(4,1) * t189 - t239;
t132 = -Icges(4,5) * t194 + t193 * t217;
t133 = Icges(4,5) * t193 + t194 * t217;
t240 = Icges(3,4) * t200;
t216 = -Icges(3,2) * t198 + t240;
t142 = -Icges(3,6) * t194 + t193 * t216;
t143 = Icges(3,6) * t193 + t194 * t216;
t241 = Icges(3,4) * t198;
t218 = Icges(3,1) * t200 - t241;
t144 = -Icges(3,5) * t194 + t193 * t218;
t145 = Icges(3,5) * t193 + t194 * t218;
t157 = Icges(4,2) * t189 + t239;
t158 = Icges(4,1) * t188 + t238;
t179 = Icges(3,2) * t200 + t241;
t180 = Icges(3,1) * t198 + t240;
t182 = -qJD(2) * t194 + V_base(5);
t183 = qJD(2) * t193 + V_base(4);
t256 = (-t157 * t188 + t158 * t189 - t179 * t198 + t180 * t200) * V_base(6) + (-t131 * t188 + t133 * t189 - t143 * t198 + t145 * t200) * t183 + (-t130 * t188 + t132 * t189 - t142 * t198 + t144 * t200) * t182;
t255 = (Icges(3,5) * t198 + Icges(4,5) * t188 + Icges(3,6) * t200 + Icges(4,6) * t189) * V_base(6) + (t270 * t193 + t268 * t194) * t183 + (t268 * t193 - t270 * t194) * t182;
t249 = pkin(2) * t198;
t248 = pkin(2) * t200;
t247 = pkin(4) * t199;
t210 = qJ(5) * t188 + t189 * t247;
t244 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t237 - pkin(4) * t233 + t193 * t210;
t243 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t236 + pkin(4) * t235 + t194 * t210;
t242 = Icges(2,4) * t193;
t231 = (-qJ(5) - rSges(6,3)) * t189 + (rSges(6,1) * t199 - rSges(6,2) * t197 + t247) * t188;
t126 = -qJ(3) * t194 + t193 * t248;
t176 = pkin(1) * t193 - pkin(5) * t194;
t230 = -t126 - t176;
t229 = qJD(4) * t188;
t228 = qJD(5) * t188;
t227 = V_base(5) * qJ(1) + V_base(1);
t223 = qJD(1) + V_base(3);
t222 = qJD(3) * t193 + t182 * t249 + t227;
t221 = pkin(3) * t189 + pkin(6) * t188;
t220 = rSges(3,1) * t200 - rSges(3,2) * t198;
t219 = rSges(4,1) * t189 - rSges(4,2) * t188;
t177 = pkin(1) * t194 + pkin(5) * t193;
t212 = -V_base(4) * qJ(1) + V_base(6) * t177 + V_base(2);
t211 = V_base(4) * t176 - t177 * V_base(5) + t223;
t209 = t183 * t126 + t211;
t127 = qJ(3) * t193 + t194 * t248;
t206 = -qJD(3) * t194 + V_base(6) * t127 + t212;
t148 = t221 * t193;
t160 = pkin(3) * t188 - pkin(6) * t189;
t205 = t182 * t160 + (-t148 + t230) * V_base(6) + t222;
t149 = t221 * t194;
t204 = t183 * t148 + (-t127 - t149) * t182 + t209;
t203 = V_base(6) * t149 + (-t160 - t249) * t183 + t206;
t190 = Icges(2,4) * t194;
t181 = t198 * rSges(3,1) + rSges(3,2) * t200;
t175 = rSges(2,1) * t194 - rSges(2,2) * t193;
t174 = rSges(2,1) * t193 + rSges(2,2) * t194;
t173 = -qJD(4) * t189 + V_base(6);
t172 = Icges(2,1) * t194 - t242;
t171 = Icges(2,1) * t193 + t190;
t170 = -Icges(2,2) * t193 + t190;
t169 = Icges(2,2) * t194 + t242;
t166 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t165 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t164 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t159 = rSges(4,1) * t188 + rSges(4,2) * t189;
t151 = t194 * t229 + t183;
t150 = t193 * t229 + t182;
t147 = t193 * rSges(3,3) + t194 * t220;
t146 = -t194 * rSges(3,3) + t193 * t220;
t137 = rSges(4,3) * t193 + t194 * t219;
t136 = -rSges(4,3) * t194 + t193 * t219;
t135 = V_base(5) * rSges(2,3) - t174 * V_base(6) + t227;
t134 = t175 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t125 = -rSges(5,3) * t189 + (rSges(5,1) * t199 - rSges(5,2) * t197) * t188;
t115 = t174 * V_base(4) - t175 * V_base(5) + t223;
t112 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t236;
t110 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t237;
t94 = t181 * t182 + (-t146 - t176) * V_base(6) + t227;
t93 = t147 * V_base(6) - t181 * t183 + t212;
t92 = t146 * t183 - t147 * t182 + t211;
t91 = t159 * t182 + (-t136 + t230) * V_base(6) + t222;
t90 = t137 * V_base(6) + (-t159 - t249) * t183 + t206;
t89 = t136 * t183 + (-t127 - t137) * t182 + t209;
t88 = -t110 * t173 + t125 * t150 + t205;
t87 = t112 * t173 - t125 * t151 + t203;
t86 = t110 * t151 - t112 * t150 + t204;
t85 = t150 * t231 - t173 * t244 + t194 * t228 + t205;
t84 = -t151 * t231 + t173 * t243 + t193 * t228 + t203;
t83 = -qJD(5) * t189 - t150 * t243 + t151 * t244 + t204;
t1 = m(1) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(2) * (t115 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t258 * t152 + t257 * t153 + t259 * t237) * t173 + (t262 * t152 + t260 * t153 + t266 * t237) * t151 + (t263 * t152 + t261 * t153 + t267 * t237) * t150) * t150 / 0.2e1 + ((t258 * t154 + t257 * t155 + t259 * t236) * t173 + (t262 * t154 + t260 * t155 + t266 * t236) * t151 + (t263 * t154 + t261 * t155 + t267 * t236) * t150) * t151 / 0.2e1 + ((-t267 * t150 - t266 * t151 - t259 * t173) * t189 + ((-t258 * t197 + t257 * t199) * t173 + (-t262 * t197 + t260 * t199) * t151 + (-t263 * t197 + t261 * t199) * t150) * t188) * t173 / 0.2e1 + (t256 * t193 - t255 * t194) * t182 / 0.2e1 + (t255 * t193 + t256 * t194) * t183 / 0.2e1 + ((-t169 * t193 + t171 * t194 + Icges(1,4)) * V_base(5) + (-t170 * t193 + t172 * t194 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t169 * t194 + t171 * t193 + Icges(1,2)) * V_base(5) + (t170 * t194 + t172 * t193 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t131 * t189 + t133 * t188 + t143 * t200 + t198 * t145) * t183 + (t130 * t189 + t132 * t188 + t142 * t200 + t198 * t144) * t182 + (t157 * t189 + t158 * t188 + t179 * t200 + t198 * t180 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t194 - Icges(2,6) * t193 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t193 + Icges(2,6) * t194 + Icges(1,6));
T = t1;
