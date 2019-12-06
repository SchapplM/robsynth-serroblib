% Calculate kinetic energy for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:07
% EndTime: 2019-12-05 15:14:09
% DurationCPUTime: 1.91s
% Computational Cost: add. (1420->310), mult. (1618->456), div. (0->0), fcn. (1512->10), ass. (0->154)
t212 = sin(pkin(8));
t214 = cos(pkin(8));
t271 = Icges(2,5) * t214 - Icges(2,6) * t212 + Icges(1,5);
t270 = Icges(2,5) * t212 + Icges(2,6) * t214 + Icges(1,6);
t211 = sin(pkin(9));
t269 = pkin(2) * t211;
t213 = cos(pkin(9));
t268 = pkin(2) * t213;
t217 = cos(qJ(4));
t267 = pkin(4) * t217;
t265 = Icges(2,4) * t212;
t264 = Icges(3,4) * t211;
t263 = Icges(3,4) * t213;
t209 = pkin(9) + qJ(3);
t203 = sin(t209);
t262 = Icges(4,4) * t203;
t204 = cos(t209);
t261 = Icges(4,4) * t204;
t260 = t203 * t212;
t259 = t203 * t214;
t210 = qJ(4) + qJ(5);
t206 = sin(t210);
t258 = t206 * t212;
t257 = t206 * t214;
t207 = cos(t210);
t256 = t207 * t212;
t255 = t207 * t214;
t216 = sin(qJ(4));
t254 = t212 * t216;
t253 = t212 * t217;
t252 = t214 * t216;
t251 = t214 * t217;
t130 = -pkin(5) * t214 + t212 * t268;
t190 = pkin(1) * t212 - qJ(2) * t214;
t249 = -t130 - t190;
t248 = qJD(4) * t203;
t247 = qJD(5) * t203;
t246 = V_base(5) * qJ(1) + V_base(1);
t242 = qJD(1) + V_base(3);
t197 = qJD(3) * t212 + V_base(4);
t241 = qJD(2) * t212 + t246;
t240 = V_base(4) * t190 + t242;
t163 = t214 * t248 + t197;
t239 = V_base(5) * t269 + t241;
t238 = pkin(3) * t204 + pkin(6) * t203;
t196 = -qJD(3) * t214 + V_base(5);
t237 = rSges(3,1) * t213 - rSges(3,2) * t211;
t236 = rSges(4,1) * t204 - rSges(4,2) * t203;
t235 = Icges(3,1) * t213 - t264;
t234 = Icges(4,1) * t204 - t262;
t233 = -Icges(3,2) * t211 + t263;
t232 = -Icges(4,2) * t203 + t261;
t231 = Icges(3,5) * t213 - Icges(3,6) * t211;
t230 = Icges(4,5) * t204 - Icges(4,6) * t203;
t192 = pkin(1) * t214 + qJ(2) * t212;
t229 = -qJD(2) * t214 + V_base(6) * t192 + V_base(2);
t162 = t212 * t248 + t196;
t228 = pkin(7) * t203 + t204 * t267;
t227 = (-Icges(4,3) * t214 + t212 * t230) * t196 + (Icges(4,3) * t212 + t214 * t230) * t197 + (Icges(4,5) * t203 + Icges(4,6) * t204) * V_base(6);
t131 = pkin(5) * t212 + t214 * t268;
t226 = V_base(4) * t130 + (-t131 - t192) * V_base(5) + t240;
t225 = (-Icges(3,3) * t214 + t212 * t231) * V_base(5) + (Icges(3,3) * t212 + t214 * t231) * V_base(4) + (Icges(3,5) * t211 + Icges(3,6) * t213) * V_base(6);
t160 = t238 * t212;
t173 = t203 * pkin(3) - t204 * pkin(6);
t224 = t196 * t173 + (-t160 + t249) * V_base(6) + t239;
t223 = V_base(6) * t131 + (-qJ(1) - t269) * V_base(4) + t229;
t161 = t238 * t214;
t222 = t197 * t160 - t161 * t196 + t226;
t221 = V_base(6) * t161 - t173 * t197 + t223;
t138 = -Icges(4,6) * t214 + t212 * t232;
t139 = Icges(4,6) * t212 + t214 * t232;
t140 = -Icges(4,5) * t214 + t212 * t234;
t141 = Icges(4,5) * t212 + t214 * t234;
t170 = Icges(4,2) * t204 + t262;
t171 = Icges(4,1) * t203 + t261;
t220 = (-t139 * t203 + t141 * t204) * t197 + (-t138 * t203 + t140 * t204) * t196 + (-t170 * t203 + t171 * t204) * V_base(6);
t154 = -Icges(3,6) * t214 + t212 * t233;
t155 = Icges(3,6) * t212 + t214 * t233;
t156 = -Icges(3,5) * t214 + t212 * t235;
t157 = Icges(3,5) * t212 + t214 * t235;
t182 = Icges(3,2) * t213 + t264;
t185 = Icges(3,1) * t211 + t263;
t219 = (-t155 * t211 + t157 * t213) * V_base(4) + (-t154 * t211 + t156 * t213) * V_base(5) + (-t182 * t211 + t185 * t213) * V_base(6);
t205 = Icges(2,4) * t214;
t193 = rSges(2,1) * t214 - rSges(2,2) * t212;
t191 = rSges(2,1) * t212 + rSges(2,2) * t214;
t189 = rSges(3,1) * t211 + rSges(3,2) * t213;
t188 = -qJD(4) * t204 + V_base(6);
t187 = Icges(2,1) * t214 - t265;
t186 = Icges(2,1) * t212 + t205;
t184 = -Icges(2,2) * t212 + t205;
t183 = Icges(2,2) * t214 + t265;
t178 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t177 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t176 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = rSges(4,1) * t203 + rSges(4,2) * t204;
t168 = V_base(6) + (-qJD(4) - qJD(5)) * t204;
t167 = t204 * t251 + t254;
t166 = -t204 * t252 + t253;
t165 = t204 * t253 - t252;
t164 = -t204 * t254 - t251;
t159 = rSges(3,3) * t212 + t214 * t237;
t158 = -rSges(3,3) * t214 + t212 * t237;
t151 = t204 * t255 + t258;
t150 = -t204 * t257 + t256;
t149 = t204 * t256 - t257;
t148 = -t204 * t258 - t255;
t145 = rSges(4,3) * t212 + t214 * t236;
t144 = -rSges(4,3) * t214 + t212 * t236;
t143 = V_base(5) * rSges(2,3) - t191 * V_base(6) + t246;
t142 = t193 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t135 = -rSges(5,3) * t204 + (rSges(5,1) * t217 - rSges(5,2) * t216) * t203;
t134 = -Icges(5,5) * t204 + (Icges(5,1) * t217 - Icges(5,4) * t216) * t203;
t133 = -Icges(5,6) * t204 + (Icges(5,4) * t217 - Icges(5,2) * t216) * t203;
t132 = -Icges(5,3) * t204 + (Icges(5,5) * t217 - Icges(5,6) * t216) * t203;
t129 = t214 * t247 + t163;
t128 = t212 * t247 + t162;
t124 = -rSges(6,3) * t204 + (rSges(6,1) * t207 - rSges(6,2) * t206) * t203;
t123 = -Icges(6,5) * t204 + (Icges(6,1) * t207 - Icges(6,4) * t206) * t203;
t122 = -Icges(6,6) * t204 + (Icges(6,4) * t207 - Icges(6,2) * t206) * t203;
t121 = -Icges(6,3) * t204 + (Icges(6,5) * t207 - Icges(6,6) * t206) * t203;
t120 = t191 * V_base(4) - t193 * V_base(5) + t242;
t119 = -pkin(7) * t204 + t203 * t267;
t118 = rSges(5,1) * t167 + rSges(5,2) * t166 + rSges(5,3) * t259;
t117 = rSges(5,1) * t165 + rSges(5,2) * t164 + rSges(5,3) * t260;
t116 = Icges(5,1) * t167 + Icges(5,4) * t166 + Icges(5,5) * t259;
t115 = Icges(5,1) * t165 + Icges(5,4) * t164 + Icges(5,5) * t260;
t114 = Icges(5,4) * t167 + Icges(5,2) * t166 + Icges(5,6) * t259;
t113 = Icges(5,4) * t165 + Icges(5,2) * t164 + Icges(5,6) * t260;
t112 = Icges(5,5) * t167 + Icges(5,6) * t166 + Icges(5,3) * t259;
t111 = Icges(5,5) * t165 + Icges(5,6) * t164 + Icges(5,3) * t260;
t110 = pkin(4) * t254 + t214 * t228;
t109 = -pkin(4) * t252 + t212 * t228;
t108 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t259;
t107 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t260;
t106 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t259;
t105 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t260;
t104 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t259;
t103 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t260;
t102 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t259;
t101 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t260;
t100 = t189 * V_base(5) + (-t158 - t190) * V_base(6) + t241;
t99 = t159 * V_base(6) + (-qJ(1) - t189) * V_base(4) + t229;
t98 = t158 * V_base(4) + (-t159 - t192) * V_base(5) + t240;
t97 = t172 * t196 + (-t144 + t249) * V_base(6) + t239;
t96 = t145 * V_base(6) - t172 * t197 + t223;
t95 = t144 * t197 - t145 * t196 + t226;
t94 = -t117 * t188 + t135 * t162 + t224;
t93 = t118 * t188 - t135 * t163 + t221;
t92 = t117 * t163 - t118 * t162 + t222;
t91 = -t107 * t168 - t109 * t188 + t119 * t162 + t124 * t128 + t224;
t90 = t108 * t168 + t110 * t188 - t119 * t163 - t124 * t129 + t221;
t89 = t107 * t129 - t108 * t128 + t109 * t163 - t110 * t162 + t222;
t1 = m(1) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(2) * (t120 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + t197 * (t212 * t227 + t214 * t220) / 0.2e1 + t196 * (t212 * t220 - t214 * t227) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t163 * ((t112 * t259 + t166 * t114 + t167 * t116) * t163 + (t111 * t259 + t113 * t166 + t115 * t167) * t162 + (t132 * t259 + t133 * t166 + t134 * t167) * t188) / 0.2e1 + t162 * ((t112 * t260 + t114 * t164 + t116 * t165) * t163 + (t111 * t260 + t164 * t113 + t165 * t115) * t162 + (t132 * t260 + t133 * t164 + t134 * t165) * t188) / 0.2e1 + t188 * ((-t111 * t162 - t112 * t163 - t132 * t188) * t204 + ((-t114 * t216 + t116 * t217) * t163 + (-t113 * t216 + t115 * t217) * t162 + (-t133 * t216 + t134 * t217) * t188) * t203) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t129 * ((t102 * t259 + t150 * t104 + t151 * t106) * t129 + (t101 * t259 + t103 * t150 + t105 * t151) * t128 + (t121 * t259 + t122 * t150 + t123 * t151) * t168) / 0.2e1 + t128 * ((t102 * t260 + t104 * t148 + t106 * t149) * t129 + (t101 * t260 + t148 * t103 + t149 * t105) * t128 + (t121 * t260 + t122 * t148 + t123 * t149) * t168) / 0.2e1 + t168 * ((-t101 * t128 - t102 * t129 - t121 * t168) * t204 + ((-t104 * t206 + t106 * t207) * t129 + (-t103 * t206 + t105 * t207) * t128 + (-t122 * t206 + t123 * t207) * t168) * t203) / 0.2e1 + (t212 * t225 + t214 * t219 + t271 * V_base(6) + (-t183 * t212 + t186 * t214 + Icges(1,4)) * V_base(5) + (-t212 * t184 + t214 * t187 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t212 * t219 - t214 * t225 + t270 * V_base(6) + (t214 * t183 + t212 * t186 + Icges(1,2)) * V_base(5) + (t184 * t214 + t187 * t212 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t139 * t204 + t141 * t203) * t197 + (t138 * t204 + t140 * t203) * t196 + (t154 * t213 + t156 * t211 + t270) * V_base(5) + (t155 * t213 + t157 * t211 + t271) * V_base(4) + (t204 * t170 + t203 * t171 + t213 * t182 + t211 * t185 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1;
T = t1;
