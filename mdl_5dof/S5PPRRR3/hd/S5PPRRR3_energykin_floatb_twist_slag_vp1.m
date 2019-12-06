% Calculate kinetic energy for
% S5PPRRR3
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:18
% EndTime: 2019-12-05 15:16:21
% DurationCPUTime: 2.41s
% Computational Cost: add. (1365->337), mult. (2654->494), div. (0->0), fcn. (2900->10), ass. (0->156)
t221 = sin(pkin(8));
t223 = cos(pkin(8));
t270 = Icges(2,5) * t223 - Icges(2,6) * t221 + Icges(1,5);
t269 = Icges(2,5) * t221 + Icges(2,6) * t223 + Icges(1,6);
t226 = cos(qJ(4));
t268 = pkin(4) * t226;
t266 = Icges(2,4) * t221;
t220 = sin(pkin(9));
t265 = Icges(3,4) * t220;
t222 = cos(pkin(9));
t264 = Icges(3,4) * t222;
t263 = t220 * t221;
t262 = t220 * t223;
t224 = sin(qJ(4));
t261 = t220 * t224;
t225 = sin(qJ(3));
t260 = t220 * t225;
t259 = t220 * t226;
t227 = cos(qJ(3));
t258 = t220 * t227;
t257 = t221 * t225;
t256 = t221 * t227;
t255 = t222 * t224;
t254 = t223 * t225;
t253 = t223 * t227;
t252 = qJD(3) * t220;
t251 = V_base(5) * qJ(1) + V_base(1);
t247 = qJD(1) + V_base(3);
t246 = t221 * t261;
t245 = t223 * t261;
t192 = t223 * t252 + V_base(4);
t191 = t221 * t252 + V_base(5);
t244 = qJD(2) * t221 + t251;
t203 = pkin(1) * t221 - qJ(2) * t223;
t243 = V_base(4) * t203 + t247;
t177 = t222 * t254 - t256;
t148 = qJD(4) * t177 + t192;
t175 = t222 * t257 + t253;
t147 = qJD(4) * t175 + t191;
t242 = pkin(2) * t222 + pkin(5) * t220;
t208 = -qJD(3) * t222 + V_base(6);
t241 = rSges(3,1) * t222 - rSges(3,2) * t220;
t240 = Icges(3,1) * t222 - t265;
t239 = -Icges(3,2) * t220 + t264;
t238 = Icges(3,5) * t222 - Icges(3,6) * t220;
t205 = pkin(1) * t223 + qJ(2) * t221;
t237 = -qJD(2) * t223 + V_base(6) * t205 + V_base(2);
t179 = qJD(4) * t260 + t208;
t182 = t242 * t221;
t207 = pkin(2) * t220 - pkin(5) * t222;
t236 = V_base(5) * t207 + (-t182 - t203) * V_base(6) + t244;
t183 = t242 * t223;
t235 = V_base(4) * t182 + (-t183 - t205) * V_base(5) + t243;
t234 = (-Icges(3,3) * t223 + t221 * t238) * V_base(5) + (Icges(3,3) * t221 + t223 * t238) * V_base(4) + (Icges(3,5) * t220 + Icges(3,6) * t222) * V_base(6);
t233 = V_base(6) * t183 + (-qJ(1) - t207) * V_base(4) + t237;
t176 = t222 * t256 - t254;
t140 = t176 * pkin(3) + t175 * pkin(6);
t184 = (pkin(3) * t227 + pkin(6) * t225) * t220;
t232 = -t140 * t208 + t191 * t184 + t236;
t178 = t222 * t253 + t257;
t141 = t178 * pkin(3) + t177 * pkin(6);
t231 = t192 * t140 - t141 * t191 + t235;
t230 = t208 * t141 - t184 * t192 + t233;
t159 = -Icges(3,6) * t223 + t221 * t239;
t160 = Icges(3,6) * t221 + t223 * t239;
t161 = -Icges(3,5) * t223 + t221 * t240;
t162 = Icges(3,5) * t221 + t223 * t240;
t196 = Icges(3,2) * t222 + t265;
t199 = Icges(3,1) * t220 + t264;
t229 = (-t160 * t220 + t162 * t222) * V_base(4) + (-t159 * t220 + t161 * t222) * V_base(5) + (-t196 * t220 + t199 * t222) * V_base(6);
t219 = qJ(4) + qJ(5);
t217 = cos(t219);
t216 = sin(t219);
t215 = Icges(2,4) * t223;
t206 = rSges(2,1) * t223 - rSges(2,2) * t221;
t204 = rSges(2,1) * t221 + rSges(2,2) * t223;
t202 = rSges(3,1) * t220 + rSges(3,2) * t222;
t201 = Icges(2,1) * t223 - t266;
t200 = Icges(2,1) * t221 + t215;
t198 = -Icges(2,2) * t221 + t215;
t197 = Icges(2,2) * t223 + t266;
t190 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t189 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t188 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t181 = t226 * t258 - t255;
t180 = -t222 * t226 - t224 * t258;
t170 = -t216 * t222 + t217 * t258;
t169 = -t216 * t258 - t217 * t222;
t168 = -rSges(4,3) * t222 + (rSges(4,1) * t227 - rSges(4,2) * t225) * t220;
t167 = -Icges(4,5) * t222 + (Icges(4,1) * t227 - Icges(4,4) * t225) * t220;
t166 = -Icges(4,6) * t222 + (Icges(4,4) * t227 - Icges(4,2) * t225) * t220;
t165 = -Icges(4,3) * t222 + (Icges(4,5) * t227 - Icges(4,6) * t225) * t220;
t164 = rSges(3,3) * t221 + t223 * t241;
t163 = -rSges(3,3) * t223 + t221 * t241;
t156 = qJD(5) * t260 + t179;
t155 = V_base(5) * rSges(2,3) - t204 * V_base(6) + t251;
t154 = t206 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t152 = t178 * t226 + t245;
t151 = -t178 * t224 + t223 * t259;
t150 = t176 * t226 + t246;
t149 = -t176 * t224 + t221 * t259;
t146 = t204 * V_base(4) - t206 * V_base(5) + t247;
t145 = t178 * t217 + t216 * t262;
t144 = -t178 * t216 + t217 * t262;
t143 = t176 * t217 + t216 * t263;
t142 = -t176 * t216 + t217 * t263;
t139 = -pkin(4) * t255 + (pkin(7) * t225 + t227 * t268) * t220;
t138 = rSges(5,1) * t181 + rSges(5,2) * t180 + rSges(5,3) * t260;
t137 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t260;
t136 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t260;
t135 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t260;
t134 = rSges(4,1) * t178 - rSges(4,2) * t177 + rSges(4,3) * t262;
t133 = rSges(4,1) * t176 - rSges(4,2) * t175 + rSges(4,3) * t263;
t132 = Icges(4,1) * t178 - Icges(4,4) * t177 + Icges(4,5) * t262;
t131 = Icges(4,1) * t176 - Icges(4,4) * t175 + Icges(4,5) * t263;
t130 = Icges(4,4) * t178 - Icges(4,2) * t177 + Icges(4,6) * t262;
t129 = Icges(4,4) * t176 - Icges(4,2) * t175 + Icges(4,6) * t263;
t128 = Icges(4,5) * t178 - Icges(4,6) * t177 + Icges(4,3) * t262;
t127 = Icges(4,5) * t176 - Icges(4,6) * t175 + Icges(4,3) * t263;
t125 = qJD(5) * t177 + t148;
t124 = qJD(5) * t175 + t147;
t122 = rSges(6,1) * t170 + rSges(6,2) * t169 + rSges(6,3) * t260;
t121 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t260;
t120 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t260;
t119 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t260;
t118 = t202 * V_base(5) + (-t163 - t203) * V_base(6) + t244;
t117 = t164 * V_base(6) + (-qJ(1) - t202) * V_base(4) + t237;
t116 = t163 * V_base(4) + (-t164 - t205) * V_base(5) + t243;
t115 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t177;
t114 = rSges(5,1) * t150 + rSges(5,2) * t149 + rSges(5,3) * t175;
t113 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t177;
t112 = Icges(5,1) * t150 + Icges(5,4) * t149 + Icges(5,5) * t175;
t111 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t177;
t110 = Icges(5,4) * t150 + Icges(5,2) * t149 + Icges(5,6) * t175;
t109 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t177;
t108 = Icges(5,5) * t150 + Icges(5,6) * t149 + Icges(5,3) * t175;
t107 = pkin(4) * t245 + pkin(7) * t177 + t178 * t268;
t106 = pkin(4) * t246 + pkin(7) * t175 + t176 * t268;
t105 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t177;
t104 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t175;
t103 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t177;
t102 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t175;
t101 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t177;
t100 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t175;
t99 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t177;
t98 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t175;
t97 = -t133 * t208 + t168 * t191 + t236;
t96 = t134 * t208 - t168 * t192 + t233;
t95 = t133 * t192 - t134 * t191 + t235;
t94 = -t114 * t179 + t138 * t147 + t232;
t93 = t115 * t179 - t138 * t148 + t230;
t92 = t114 * t148 - t115 * t147 + t231;
t91 = -t104 * t156 - t106 * t179 + t122 * t124 + t139 * t147 + t232;
t90 = t105 * t156 + t107 * t179 - t122 * t125 - t139 * t148 + t230;
t89 = t104 * t125 - t105 * t124 + t106 * t148 - t107 * t147 + t231;
t1 = m(1) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(2) * (t146 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(3) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + t192 * ((t128 * t262 - t130 * t177 + t132 * t178) * t192 + (t127 * t262 - t129 * t177 + t131 * t178) * t191 + (t165 * t262 - t166 * t177 + t167 * t178) * t208) / 0.2e1 + t191 * ((t128 * t263 - t130 * t175 + t132 * t176) * t192 + (t127 * t263 - t175 * t129 + t176 * t131) * t191 + (t165 * t263 - t166 * t175 + t167 * t176) * t208) / 0.2e1 + t208 * ((-t127 * t191 - t128 * t192 - t165 * t208) * t222 + ((-t130 * t225 + t132 * t227) * t192 + (-t129 * t225 + t131 * t227) * t191 + (-t166 * t225 + t167 * t227) * t208) * t220) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t148 * ((t177 * t109 + t151 * t111 + t152 * t113) * t148 + (t108 * t177 + t110 * t151 + t112 * t152) * t147 + (t135 * t177 + t136 * t151 + t137 * t152) * t179) / 0.2e1 + t147 * ((t109 * t175 + t111 * t149 + t113 * t150) * t148 + (t175 * t108 + t149 * t110 + t150 * t112) * t147 + (t135 * t175 + t136 * t149 + t137 * t150) * t179) / 0.2e1 + t179 * ((t109 * t260 + t111 * t180 + t113 * t181) * t148 + (t108 * t260 + t110 * t180 + t112 * t181) * t147 + (t135 * t260 + t180 * t136 + t181 * t137) * t179) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t125 * ((t144 * t101 + t145 * t103 + t177 * t99) * t125 + (t100 * t144 + t102 * t145 + t177 * t98) * t124 + (t119 * t177 + t120 * t144 + t121 * t145) * t156) / 0.2e1 + t124 * ((t101 * t142 + t103 * t143 + t175 * t99) * t125 + (t142 * t100 + t143 * t102 + t175 * t98) * t124 + (t119 * t175 + t120 * t142 + t121 * t143) * t156) / 0.2e1 + t156 * ((t101 * t169 + t103 * t170 + t260 * t99) * t125 + (t100 * t169 + t102 * t170 + t260 * t98) * t124 + (t119 * t260 + t169 * t120 + t170 * t121) * t156) / 0.2e1 + (t234 * t221 + t229 * t223 + t270 * V_base(6) + (-t197 * t221 + t200 * t223 + Icges(1,4)) * V_base(5) + (-t198 * t221 + t201 * t223 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t229 * t221 - t234 * t223 + t269 * V_base(6) + (t197 * t223 + t200 * t221 + Icges(1,2)) * V_base(5) + (t198 * t223 + t201 * t221 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t196 * t222 + t199 * t220 + Icges(1,3) + Icges(2,3)) * V_base(6) + (t159 * t222 + t161 * t220 + t269) * V_base(5) + (t160 * t222 + t162 * t220 + t270) * V_base(4)) * V_base(6) / 0.2e1;
T = t1;
