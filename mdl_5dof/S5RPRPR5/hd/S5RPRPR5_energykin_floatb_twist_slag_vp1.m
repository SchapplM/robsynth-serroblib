% Calculate kinetic energy for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:49
% EndTime: 2022-01-23 09:24:52
% DurationCPUTime: 3.44s
% Computational Cost: add. (1490->350), mult. (1916->471), div. (0->0), fcn. (1866->10), ass. (0->162)
t285 = -Icges(5,3) - Icges(4,3);
t228 = qJ(3) + pkin(9);
t217 = cos(t228);
t230 = cos(pkin(8));
t235 = cos(qJ(1));
t216 = sin(t228);
t233 = sin(qJ(1));
t265 = t233 * t216;
t165 = -t217 * t235 - t230 * t265;
t264 = t233 * t217;
t166 = -t216 * t235 + t230 * t264;
t234 = cos(qJ(3));
t261 = t234 * t235;
t232 = sin(qJ(3));
t263 = t233 * t232;
t177 = -t230 * t263 - t261;
t262 = t233 * t234;
t268 = t232 * t235;
t178 = t230 * t262 - t268;
t229 = sin(pkin(8));
t271 = t229 * t233;
t284 = Icges(4,5) * t178 + Icges(5,5) * t166 + Icges(4,6) * t177 + Icges(5,6) * t165 - t285 * t271;
t269 = t230 * t235;
t167 = -t216 * t269 + t264;
t168 = t217 * t269 + t265;
t179 = -t230 * t268 + t262;
t180 = t230 * t261 + t263;
t270 = t229 * t235;
t283 = Icges(4,5) * t180 + Icges(5,5) * t168 + Icges(4,6) * t179 + Icges(5,6) * t167 - t285 * t270;
t282 = t285 * t230 + (Icges(4,5) * t234 + Icges(5,5) * t217 - Icges(4,6) * t232 - Icges(5,6) * t216) * t229;
t277 = pkin(3) * t232;
t276 = pkin(3) * t234;
t275 = pkin(2) * t230 + pkin(1);
t231 = qJ(4) + pkin(6);
t274 = Icges(2,4) * t233;
t273 = Icges(3,4) * t229;
t272 = Icges(3,4) * t230;
t218 = qJ(5) + t228;
t213 = sin(t218);
t267 = t233 * t213;
t214 = cos(t218);
t266 = t233 * t214;
t169 = t229 * t231 + t230 * t276 + t275;
t192 = pkin(6) * t229 + t275;
t260 = t169 - t192;
t215 = qJ(2) + t277;
t259 = pkin(4) * t216 - t215 + t277;
t221 = t233 * qJ(2);
t225 = t235 * pkin(1);
t207 = t225 + t221;
t258 = qJD(3) * t229;
t257 = qJD(4) * t229;
t256 = qJD(5) * t229;
t224 = t233 * pkin(1);
t205 = -qJ(2) * t235 + t224;
t255 = V_base(4) * t205 + V_base(3);
t254 = V_base(5) * pkin(5) + V_base(1);
t251 = -pkin(2) - t276;
t191 = t235 * t258 + V_base(4);
t190 = t233 * t258 + V_base(5);
t219 = V_base(6) + qJD(1);
t250 = qJD(2) * t233 + t254;
t249 = rSges(3,1) * t230 - rSges(3,2) * t229;
t248 = Icges(3,1) * t230 - t273;
t247 = -Icges(3,2) * t229 + t272;
t246 = Icges(3,5) * t230 - Icges(3,6) * t229;
t245 = -qJD(2) * t235 + t219 * t207 + V_base(2);
t185 = pkin(4) * t217 - t251;
t227 = -pkin(7) - t231;
t244 = t185 * t230 - t227 * t229 - t169;
t173 = t192 * t233 - t224;
t198 = pkin(2) * t229 - pkin(6) * t230;
t243 = V_base(5) * t198 + (-t173 - t205) * t219 + t250;
t174 = t192 * t235 - t225;
t242 = V_base(4) * t173 + (-t174 - t207) * V_base(5) + t255;
t241 = (-Icges(3,3) * t235 + t233 * t246) * V_base(5) + (Icges(3,3) * t233 + t235 * t246) * V_base(4) + (Icges(3,5) * t229 + Icges(3,6) * t230) * t219;
t172 = t229 * t276 + (pkin(6) - t231) * t230;
t240 = t190 * t172 + t235 * t257 + t243;
t239 = t219 * t174 + (-pkin(5) - t198) * V_base(4) + t245;
t124 = (qJ(2) - t215) * t235 + t260 * t233;
t238 = -qJD(4) * t230 + t191 * t124 + t242;
t125 = t215 * t233 + t235 * t260 - t221;
t197 = -qJD(3) * t230 + t219;
t237 = t197 * t125 + t233 * t257 + t239;
t160 = -Icges(3,6) * t235 + t233 * t247;
t161 = Icges(3,6) * t233 + t235 * t247;
t162 = -Icges(3,5) * t235 + t233 * t248;
t163 = Icges(3,5) * t233 + t235 * t248;
t194 = Icges(3,2) * t230 + t273;
t195 = Icges(3,1) * t229 + t272;
t236 = (-t161 * t229 + t163 * t230) * V_base(4) + (-t160 * t229 + t162 * t230) * V_base(5) + (-t194 * t229 + t195 * t230) * t219;
t223 = Icges(2,4) * t235;
t208 = rSges(2,1) * t235 - t233 * rSges(2,2);
t206 = t233 * rSges(2,1) + rSges(2,2) * t235;
t204 = Icges(2,1) * t235 - t274;
t203 = Icges(2,1) * t233 + t223;
t202 = -Icges(2,2) * t233 + t223;
t201 = Icges(2,2) * t235 + t274;
t200 = Icges(2,5) * t235 - Icges(2,6) * t233;
t199 = Icges(2,5) * t233 + Icges(2,6) * t235;
t196 = rSges(3,1) * t229 + rSges(3,2) * t230;
t188 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t187 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t186 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t181 = (-qJD(3) - qJD(5)) * t230 + t219;
t176 = t235 * t256 + t191;
t175 = t233 * t256 + t190;
t171 = t233 * rSges(3,3) + t235 * t249;
t170 = -rSges(3,3) * t235 + t233 * t249;
t164 = -rSges(4,3) * t230 + (rSges(4,1) * t234 - rSges(4,2) * t232) * t229;
t157 = -Icges(4,5) * t230 + (Icges(4,1) * t234 - Icges(4,4) * t232) * t229;
t156 = -Icges(4,6) * t230 + (Icges(4,4) * t234 - Icges(4,2) * t232) * t229;
t153 = t214 * t269 + t267;
t152 = -t213 * t269 + t266;
t151 = -t213 * t235 + t230 * t266;
t150 = -t214 * t235 - t230 * t267;
t148 = -rSges(5,3) * t230 + (rSges(5,1) * t217 - rSges(5,2) * t216) * t229;
t147 = -Icges(5,5) * t230 + (Icges(5,1) * t217 - Icges(5,4) * t216) * t229;
t146 = -Icges(5,6) * t230 + (Icges(5,4) * t217 - Icges(5,2) * t216) * t229;
t144 = V_base(5) * rSges(2,3) - t206 * t219 + t254;
t143 = t208 * t219 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t142 = -rSges(6,3) * t230 + (rSges(6,1) * t214 - rSges(6,2) * t213) * t229;
t141 = -Icges(6,5) * t230 + (Icges(6,1) * t214 - Icges(6,4) * t213) * t229;
t140 = -Icges(6,6) * t230 + (Icges(6,4) * t214 - Icges(6,2) * t213) * t229;
t139 = -Icges(6,3) * t230 + (Icges(6,5) * t214 - Icges(6,6) * t213) * t229;
t138 = t206 * V_base(4) - t208 * V_base(5) + V_base(3);
t136 = (t227 + t231) * t230 + (t185 + t251) * t229;
t135 = t180 * rSges(4,1) + t179 * rSges(4,2) + rSges(4,3) * t270;
t134 = rSges(4,1) * t178 + rSges(4,2) * t177 + rSges(4,3) * t271;
t133 = Icges(4,1) * t180 + Icges(4,4) * t179 + Icges(4,5) * t270;
t132 = Icges(4,1) * t178 + Icges(4,4) * t177 + Icges(4,5) * t271;
t131 = Icges(4,4) * t180 + Icges(4,2) * t179 + Icges(4,6) * t270;
t130 = Icges(4,4) * t178 + Icges(4,2) * t177 + Icges(4,6) * t271;
t127 = t168 * rSges(5,1) + t167 * rSges(5,2) + rSges(5,3) * t270;
t126 = rSges(5,1) * t166 + rSges(5,2) * t165 + rSges(5,3) * t271;
t123 = Icges(5,1) * t168 + Icges(5,4) * t167 + Icges(5,5) * t270;
t122 = Icges(5,1) * t166 + Icges(5,4) * t165 + Icges(5,5) * t271;
t121 = Icges(5,4) * t168 + Icges(5,2) * t167 + Icges(5,6) * t270;
t120 = Icges(5,4) * t166 + Icges(5,2) * t165 + Icges(5,6) * t271;
t117 = t153 * rSges(6,1) + t152 * rSges(6,2) + rSges(6,3) * t270;
t116 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t271;
t115 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t270;
t114 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t271;
t113 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t270;
t112 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t271;
t111 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t270;
t110 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t271;
t107 = t196 * V_base(5) + (-t170 - t205) * t219 + t250;
t106 = t219 * t171 + (-pkin(5) - t196) * V_base(4) + t245;
t105 = t170 * V_base(4) + (-t171 - t207) * V_base(5) + t255;
t104 = t233 * t259 + t235 * t244 + t207;
t103 = t224 + (-qJ(2) - t259) * t235 + t244 * t233;
t102 = -t134 * t197 + t164 * t190 + t243;
t101 = t197 * t135 - t191 * t164 + t239;
t100 = t134 * t191 - t135 * t190 + t242;
t99 = t148 * t190 + (-t124 - t126) * t197 + t240;
t98 = t197 * t127 + (-t148 - t172) * t191 + t237;
t97 = t126 * t191 + (-t125 - t127) * t190 + t238;
t96 = -t116 * t181 + t136 * t190 + t142 * t175 + (-t103 - t124) * t197 + t240;
t95 = t197 * t104 + t181 * t117 - t176 * t142 + (-t136 - t172) * t191 + t237;
t94 = t103 * t191 + t116 * t176 - t117 * t175 + (-t104 - t125) * t190 + t238;
t1 = m(1) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(2) * (t138 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(5) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t176 * ((t111 * t270 + t152 * t113 + t153 * t115) * t176 + (t110 * t270 + t152 * t112 + t153 * t114) * t175 + (t139 * t270 + t152 * t140 + t153 * t141) * t181) / 0.2e1 + t175 * ((t111 * t271 + t113 * t150 + t115 * t151) * t176 + (t110 * t271 + t112 * t150 + t114 * t151) * t175 + (t139 * t271 + t140 * t150 + t141 * t151) * t181) / 0.2e1 + t181 * ((-t110 * t175 - t111 * t176 - t139 * t181) * t230 + ((-t113 * t213 + t115 * t214) * t176 + (-t112 * t213 + t114 * t214) * t175 + (-t140 * t213 + t141 * t214) * t181) * t229) / 0.2e1 + ((t146 * t165 + t147 * t166 + t156 * t177 + t157 * t178 + t282 * t271) * t197 + (t121 * t165 + t123 * t166 + t131 * t177 + t133 * t178 + t283 * t271) * t191 + (t120 * t165 + t122 * t166 + t130 * t177 + t132 * t178 + t284 * t271) * t190) * t190 / 0.2e1 + ((t167 * t146 + t168 * t147 + t179 * t156 + t180 * t157 + t282 * t270) * t197 + (t167 * t121 + t168 * t123 + t179 * t131 + t180 * t133 + t283 * t270) * t191 + (t167 * t120 + t168 * t122 + t179 * t130 + t180 * t132 + t284 * t270) * t190) * t191 / 0.2e1 + ((-t284 * t190 - t283 * t191 - t282 * t197) * t230 + ((-t146 * t216 + t147 * t217 - t156 * t232 + t157 * t234) * t197 + (-t121 * t216 + t123 * t217 - t131 * t232 + t133 * t234) * t191 + (-t120 * t216 + t122 * t217 - t130 * t232 + t132 * t234) * t190) * t229) * t197 / 0.2e1 + ((t160 * t230 + t162 * t229 + t199) * V_base(5) + (t161 * t230 + t163 * t229 + t200) * V_base(4) + (t194 * t230 + t195 * t229 + Icges(2,3)) * t219) * t219 / 0.2e1 + (t200 * t219 + t233 * t241 + t235 * t236 + (-t233 * t201 + t203 * t235 + Icges(1,4)) * V_base(5) + (-t233 * t202 + t204 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t199 * t219 + t233 * t236 - t235 * t241 + (t201 * t235 + t233 * t203 + Icges(1,2)) * V_base(5) + (t202 * t235 + t233 * t204 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
