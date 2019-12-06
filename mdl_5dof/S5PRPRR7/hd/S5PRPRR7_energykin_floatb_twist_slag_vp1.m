% Calculate kinetic energy for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:39
% EndTime: 2019-12-05 15:59:42
% DurationCPUTime: 3.57s
% Computational Cost: add. (957->284), mult. (1625->418), div. (0->0), fcn. (1519->8), ass. (0->137)
t275 = Icges(3,4) + Icges(4,6);
t274 = Icges(3,1) + Icges(4,2);
t273 = -Icges(3,2) - Icges(4,3);
t206 = cos(qJ(2));
t272 = t275 * t206;
t204 = sin(qJ(2));
t271 = t275 * t204;
t270 = Icges(4,4) - Icges(3,5);
t269 = Icges(4,5) - Icges(3,6);
t268 = t273 * t204 + t272;
t267 = t274 * t206 - t271;
t266 = Icges(4,1) + Icges(3,3);
t201 = sin(pkin(8));
t202 = cos(pkin(8));
t265 = t268 * t201 + t269 * t202;
t264 = -t269 * t201 + t268 * t202;
t263 = t267 * t201 + t270 * t202;
t262 = -t270 * t201 + t267 * t202;
t261 = t273 * t206 - t271;
t260 = t274 * t204 + t272;
t259 = t269 * t204 - t270 * t206;
t186 = -qJD(2) * t202 + V_base(5);
t187 = qJD(2) * t201 + V_base(4);
t256 = (t261 * t204 + t260 * t206) * V_base(6) + (-t264 * t204 + t262 * t206) * t187 + (-t265 * t204 + t263 * t206) * t186;
t255 = (-t270 * t204 - t269 * t206) * V_base(6) + (t266 * t201 + t259 * t202) * t187 + (t259 * t201 - t266 * t202) * t186;
t252 = t204 * pkin(6);
t205 = cos(qJ(4));
t251 = pkin(4) * t205;
t249 = Icges(2,4) * t201;
t244 = t201 * t204;
t243 = t201 * t206;
t242 = t202 * t204;
t241 = t202 * t206;
t203 = sin(qJ(4));
t240 = t203 * t204;
t239 = t204 * t205;
t226 = pkin(2) * t206 + qJ(3) * t204;
t156 = t226 * t201;
t175 = pkin(1) * t201 - pkin(5) * t202;
t238 = -t156 - t175;
t237 = qJD(3) * t204;
t236 = qJD(4) * t206;
t235 = qJD(5) * t206;
t234 = V_base(5) * qJ(1) + V_base(1);
t230 = qJD(1) + V_base(3);
t188 = qJD(4) * t204 + V_base(6);
t150 = t202 * t236 + t187;
t183 = pkin(2) * t204 - qJ(3) * t206;
t229 = t186 * t183 + t202 * t237 + t234;
t228 = rSges(3,1) * t206 - rSges(3,2) * t204;
t227 = -rSges(4,2) * t206 + rSges(4,3) * t204;
t149 = t201 * t236 + t186;
t176 = pkin(1) * t202 + pkin(5) * t201;
t219 = -V_base(4) * qJ(1) + V_base(6) * t176 + V_base(2);
t218 = pkin(4) * t240 + pkin(7) * t206;
t217 = V_base(4) * t175 - t176 * V_base(5) + t230;
t157 = t226 * t202;
t216 = V_base(6) * t157 + t201 * t237 + t219;
t213 = -qJD(3) * t206 + t187 * t156 + t217;
t163 = -t202 * pkin(3) + pkin(6) * t243;
t212 = t186 * t252 + (-t163 + t238) * V_base(6) + t229;
t162 = t201 * pkin(3) + pkin(6) * t241;
t211 = V_base(6) * t162 + (-t183 - t252) * t187 + t216;
t210 = t187 * t163 + (-t157 - t162) * t186 + t213;
t200 = qJ(4) + qJ(5);
t198 = cos(t200);
t197 = sin(t200);
t195 = Icges(2,4) * t202;
t185 = rSges(3,1) * t204 + rSges(3,2) * t206;
t184 = -rSges(4,2) * t204 - rSges(4,3) * t206;
t174 = rSges(2,1) * t202 - rSges(2,2) * t201;
t173 = rSges(2,1) * t201 + rSges(2,2) * t202;
t172 = Icges(2,1) * t202 - t249;
t171 = Icges(2,1) * t201 + t195;
t170 = -Icges(2,2) * t201 + t195;
t169 = Icges(2,2) * t202 + t249;
t166 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t165 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t164 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t160 = qJD(5) * t204 + t188;
t154 = t201 * t240 - t202 * t205;
t153 = t201 * t239 + t202 * t203;
t152 = t201 * t205 + t202 * t240;
t151 = -t201 * t203 + t202 * t239;
t148 = -t206 * t203 * pkin(4) + pkin(7) * t204;
t145 = rSges(5,3) * t204 + (-rSges(5,1) * t203 - rSges(5,2) * t205) * t206;
t144 = Icges(5,5) * t204 + (-Icges(5,1) * t203 - Icges(5,4) * t205) * t206;
t143 = Icges(5,6) * t204 + (-Icges(5,4) * t203 - Icges(5,2) * t205) * t206;
t142 = Icges(5,3) * t204 + (-Icges(5,5) * t203 - Icges(5,6) * t205) * t206;
t141 = t197 * t244 - t198 * t202;
t140 = t197 * t202 + t198 * t244;
t139 = t197 * t242 + t198 * t201;
t138 = -t197 * t201 + t198 * t242;
t137 = -rSges(4,1) * t202 + t201 * t227;
t136 = rSges(4,1) * t201 + t202 * t227;
t135 = rSges(3,3) * t201 + t202 * t228;
t134 = -rSges(3,3) * t202 + t201 * t228;
t120 = rSges(6,3) * t204 + (-rSges(6,1) * t197 - rSges(6,2) * t198) * t206;
t119 = Icges(6,5) * t204 + (-Icges(6,1) * t197 - Icges(6,4) * t198) * t206;
t118 = Icges(6,6) * t204 + (-Icges(6,4) * t197 - Icges(6,2) * t198) * t206;
t117 = Icges(6,3) * t204 + (-Icges(6,5) * t197 - Icges(6,6) * t198) * t206;
t116 = t202 * t235 + t150;
t115 = t201 * t235 + t149;
t113 = V_base(5) * rSges(2,3) - t173 * V_base(6) + t234;
t112 = t174 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t111 = t173 * V_base(4) - t174 * V_base(5) + t230;
t110 = t201 * t218 - t202 * t251;
t109 = t201 * t251 + t202 * t218;
t108 = rSges(5,1) * t154 + rSges(5,2) * t153 + rSges(5,3) * t243;
t107 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t241;
t106 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t243;
t105 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t241;
t104 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t243;
t103 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t241;
t102 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t243;
t101 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t241;
t100 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t243;
t99 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t241;
t98 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t243;
t97 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t241;
t96 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t243;
t95 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t241;
t94 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t243;
t93 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t241;
t92 = t185 * t186 + (-t134 - t175) * V_base(6) + t234;
t91 = t135 * V_base(6) - t185 * t187 + t219;
t90 = t134 * t187 - t135 * t186 + t217;
t89 = t184 * t186 + (-t137 + t238) * V_base(6) + t229;
t88 = t136 * V_base(6) + (-t183 - t184) * t187 + t216;
t87 = t137 * t187 + (-t136 - t157) * t186 + t213;
t86 = -t108 * t188 + t145 * t149 + t212;
t85 = t107 * t188 - t145 * t150 + t211;
t84 = -t107 * t149 + t108 * t150 + t210;
t83 = -t100 * t160 - t110 * t188 + t115 * t120 + t148 * t149 + t212;
t82 = t109 * t188 - t116 * t120 - t148 * t150 + t160 * t99 + t211;
t81 = t100 * t116 - t109 * t149 + t110 * t150 - t115 * t99 + t210;
t1 = m(1) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(2) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(3) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t150 * ((t101 * t241 + t103 * t151 + t105 * t152) * t150 + (t102 * t241 + t104 * t151 + t106 * t152) * t149 + (t142 * t241 + t143 * t151 + t144 * t152) * t188) / 0.2e1 + t149 * ((t101 * t243 + t103 * t153 + t105 * t154) * t150 + (t102 * t243 + t104 * t153 + t106 * t154) * t149 + (t142 * t243 + t143 * t153 + t144 * t154) * t188) / 0.2e1 + t188 * ((t101 * t150 + t102 * t149 + t142 * t188) * t204 + ((-t103 * t205 - t105 * t203) * t150 + (-t104 * t205 - t106 * t203) * t149 + (-t143 * t205 - t144 * t203) * t188) * t206) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t116 * ((t138 * t95 + t139 * t97 + t93 * t241) * t116 + (t138 * t96 + t139 * t98 + t241 * t94) * t115 + (t117 * t241 + t118 * t138 + t119 * t139) * t160) / 0.2e1 + t115 * ((t140 * t95 + t141 * t97 + t243 * t93) * t116 + (t140 * t96 + t141 * t98 + t94 * t243) * t115 + (t117 * t243 + t118 * t140 + t119 * t141) * t160) / 0.2e1 + t160 * ((t94 * t115 + t93 * t116 + t117 * t160) * t204 + ((-t197 * t97 - t198 * t95) * t116 + (-t197 * t98 - t198 * t96) * t115 + (-t118 * t198 - t119 * t197) * t160) * t206) / 0.2e1 + (t256 * t201 - t255 * t202) * t186 / 0.2e1 + (t255 * t201 + t256 * t202) * t187 / 0.2e1 + ((-t169 * t201 + t171 * t202 + Icges(1,4)) * V_base(5) + (-t170 * t201 + t172 * t202 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t169 * t202 + t171 * t201 + Icges(1,2)) * V_base(5) + (t170 * t202 + t172 * t201 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t262 * t204 + t264 * t206) * t187 + (t263 * t204 + t265 * t206) * t186 + (t260 * t204 - t261 * t206 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t202 - Icges(2,6) * t201 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t201 + Icges(2,6) * t202 + Icges(1,6));
T = t1;
