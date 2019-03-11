% Calculate kinetic energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:27
% EndTime: 2019-03-09 01:31:29
% DurationCPUTime: 2.29s
% Computational Cost: add. (1466->298), mult. (1250->396), div. (0->0), fcn. (1030->10), ass. (0->156)
t285 = Icges(3,4) + Icges(4,6);
t284 = Icges(3,1) + Icges(4,2);
t283 = -Icges(4,4) + Icges(3,5);
t282 = Icges(4,5) - Icges(3,6);
t281 = Icges(3,2) + Icges(4,3);
t202 = qJ(1) + pkin(9);
t195 = cos(t202);
t280 = t285 * t195;
t193 = sin(t202);
t279 = t285 * t193;
t278 = t284 * t193 + t280;
t277 = t284 * t195 - t279;
t204 = cos(pkin(10));
t203 = sin(pkin(10));
t263 = Icges(5,4) * t203;
t231 = Icges(5,2) * t204 + t263;
t117 = Icges(5,6) * t193 - t195 * t231;
t262 = Icges(5,4) * t204;
t233 = Icges(5,1) * t203 + t262;
t119 = Icges(5,5) * t193 - t195 * t233;
t276 = t117 * t204 + t119 * t203 - t281 * t195 - t279;
t116 = Icges(5,6) * t195 + t193 * t231;
t118 = Icges(5,5) * t195 + t193 * t233;
t275 = t116 * t204 + t118 * t203 + t281 * t193 - t280;
t207 = sin(qJ(1));
t209 = cos(qJ(1));
t274 = Icges(2,5) * t207 + Icges(2,6) * t209 + t283 * t193 - t282 * t195;
t273 = Icges(2,5) * t209 - Icges(2,6) * t207 + t282 * t193 + t283 * t195;
t201 = pkin(10) + qJ(5);
t194 = cos(t201);
t192 = sin(t201);
t261 = Icges(6,4) * t192;
t230 = Icges(6,2) * t194 + t261;
t105 = Icges(6,6) * t195 + t193 * t230;
t106 = Icges(6,6) * t193 - t195 * t230;
t260 = Icges(6,4) * t194;
t232 = Icges(6,1) * t192 + t260;
t107 = Icges(6,5) * t195 + t193 * t232;
t108 = Icges(6,5) * t193 - t195 * t232;
t147 = -Icges(6,2) * t192 + t260;
t152 = Icges(6,1) * t194 - t261;
t171 = qJD(5) * t193 + V_base(5);
t172 = qJD(5) * t195 + V_base(4);
t196 = V_base(6) + qJD(1);
t272 = (t105 * t194 + t107 * t192) * t172 + (t106 * t194 + t108 * t192) * t171 + (t147 * t194 + t152 * t192) * t196;
t270 = pkin(1) * t207;
t269 = pkin(1) * t209;
t268 = pkin(4) * t203;
t267 = pkin(4) * t204;
t266 = -pkin(6) - qJ(2);
t265 = Icges(2,4) * t207;
t257 = qJ(4) * t193;
t256 = qJ(4) * t195;
t255 = t193 * t194;
t206 = sin(qJ(6));
t254 = t193 * t206;
t208 = cos(qJ(6));
t253 = t193 * t208;
t252 = t194 * t195;
t251 = t195 * t206;
t250 = t195 * t208;
t248 = qJD(6) * t194;
t247 = -pkin(3) + t266;
t246 = t196 * t269 + V_base(2);
t245 = V_base(5) * pkin(6) + V_base(1);
t155 = pkin(2) * t193 - qJ(3) * t195;
t242 = -t155 - t270;
t159 = pkin(2) * t195 + qJ(3) * t193;
t241 = -t159 - t269;
t240 = V_base(5) * qJ(2) + t245;
t239 = V_base(4) * t270 + qJD(2) + V_base(3);
t238 = qJD(3) * t193 + t240;
t237 = pkin(5) * t192 - pkin(8) * t194;
t236 = V_base(4) * t155 + t239;
t235 = rSges(5,1) * t203 + rSges(5,2) * t204;
t234 = rSges(6,1) * t192 + rSges(6,2) * t194;
t229 = Icges(5,5) * t203 + Icges(5,6) * t204;
t228 = Icges(6,5) * t192 + Icges(6,6) * t194;
t169 = -Icges(5,2) * t203 + t262;
t170 = Icges(5,1) * t204 - t263;
t222 = t169 * t204 + t170 * t203;
t221 = V_base(4) * t257 + t236;
t220 = t242 - t257;
t219 = t241 - t256;
t218 = -qJD(3) * t195 + t196 * t159 + t246;
t217 = V_base(5) * pkin(3) + qJD(4) * t195 + t238;
t128 = pkin(7) * t193 - t195 * t268;
t216 = -t128 + t220;
t215 = V_base(5) * t267 + t217;
t214 = (Icges(6,3) * t195 + t193 * t228) * t172 + (Icges(6,3) * t193 - t195 * t228) * t171 + (Icges(6,5) * t194 - Icges(6,6) * t192) * t196;
t213 = qJD(4) * t193 + t196 * t256 + t218;
t212 = (Icges(5,3) * t195 + t193 * t229) * V_base(4) + (Icges(5,3) * t193 - t195 * t229) * V_base(5) + (Icges(5,5) * t204 - Icges(5,6) * t203) * t196;
t129 = pkin(7) * t195 + t193 * t268;
t211 = V_base(4) * t128 + (-t129 + t219) * V_base(5) + t221;
t210 = t196 * t129 + (t247 - t267) * V_base(4) + t213;
t198 = Icges(2,4) * t209;
t181 = rSges(2,1) * t209 - t207 * rSges(2,2);
t180 = t207 * rSges(2,1) + rSges(2,2) * t209;
t179 = Icges(2,1) * t209 - t265;
t178 = Icges(2,1) * t207 + t198;
t177 = -Icges(2,2) * t207 + t198;
t176 = Icges(2,2) * t209 + t265;
t173 = rSges(5,1) * t204 - rSges(5,2) * t203;
t167 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t166 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t165 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t163 = qJD(6) * t192 + t196;
t162 = pkin(5) * t194 + pkin(8) * t192;
t161 = rSges(3,1) * t195 - rSges(3,2) * t193;
t160 = -rSges(4,2) * t195 + rSges(4,3) * t193;
t158 = rSges(6,1) * t194 - rSges(6,2) * t192;
t157 = rSges(3,1) * t193 + rSges(3,2) * t195;
t156 = -rSges(4,2) * t193 - rSges(4,3) * t195;
t137 = -t192 * t250 + t254;
t136 = t192 * t251 + t253;
t135 = t192 * t253 + t251;
t134 = -t192 * t254 + t250;
t133 = -t193 * t248 + t172;
t132 = t195 * t248 + t171;
t131 = t237 * t195;
t130 = t237 * t193;
t127 = rSges(7,3) * t192 + (rSges(7,1) * t208 - rSges(7,2) * t206) * t194;
t126 = V_base(5) * rSges(2,3) - t180 * t196 + t245;
t125 = t181 * t196 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t124 = Icges(7,5) * t192 + (Icges(7,1) * t208 - Icges(7,4) * t206) * t194;
t123 = Icges(7,6) * t192 + (Icges(7,4) * t208 - Icges(7,2) * t206) * t194;
t122 = Icges(7,3) * t192 + (Icges(7,5) * t208 - Icges(7,6) * t206) * t194;
t121 = rSges(5,3) * t193 - t195 * t235;
t120 = rSges(5,3) * t195 + t193 * t235;
t113 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t111 = rSges(6,3) * t193 - t195 * t234;
t110 = rSges(6,3) * t195 + t193 * t234;
t102 = V_base(5) * rSges(3,3) + (-t157 - t270) * t196 + t240;
t101 = t161 * t196 + (-rSges(3,3) + t266) * V_base(4) + t246;
t100 = V_base(4) * t157 + (-t161 - t269) * V_base(5) + t239;
t99 = rSges(7,1) * t137 + rSges(7,2) * t136 + rSges(7,3) * t252;
t98 = rSges(7,1) * t135 + rSges(7,2) * t134 - rSges(7,3) * t255;
t97 = Icges(7,1) * t137 + Icges(7,4) * t136 + Icges(7,5) * t252;
t96 = Icges(7,1) * t135 + Icges(7,4) * t134 - Icges(7,5) * t255;
t95 = Icges(7,4) * t137 + Icges(7,2) * t136 + Icges(7,6) * t252;
t94 = Icges(7,4) * t135 + Icges(7,2) * t134 - Icges(7,6) * t255;
t93 = Icges(7,5) * t137 + Icges(7,6) * t136 + Icges(7,3) * t252;
t92 = Icges(7,5) * t135 + Icges(7,6) * t134 - Icges(7,3) * t255;
t91 = V_base(5) * rSges(4,1) + (-t156 + t242) * t196 + t238;
t90 = t160 * t196 + (-rSges(4,1) + t266) * V_base(4) + t218;
t89 = V_base(4) * t156 + (-t160 + t241) * V_base(5) + t236;
t88 = t173 * V_base(5) + (-t121 + t220) * t196 + t217;
t87 = t120 * t196 + (-t173 + t247) * V_base(4) + t213;
t86 = V_base(4) * t121 + (-t120 + t219) * V_base(5) + t221;
t85 = t158 * t171 + (-t111 + t216) * t196 + t215;
t84 = t110 * t196 - t158 * t172 + t210;
t83 = -t171 * t110 + t172 * t111 + t211;
t82 = t127 * t132 + t162 * t171 - t163 * t99 + (t131 + t216) * t196 + t215;
t81 = -t127 * t133 + t130 * t196 - t162 * t172 + t163 * t98 + t210;
t80 = -t171 * t130 - t172 * t131 - t132 * t98 + t133 * t99 + t211;
t1 = t163 * ((t122 * t163 + t93 * t132 + t92 * t133) * t192 + ((-t206 * t94 + t208 * t96) * t133 + (-t206 * t95 + t208 * t97) * t132 + (-t123 * t206 + t124 * t208) * t163) * t194) / 0.2e1 + t132 * ((t136 * t94 + t137 * t96 + t92 * t252) * t133 + (t136 * t95 + t137 * t97 + t93 * t252) * t132 + (t122 * t252 + t123 * t136 + t124 * t137) * t163) / 0.2e1 + t133 * ((t134 * t94 + t135 * t96 - t92 * t255) * t133 + (t134 * t95 + t135 * t97 - t255 * t93) * t132 + (-t122 * t255 + t123 * t134 + t124 * t135) * t163) / 0.2e1 + t172 * (t272 * t193 + t214 * t195) / 0.2e1 + t171 * (t214 * t193 - t272 * t195) / 0.2e1 + m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(2) * (t113 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + ((-t105 * t192 + t107 * t194) * t172 + (-t106 * t192 + t108 * t194) * t171 + (-t117 * t203 + t119 * t204 + t274) * V_base(5) + (-t116 * t203 + t118 * t204 + t273) * V_base(4) + (-t192 * t147 + t194 * t152 - t203 * t169 + t204 * t170 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t196) * t196 / 0.2e1 + (t212 * t195 + (t222 * t193 + t273) * t196 + (-t207 * t176 + t178 * t209 + t276 * t193 + t195 * t278 + Icges(1,4)) * V_base(5) + (-t207 * t177 + t209 * t179 + t275 * t193 + t277 * t195 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t212 * t193 + (-t222 * t195 + t274) * t196 + (t209 * t176 + t207 * t178 + t278 * t193 - t276 * t195 + Icges(1,2)) * V_base(5) + (t177 * t209 + t207 * t179 + t193 * t277 - t195 * t275 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
