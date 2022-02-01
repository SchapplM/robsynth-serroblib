% Calculate time derivative of joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:41
% EndTime: 2022-01-20 12:01:51
% DurationCPUTime: 3.77s
% Computational Cost: add. (12576->403), mult. (7564->563), div. (0->0), fcn. (5570->10), ass. (0->233)
t183 = qJ(1) + qJ(2);
t178 = qJ(3) + t183;
t170 = sin(t178);
t171 = cos(t178);
t186 = cos(qJ(4));
t249 = qJD(4) * t186;
t181 = qJD(1) + qJD(2);
t173 = qJD(3) + t181;
t184 = sin(qJ(4));
t263 = t173 * t184;
t298 = -t170 * t249 - t171 * t263;
t277 = rSges(5,2) * t184;
t279 = rSges(5,1) * t186;
t297 = -t277 + t279;
t273 = Icges(5,4) * t186;
t211 = -Icges(5,2) * t184 + t273;
t202 = t211 * t171;
t95 = Icges(5,6) * t170 + t202;
t274 = Icges(5,4) * t184;
t213 = Icges(5,1) * t186 - t274;
t204 = t213 * t171;
t97 = Icges(5,5) * t170 + t204;
t215 = t184 * t95 - t186 * t97;
t296 = t215 * t170;
t94 = -Icges(5,6) * t171 + t211 * t170;
t96 = -Icges(5,5) * t171 + t213 * t170;
t217 = t184 * t94 - t186 * t96;
t295 = t217 * t171;
t182 = qJ(4) + qJ(5);
t174 = sin(t182);
t176 = cos(t182);
t271 = Icges(6,4) * t176;
t210 = -Icges(6,2) * t174 + t271;
t201 = t210 * t171;
t85 = Icges(6,6) * t170 + t201;
t272 = Icges(6,4) * t174;
t212 = Icges(6,1) * t176 - t272;
t203 = t212 * t171;
t87 = Icges(6,5) * t170 + t203;
t221 = t174 * t85 - t176 * t87;
t294 = t221 * t170;
t84 = -Icges(6,6) * t171 + t210 * t170;
t86 = -Icges(6,5) * t171 + t212 * t170;
t222 = t174 * t84 - t176 * t86;
t293 = t222 * t171;
t158 = t170 * rSges(6,3);
t278 = rSges(6,1) * t176;
t292 = t171 * t278 + t158;
t180 = qJD(4) + qJD(5);
t132 = Icges(6,2) * t176 + t272;
t133 = Icges(6,1) * t174 + t271;
t207 = t132 * t174 - t133 * t176;
t208 = Icges(6,5) * t176 - Icges(6,6) * t174;
t291 = t207 * t173 + t208 * t180;
t154 = Icges(5,2) * t186 + t274;
t155 = Icges(5,1) * t184 + t273;
t206 = t154 * t184 - t155 * t186;
t209 = Icges(5,5) * t186 - Icges(5,6) * t184;
t290 = t209 * qJD(4) + t206 * t173;
t289 = 2 * m(3);
t288 = 2 * m(4);
t287 = 2 * m(5);
t286 = 2 * m(6);
t285 = t170 / 0.2e1;
t284 = -t171 / 0.2e1;
t143 = t297 * qJD(4);
t283 = m(5) * t143;
t157 = rSges(5,1) * t184 + rSges(5,2) * t186;
t282 = m(5) * t157;
t175 = sin(t183);
t281 = pkin(2) * t175;
t163 = t170 * pkin(8);
t185 = sin(qJ(1));
t280 = t185 * pkin(1);
t276 = rSges(6,2) * t174;
t146 = t170 * t276;
t254 = t171 * rSges(6,3) + t146;
t267 = t170 * t176;
t88 = rSges(6,1) * t267 - t254;
t242 = t171 * t276;
t89 = -t242 + t292;
t50 = t170 * t88 + t171 * t89;
t275 = pkin(1) * qJD(1);
t159 = t170 * rSges(5,3);
t270 = t132 * t180;
t269 = t133 * t180;
t268 = t170 * t173;
t266 = t170 * t186;
t265 = t171 * t173;
t188 = -pkin(9) - pkin(8);
t264 = t171 * t188;
t262 = t173 * t188;
t261 = t174 * t180;
t260 = t175 * t181;
t259 = t176 * t180;
t177 = cos(t183);
t258 = t177 * t181;
t257 = rSges(6,3) * t265 + t173 * t146;
t241 = t170 * t263;
t256 = rSges(5,2) * t241 + rSges(5,3) * t265;
t250 = qJD(4) * t184;
t237 = t170 * t250;
t255 = -pkin(4) * t237 - t170 * t262;
t253 = t171 * rSges(5,3) + t170 * t277;
t252 = -t171 * pkin(3) - t163;
t251 = t170 ^ 2 + t171 ^ 2;
t245 = rSges(6,2) * t259;
t239 = -t173 * t242 + (-t261 * rSges(6,1) - t245) * t170;
t248 = t170 * (t173 * t292 + t239) + t171 * (-t171 * t245 + (-t171 * t261 - t173 * t267) * rSges(6,1) + t257) + t88 * t265;
t247 = pkin(2) * t258;
t246 = pkin(2) * t260;
t187 = cos(qJ(1));
t244 = t187 * t275;
t243 = t185 * t275;
t238 = -rSges(5,1) * t237 + t298 * rSges(5,2);
t235 = t171 * t250;
t234 = t171 * t249;
t233 = t268 / 0.2e1;
t232 = t265 / 0.2e1;
t231 = -pkin(3) - t279;
t134 = rSges(6,1) * t174 + rSges(6,2) * t176;
t230 = -pkin(4) * t184 - t134;
t197 = Icges(6,6) * t173 - t270;
t44 = t197 * t171 - t210 * t268;
t229 = t180 * t87 + t44;
t45 = t197 * t170 + t173 * t201;
t228 = t180 * t86 + t45;
t198 = Icges(6,5) * t173 - t269;
t46 = t198 * t171 - t212 * t268;
t227 = -t180 * t85 + t46;
t47 = t198 * t170 + t173 * t203;
t226 = t180 * t84 - t47;
t172 = pkin(4) * t186 + pkin(3);
t225 = -t172 - t278;
t136 = t177 * rSges(3,1) - rSges(3,2) * t175;
t124 = t171 * rSges(4,1) - rSges(4,2) * t170;
t224 = -t170 * t188 + t171 * t172;
t116 = -rSges(3,1) * t258 + rSges(3,2) * t260;
t107 = -rSges(4,1) * t265 + rSges(4,2) * t268;
t169 = pkin(2) * t177;
t110 = t124 + t169;
t82 = -Icges(6,3) * t171 + t208 * t170;
t20 = -t222 * t170 - t171 * t82;
t199 = t208 * t171;
t83 = Icges(6,3) * t170 + t199;
t21 = -t171 * t83 - t294;
t22 = t170 * t82 - t293;
t23 = t170 * t83 - t221 * t171;
t131 = Icges(6,5) * t174 + Icges(6,6) * t176;
t196 = Icges(6,3) * t173 - t131 * t180;
t42 = t196 * t171 - t208 * t268;
t43 = t196 * t170 + t173 * t199;
t223 = -t171 * ((t171 * t43 + (t21 + t293) * t173) * t171 + (t20 * t173 + (-t174 * t44 + t176 * t46 - t85 * t259 - t87 * t261) * t170 + (-t221 * t173 + t228 * t174 + t226 * t176 - t42) * t171) * t170) + t170 * ((t170 * t42 + (t22 + t294) * t173) * t170 + (t23 * t173 + (t174 * t45 - t176 * t47 + t84 * t259 + t86 * t261) * t171 + (-t222 * t173 - t229 * t174 + t227 * t176 - t43) * t170) * t171) + (t170 * t21 - t171 * t20) * t268 + (t170 * t23 - t171 * t22) * t265;
t135 = -rSges(3,1) * t175 - rSges(3,2) * t177;
t123 = -rSges(4,1) * t170 - rSges(4,2) * t171;
t218 = t184 * t96 + t186 * t94;
t216 = t184 * t97 + t186 * t95;
t214 = t225 * t170;
t153 = Icges(5,5) * t184 + Icges(5,6) * t186;
t99 = t297 * t171 + t159;
t115 = t135 * t181;
t106 = t123 * t173;
t112 = t210 * t180;
t113 = t212 * t180;
t190 = t131 * t173 + (t113 - t270) * t176 + (-t112 - t269) * t174;
t205 = (t291 * t170 + t190 * t171 + t227 * t174 + t229 * t176) * t285 + (t190 * t170 - t171 * t291 - t226 * t174 + t228 * t176) * t284 + (-t131 * t171 - t207 * t170 + t174 * t86 + t176 * t84) * t233 + (t131 * t170 - t207 * t171 + t174 * t87 + t176 * t85) * t232;
t200 = t209 * t171;
t109 = t123 - t281;
t73 = t99 - t252;
t81 = t107 - t247;
t164 = t171 * pkin(8);
t72 = t231 * t170 + t164 + t253;
t70 = t169 + t73;
t68 = t224 + t89;
t195 = Icges(5,5) * t173 - t155 * qJD(4);
t194 = Icges(5,6) * t173 - t154 * qJD(4);
t193 = Icges(5,3) * t173 - t153 * qJD(4);
t64 = t169 + t68;
t67 = t214 + t254 - t264;
t80 = t106 - t246;
t69 = t72 - t281;
t63 = t67 - t281;
t138 = t211 * qJD(4);
t139 = t213 * qJD(4);
t189 = -t138 * t184 + t139 * t186 + t153 * t173 + (-t154 * t186 - t155 * t184) * qJD(4);
t192 = t205 + (-t215 * qJD(4) + t290 * t170 + t189 * t171 + t184 * (t195 * t171 - t213 * t268) + t186 * (t194 * t171 - t211 * t268)) * t285 + (-t217 * qJD(4) + t189 * t170 - t171 * t290 + t184 * (t195 * t170 + t173 * t204) + t186 * (t194 * t170 + t173 * t202)) * t284 + (-t153 * t171 - t206 * t170 + t218) * t233 + (t153 * t170 - t206 * t171 + t216) * t232;
t191 = t176 * t112 + t174 * t113 - t132 * t261 + t133 * t259 + t186 * t138 + t184 * t139 - t154 * t250 + t155 * t249;
t35 = (t231 * t171 + (-rSges(5,3) - pkin(8)) * t170) * t173 - t238;
t25 = (t225 * t171 - t158) * t173 - t239 - t255;
t33 = t35 - t247;
t19 = t25 - t247;
t144 = pkin(8) * t265;
t34 = -rSges(5,2) * t234 - pkin(3) * t268 + t144 + (-t173 * t266 - t235) * rSges(5,1) + t256;
t32 = t34 - t246;
t24 = t173 * t214 + (-pkin(4) * t250 - t134 * t180 - t262) * t171 + t257;
t18 = t24 - t246;
t179 = t187 * pkin(1);
t118 = t136 + t179;
t117 = t135 - t280;
t114 = (-t276 + t278) * t180;
t103 = t116 - t244;
t102 = t115 - t243;
t101 = t110 + t179;
t100 = t109 - t280;
t98 = rSges(5,1) * t266 - t253;
t93 = Icges(5,3) * t170 + t200;
t92 = -Icges(5,3) * t171 + t209 * t170;
t91 = t230 * t171;
t90 = t230 * t170;
t79 = t224 + t252;
t78 = t264 + t164 + t170 * (-pkin(3) + t172);
t75 = t81 - t244;
t74 = t80 - t243;
t66 = t179 + t70;
t65 = t69 - t280;
t62 = t179 + t64;
t61 = t63 - t280;
t56 = t193 * t170 + t173 * t200;
t55 = t193 * t171 - t209 * t268;
t49 = t298 * pkin(4) - t114 * t170 - t134 * t265;
t48 = t134 * t268 - t114 * t171 + (-t234 + t241) * pkin(4);
t31 = t33 - t244;
t30 = t32 - t243;
t29 = t170 * t93 - t215 * t171;
t28 = t170 * t92 - t295;
t27 = -t171 * t93 - t296;
t26 = -t217 * t170 - t171 * t92;
t17 = t19 - t244;
t16 = t18 - t243;
t15 = t170 * t78 + t171 * t79 + t50;
t8 = -t89 * t268 + t248;
t3 = t170 * t255 + t171 * (-pkin(4) * t235 - t144) + ((t78 - t264) * t171 + (-t79 - t89 - t163) * t170) * t173 + t248;
t1 = [(t16 * t62 + t17 * t61) * t286 + (t30 * t66 + t31 * t65) * t287 + (t100 * t75 + t101 * t74) * t288 + (t102 * t118 + t103 * t117) * t289 + t191; m(6) * (t16 * t64 + t17 * t63 + t18 * t62 + t19 * t61) + m(5) * (t30 * t70 + t31 * t69 + t32 * t66 + t33 * t65) + m(4) * (t100 * t81 + t101 * t80 + t109 * t75 + t110 * t74) + m(3) * (t102 * t136 + t103 * t135 + t115 * t118 + t116 * t117) + t191; (t18 * t64 + t19 * t63) * t286 + (t32 * t70 + t33 * t69) * t287 + (t109 * t81 + t110 * t80) * t288 + (t115 * t136 + t116 * t135) * t289 + t191; m(6) * (t16 * t68 + t17 * t67 + t24 * t62 + t25 * t61) + m(5) * (t30 * t73 + t31 * t72 + t34 * t66 + t35 * t65) + m(4) * (t100 * t107 + t101 * t106 + t123 * t75 + t124 * t74) + t191; m(6) * (t18 * t68 + t19 * t67 + t24 * t64 + t25 * t63) + m(5) * (t32 * t73 + t33 * t72 + t34 * t70 + t35 * t69) + m(4) * (t106 * t110 + t107 * t109 + t123 * t81 + t124 * t80) + t191; (t24 * t68 + t25 * t67) * t286 + (t34 * t73 + t35 * t72) * t287 + (t106 * t124 + t107 * t123) * t288 + t191; m(6) * (t16 * t90 + t17 * t91 + t48 * t61 + t49 * t62) + ((-t173 * t66 - t31) * t171 + (t173 * t65 - t30) * t170) * t282 + (-t170 * t66 - t171 * t65) * t283 + t192; ((-t173 * t70 - t33) * t171 + (t173 * t69 - t32) * t170) * t282 + m(6) * (t18 * t90 + t19 * t91 + t48 * t63 + t49 * t64) + (-t170 * t70 - t171 * t69) * t283 + t192; (-t170 * t73 - t171 * t72) * t283 + ((-t173 * t73 - t35) * t171 + (t173 * t72 - t34) * t170) * t282 + t192 + m(6) * (t24 * t90 + t25 * t91 + t48 * t67 + t49 * t68); ((t170 * t98 + t171 * t99) * (((-t99 + t159) * t173 + t238) * t170 + (-qJD(4) * t157 * t171 + t173 * t98 + t256) * t171) + t251 * t157 * t143) * t287 + (t170 * t29 - t171 * t28) * t265 + t170 * ((t170 * t55 + (t28 + t296) * t173) * t170 + (t29 * t173 + (t249 * t94 + t250 * t96) * t171 + (-t216 * qJD(4) - t217 * t173 - t56) * t170) * t171) + (t170 * t27 - t171 * t26) * t268 - t171 * ((t171 * t56 + (t27 + t295) * t173) * t171 + (t26 * t173 + (-t249 * t95 - t250 * t97) * t170 + (t218 * qJD(4) - t215 * t173 - t55) * t171) * t170) + (t15 * t3 + t48 * t91 + t49 * t90) * t286 + t223; m(6) * ((-t170 * t62 - t171 * t61) * t114 + ((-t173 * t62 - t17) * t171 + (t173 * t61 - t16) * t170) * t134) + t205; m(6) * ((-t170 * t64 - t171 * t63) * t114 + ((-t173 * t64 - t19) * t171 + (t173 * t63 - t18) * t170) * t134) + t205; m(6) * ((-t170 * t68 - t171 * t67) * t114 + ((-t173 * t68 - t25) * t171 + (t173 * t67 - t24) * t170) * t134) + t205; m(6) * (t8 * t15 + t50 * t3 + (-t170 * t90 - t171 * t91) * t114 + ((-t173 * t90 - t48) * t171 + (t173 * t91 - t49) * t170) * t134) + t223; (t114 * t134 * t251 + t50 * t8) * t286 + t223;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
