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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:13:03
% EndTime: 2020-01-03 12:13:11
% DurationCPUTime: 3.50s
% Computational Cost: add. (12576->418), mult. (7564->581), div. (0->0), fcn. (5570->10), ass. (0->236)
t190 = qJ(1) + qJ(2);
t184 = qJ(3) + t190;
t175 = sin(t184);
t191 = sin(qJ(4));
t193 = cos(qJ(4));
t176 = cos(t184);
t286 = Icges(5,4) * t193;
t220 = -Icges(5,2) * t191 + t286;
t95 = -Icges(5,6) * t175 - t176 * t220;
t287 = Icges(5,4) * t191;
t222 = Icges(5,1) * t193 - t287;
t97 = -Icges(5,5) * t175 - t176 * t222;
t223 = t191 * t95 - t193 * t97;
t310 = t175 * t223;
t189 = qJ(4) + qJ(5);
t180 = sin(t189);
t182 = cos(t189);
t284 = Icges(6,4) * t182;
t219 = -Icges(6,2) * t180 + t284;
t85 = -Icges(6,6) * t175 - t176 * t219;
t285 = Icges(6,4) * t180;
t221 = Icges(6,1) * t182 - t285;
t87 = -Icges(6,5) * t175 - t176 * t221;
t229 = t180 * t85 - t182 * t87;
t309 = t175 * t229;
t124 = t175 * rSges(4,1) + t176 * rSges(4,2);
t181 = sin(t190);
t183 = cos(t190);
t133 = t181 * rSges(3,1) + t183 * rSges(3,2);
t188 = qJD(1) + qJD(2);
t179 = qJD(3) + t188;
t131 = Icges(6,1) * t180 + t284;
t187 = qJD(4) + qJD(5);
t274 = t131 * t187;
t308 = -Icges(6,5) * t179 + t274;
t130 = Icges(6,2) * t182 + t285;
t275 = t130 * t187;
t307 = -Icges(6,6) * t179 + t275;
t129 = Icges(6,5) * t180 + Icges(6,6) * t182;
t306 = -Icges(6,3) * t179 + t129 * t187;
t153 = Icges(5,5) * t191 + Icges(5,6) * t193;
t305 = -Icges(5,3) * t179 + qJD(4) * t153;
t154 = Icges(5,2) * t193 + t287;
t304 = -Icges(5,6) * t179 + qJD(4) * t154;
t155 = Icges(5,1) * t191 + t286;
t303 = -Icges(5,5) * t179 + qJD(4) * t155;
t136 = t220 * qJD(4);
t137 = t222 * qJD(4);
t302 = t136 * t191 - t137 * t193 - t153 * t179 + (t154 * t193 + t155 * t191) * qJD(4);
t112 = t219 * t187;
t113 = t221 * t187;
t301 = t180 * (t112 + t274) + t182 * (-t113 + t275) - t129 * t179;
t300 = 2 * m(3);
t299 = 2 * m(4);
t298 = 2 * m(5);
t297 = 2 * m(6);
t296 = -t175 / 0.2e1;
t295 = -t176 / 0.2e1;
t290 = rSges(5,2) * t191;
t292 = rSges(5,1) * t193;
t143 = (-t290 + t292) * qJD(4);
t294 = m(5) * t143;
t158 = rSges(5,1) * t191 + rSges(5,2) * t193;
t293 = m(5) * t158;
t291 = rSges(6,1) * t182;
t289 = rSges(6,2) * t180;
t288 = pkin(1) * qJD(1);
t114 = (-t289 + t291) * t187;
t277 = t114 * t175;
t273 = t175 * t179;
t272 = t175 * t182;
t271 = t175 * t193;
t195 = -pkin(9) - pkin(8);
t270 = t175 * t195;
t269 = t176 * t179;
t268 = t176 * t180;
t267 = t179 * t191;
t266 = t180 * t187;
t265 = t181 * t188;
t264 = t182 * t187;
t263 = t183 * t188;
t250 = t175 * t289;
t262 = rSges(6,3) * t269 + t179 * t250;
t150 = t176 * t291;
t261 = rSges(6,3) * t273 + t150 * t179;
t247 = t175 * t267;
t260 = rSges(5,2) * t247 + rSges(5,3) * t269;
t152 = t176 * t292;
t259 = rSges(5,3) * t273 + t152 * t179;
t258 = -pkin(3) * t269 - pkin(8) * t273;
t177 = pkin(4) * t193 + pkin(3);
t257 = t175 * t177 + t176 * t195;
t256 = pkin(3) * t176 + pkin(8) * t175;
t255 = t175 ^ 2 + t176 ^ 2;
t254 = qJD(4) * t191;
t253 = qJD(4) * t193;
t246 = t179 * t268;
t88 = rSges(6,1) * t272 - rSges(6,3) * t176 - t250;
t89 = rSges(6,2) * t268 - rSges(6,3) * t175 - t150;
t252 = t175 * (-t175 * rSges(6,1) * t266 + (-t175 * t264 - t246) * rSges(6,2) + t261) + t89 * t273 + t88 * t269;
t251 = pkin(2) * t265;
t192 = sin(qJ(1));
t249 = t192 * t288;
t132 = rSges(6,1) * t180 + rSges(6,2) * t182;
t248 = t132 * t273;
t173 = pkin(2) * t181;
t109 = t173 + t124;
t84 = -Icges(6,6) * t176 + t175 * t219;
t86 = -Icges(6,5) * t176 + t175 * t221;
t230 = t180 * t84 - t182 * t86;
t217 = Icges(6,5) * t182 - Icges(6,6) * t180;
t82 = -Icges(6,3) * t176 + t175 * t217;
t20 = -t175 * t230 - t176 * t82;
t83 = -Icges(6,3) * t175 - t176 * t217;
t21 = -t176 * t83 - t309;
t212 = t230 * t176;
t22 = -t175 * t82 + t212;
t23 = -t175 * t83 + t176 * t229;
t207 = t221 * t179;
t45 = t175 * t207 + t176 * t308;
t236 = -t187 * t85 + t45;
t205 = t219 * t179;
t43 = t175 * t205 + t176 * t307;
t238 = t187 * t87 + t43;
t203 = t217 * t179;
t41 = t175 * t203 + t176 * t306;
t42 = -t175 * t306 + t176 * t203;
t44 = -t175 * t307 + t176 * t205;
t46 = -t175 * t308 + t176 * t207;
t245 = -t175 * ((t175 * t41 + (t22 + t309) * t179) * t175 + (-t23 * t179 + (t179 * t82 - t180 * t44 + t182 * t46 - t264 * t84 - t266 * t86) * t176 + (t42 + t236 * t182 - t238 * t180 + (t230 + t83) * t179) * t175) * t176) + (-t175 * t21 - t176 * t20) * t273;
t235 = t187 * t84 - t46;
t237 = t187 * t86 + t44;
t2 = (t176 * t42 + (-t21 + t212) * t179) * t176 + (t20 * t179 + (-t179 * t83 + t180 * t43 - t182 * t45 + t264 * t85 + t266 * t87) * t175 + (t41 + t235 * t182 + t237 * t180 + (t229 - t82) * t179) * t176) * t175;
t5 = -t175 * t23 - t176 * t22;
t244 = -t179 * t5 - t2;
t243 = t175 * t254;
t242 = t176 * t253;
t241 = t273 / 0.2e1;
t240 = -t269 / 0.2e1;
t239 = pkin(4) * t191 + t132;
t134 = t183 * rSges(3,1) - rSges(3,2) * t181;
t125 = rSges(4,1) * t176 - rSges(4,2) * t175;
t232 = -t176 * t177 + t270;
t116 = rSges(3,1) * t263 - rSges(3,2) * t265;
t107 = rSges(4,1) * t269 - rSges(4,2) * t273;
t231 = rSges(5,1) * t271 - t175 * t290;
t174 = pkin(2) * t183;
t110 = t125 + t174;
t94 = -Icges(5,6) * t176 + t175 * t220;
t96 = -Icges(5,5) * t176 + t175 * t222;
t226 = t191 * t96 + t193 * t94;
t225 = t191 * t94 - t193 * t96;
t224 = t191 * t97 + t193 * t95;
t218 = Icges(5,5) * t193 - Icges(5,6) * t191;
t216 = t130 * t180 - t131 * t182;
t214 = t154 * t191 - t155 * t193;
t157 = pkin(2) * t263;
t81 = t107 + t157;
t99 = -rSges(5,3) * t175 + t176 * t290 - t152;
t213 = -pkin(4) * t254 - t179 * t195;
t115 = t133 * t188;
t106 = t124 * t179;
t211 = t225 * t176;
t200 = -t179 * t216 - t187 * t217;
t210 = (t200 * t175 + t176 * t301 + t180 * t236 + t182 * t238) * t296 + (-t175 * t301 + t200 * t176 - t180 * t235 + t182 * t237) * t295 + (-t129 * t176 - t175 * t216 + t180 * t86 + t182 * t84) * t241 + (-t129 * t175 + t176 * t216 + t180 * t87 + t182 * t85) * t240;
t209 = qJD(4) * t158;
t208 = t222 * t179;
t206 = t220 * t179;
t204 = t218 * t179;
t202 = -t175 * t253 - t176 * t267;
t74 = -t99 + t256;
t199 = -qJD(4) * t218 - t179 * t214;
t70 = t174 + t74;
t67 = t88 + t257;
t68 = -t232 - t89;
t166 = t175 * pkin(3);
t73 = t166 + (-rSges(5,3) - pkin(8)) * t176 + t231;
t63 = t173 + t67;
t64 = t174 + t68;
t69 = t173 + t73;
t80 = -t106 - t251;
t198 = -t132 * t187 + t213;
t197 = t210 + (-qJD(4) * t223 + t199 * t175 + t302 * t176 + t191 * (t175 * t208 + t176 * t303) + t193 * (t175 * t206 + t176 * t304)) * t296 + (-qJD(4) * t225 - t175 * t302 + t199 * t176 + t191 * (-t175 * t303 + t176 * t208) + t193 * (-t175 * t304 + t176 * t206)) * t295 + (-t153 * t176 - t175 * t214 + t226) * t241 + (-t153 * t175 + t176 * t214 + t224) * t240;
t196 = t112 * t182 + t180 * t113 - t130 * t266 + t131 * t264 + t136 * t193 + t137 * t191 - t154 * t254 + t155 * t253;
t35 = -rSges(5,1) * t243 + rSges(5,2) * t202 - t258 + t259;
t33 = t157 + t35;
t145 = pkin(8) * t269;
t34 = -rSges(5,2) * t242 - pkin(3) * t273 + t145 + (-t176 * t254 - t179 * t271) * rSges(5,1) + t260;
t119 = t177 * t269;
t25 = -rSges(6,2) * t246 + t175 * t198 + t119 + t261;
t32 = t34 - t251;
t19 = t157 + t25;
t24 = (-t177 - t291) * t273 + t198 * t176 + t262;
t18 = t24 - t251;
t194 = cos(qJ(1));
t186 = t194 * pkin(1);
t185 = t192 * pkin(1);
t178 = t194 * t288;
t118 = t134 + t186;
t117 = t185 + t133;
t103 = t116 + t178;
t102 = -t115 - t249;
t101 = t110 + t186;
t100 = t185 + t109;
t98 = -rSges(5,3) * t176 + t231;
t93 = -Icges(5,3) * t175 - t176 * t218;
t92 = -Icges(5,3) * t176 + t175 * t218;
t91 = t239 * t176;
t90 = t239 * t175;
t79 = t232 + t256;
t78 = pkin(8) * t176 - t166 + t257;
t77 = t175 * t88;
t76 = t178 + t81;
t75 = t80 - t249;
t66 = t186 + t70;
t65 = t185 + t69;
t62 = t186 + t64;
t61 = t185 + t63;
t56 = -t175 * t305 + t176 * t204;
t55 = t175 * t204 + t176 * t305;
t50 = t176 * rSges(6,2) * t264 + (t176 * t266 + t179 * t272) * rSges(6,1) - t262;
t49 = -t176 * t89 + t77;
t48 = pkin(4) * t202 - t132 * t269 - t277;
t47 = -t248 + t114 * t176 + (t242 - t247) * pkin(4);
t31 = t178 + t33;
t30 = t32 - t249;
t29 = -t175 * t93 + t176 * t223;
t28 = -t175 * t92 + t211;
t27 = -t176 * t93 - t310;
t26 = -t175 * t225 - t176 * t92;
t17 = t178 + t19;
t16 = t18 - t249;
t15 = t175 * t78 + t77 + (-t79 - t89) * t176;
t8 = -t176 * t50 + t252;
t3 = (t176 * t213 + t179 * t78 - t145 - t50) * t176 + (-pkin(4) * t243 + t119 + (t79 - t176 * (-pkin(3) + t177) - t270) * t179 + t258) * t175 + t252;
t1 = [(t102 * t118 + t103 * t117) * t300 + (t100 * t76 + t101 * t75) * t299 + (t30 * t66 + t31 * t65) * t298 + (t16 * t62 + t17 * t61) * t297 + t196; m(3) * (t102 * t134 + t103 * t133 - t115 * t118 + t116 * t117) + m(4) * (t100 * t81 + t101 * t80 + t109 * t76 + t110 * t75) + m(5) * (t30 * t70 + t31 * t69 + t32 * t66 + t33 * t65) + m(6) * (t16 * t64 + t17 * t63 + t18 * t62 + t19 * t61) + t196; (t18 * t64 + t19 * t63) * t297 + (t32 * t70 + t33 * t69) * t298 + (t109 * t81 + t110 * t80) * t299 + (-t115 * t134 + t116 * t133) * t300 + t196; m(4) * (t100 * t107 - t101 * t106 + t124 * t76 + t125 * t75) + m(5) * (t30 * t74 + t31 * t73 + t34 * t66 + t35 * t65) + m(6) * (t16 * t68 + t17 * t67 + t24 * t62 + t25 * t61) + t196; m(6) * (t18 * t68 + t19 * t67 + t24 * t64 + t25 * t63) + m(5) * (t32 * t74 + t33 * t73 + t34 * t70 + t35 * t69) + m(4) * (-t106 * t110 + t107 * t109 + t124 * t81 + t125 * t80) + t196; (t24 * t68 + t25 * t67) * t297 + (t34 * t74 + t35 * t73) * t298 + (-t106 * t125 + t107 * t124) * t299 + t196; t197 + m(6) * (-t16 * t90 + t17 * t91 + t47 * t61 + t48 * t62) + ((-t179 * t66 + t31) * t176 + (-t179 * t65 - t30) * t175) * t293 + (-t175 * t66 + t176 * t65) * t294; t197 + m(6) * (-t18 * t90 + t19 * t91 + t47 * t63 + t48 * t64) + ((-t179 * t70 + t33) * t176 + (-t179 * t69 - t32) * t175) * t293 + (-t175 * t70 + t176 * t69) * t294; t197 + ((-t179 * t74 + t35) * t176 + (-t179 * t73 - t34) * t175) * t293 + (-t175 * t74 + t176 * t73) * t294 + m(6) * (-t24 * t90 + t25 * t91 + t47 * t67 + t48 * t68); ((t175 * t98 - t176 * t99) * ((-t176 * t209 + t179 * t98 + t260) * t176 + (-t175 * t209 + (t99 + (-t290 - t292) * t176) * t179 + t259) * t175) + t255 * t158 * t143) * t298 + (-t175 * t27 - t176 * t26) * t273 - t176 * ((t176 * t56 + (-t27 + t211) * t179) * t176 + (t26 * t179 + (t253 * t95 + t254 * t97) * t175 + (t226 * qJD(4) + t179 * t223 + t55) * t176) * t175) - t175 * ((t175 * t55 + (t28 + t310) * t179) * t175 + (-t29 * t179 + (-t253 * t94 - t254 * t96) * t176 + (-t224 * qJD(4) + t179 * t225 + t56) * t175) * t176) + (t15 * t3 + t47 * t91 - t48 * t90) * t297 - t176 * t2 + t245 + (t175 * t29 + t176 * t28 - t5) * t269; m(6) * ((-t175 * t62 + t176 * t61) * t114 + ((-t179 * t62 + t17) * t176 + (-t179 * t61 - t16) * t175) * t132) + t210; m(6) * ((-t175 * t64 + t176 * t63) * t114 + ((-t179 * t64 + t19) * t176 + (-t179 * t63 - t18) * t175) * t132) + t210; m(6) * ((-t175 * t68 + t176 * t67) * t114 + ((-t179 * t68 + t25) * t176 + (-t179 * t67 - t24) * t175) * t132) + t210; m(6) * (-t132 * t175 * t48 + t8 * t15 - t248 * t91 + t90 * t277 + t49 * t3) + (m(6) * (t114 * t91 + (t179 * t90 + t47) * t132) + t244) * t176 + t245; (t114 * t132 * t255 + t49 * t8) * t297 + t244 * t176 + t245;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
