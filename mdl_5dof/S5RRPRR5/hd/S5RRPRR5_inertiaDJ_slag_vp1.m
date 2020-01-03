% Calculate time derivative of joint inertia matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:32
% EndTime: 2020-01-03 12:03:40
% DurationCPUTime: 3.99s
% Computational Cost: add. (9454->393), mult. (6258->552), div. (0->0), fcn. (4706->10), ass. (0->232)
t286 = rSges(4,2) * sin(pkin(9));
t314 = pkin(2) - t286;
t313 = rSges(4,3) + qJ(3);
t182 = pkin(9) + qJ(4);
t175 = qJ(5) + t182;
t165 = sin(t175);
t284 = rSges(6,2) * t165;
t166 = cos(t175);
t287 = rSges(6,1) * t166;
t312 = -t284 + t287;
t173 = sin(t182);
t285 = rSges(5,2) * t173;
t174 = cos(t182);
t289 = rSges(5,1) * t174;
t311 = -t285 + t289;
t184 = qJD(1) + qJD(2);
t183 = qJD(4) + qJD(5);
t104 = t312 * t183;
t185 = qJ(1) + qJ(2);
t177 = cos(t185);
t176 = sin(t185);
t259 = t176 * t184;
t240 = t173 * t259;
t283 = rSges(6,2) * t166;
t288 = rSges(6,1) * t165;
t122 = t283 + t288;
t241 = t122 * t259;
t247 = qJD(4) * t174;
t38 = -t241 + t104 * t177 + (t247 * t177 - t240) * pkin(4);
t235 = pkin(4) * t173 + t122;
t78 = t235 * t176;
t310 = -t184 * t78 - t38;
t275 = Icges(5,4) * t174;
t214 = -Icges(5,2) * t173 + t275;
t93 = -Icges(5,6) * t176 - t177 * t214;
t276 = Icges(5,4) * t173;
t216 = Icges(5,1) * t174 - t276;
t95 = -Icges(5,5) * t176 - t177 * t216;
t217 = t173 * t93 - t174 * t95;
t309 = t176 * t217;
t273 = Icges(6,4) * t166;
t213 = -Icges(6,2) * t165 + t273;
t83 = -Icges(6,6) * t176 - t177 * t213;
t274 = Icges(6,4) * t165;
t215 = Icges(6,1) * t166 - t274;
t85 = -Icges(6,5) * t176 - t177 * t215;
t219 = t165 * t83 - t166 * t85;
t308 = t176 * t219;
t188 = -pkin(7) - qJ(3);
t181 = -pkin(8) + t188;
t249 = t181 - t188;
t307 = t249 * t176;
t132 = t176 * rSges(3,1) + t177 * rSges(3,2);
t121 = Icges(6,1) * t165 + t273;
t263 = t121 * t183;
t306 = -Icges(6,5) * t184 + t263;
t120 = Icges(6,2) * t166 + t274;
t264 = t120 * t183;
t305 = -Icges(6,6) * t184 + t264;
t119 = Icges(6,5) * t165 + Icges(6,6) * t166;
t304 = -Icges(6,3) * t184 + t119 * t183;
t128 = Icges(5,5) * t173 + Icges(5,6) * t174;
t303 = -Icges(5,3) * t184 + qJD(4) * t128;
t129 = Icges(5,2) * t174 + t276;
t302 = -Icges(5,6) * t184 + qJD(4) * t129;
t130 = Icges(5,1) * t173 + t275;
t301 = -Icges(5,5) * t184 + qJD(4) * t130;
t116 = t214 * qJD(4);
t117 = t216 * qJD(4);
t300 = t116 * t173 - t117 * t174 - t128 * t184 + (t129 * t174 + t130 * t173) * qJD(4);
t102 = t213 * t183;
t103 = t215 * t183;
t299 = t165 * (t102 + t263) + t166 * (-t103 + t264) - t119 * t184;
t298 = 2 * m(3);
t297 = 2 * m(4);
t296 = 2 * m(5);
t295 = 2 * m(6);
t294 = -t176 / 0.2e1;
t293 = -t177 / 0.2e1;
t118 = t311 * qJD(4);
t292 = m(5) * t118;
t131 = rSges(5,1) * t173 + rSges(5,2) * t174;
t291 = m(5) * t131;
t187 = cos(pkin(9));
t167 = t187 * pkin(3) + pkin(2);
t290 = rSges(4,1) * t187;
t282 = pkin(1) * qJD(1);
t211 = Icges(6,5) * t166 - Icges(6,6) * t165;
t80 = -Icges(6,3) * t177 + t176 * t211;
t280 = t184 * t80;
t81 = -Icges(6,3) * t176 - t177 * t211;
t279 = t184 * t81;
t266 = t104 * t176;
t262 = t165 * t183;
t261 = t166 * t183;
t260 = t176 * t183;
t258 = t177 * t183;
t257 = t177 * t184;
t136 = pkin(4) * t174 + t167;
t256 = t176 * t136 + t177 * t181;
t244 = t184 * t284;
t255 = -rSges(6,3) * t257 - t176 * t244;
t138 = t177 * t287;
t254 = rSges(6,3) * t259 + t184 * t138;
t253 = rSges(5,2) * t240 + rSges(5,3) * t257;
t142 = t177 * t289;
t252 = rSges(5,3) * t259 + t184 * t142;
t251 = -t176 * t167 - t177 * t188;
t250 = t176 ^ 2 + t177 ^ 2;
t248 = qJD(4) * t173;
t86 = -rSges(6,3) * t177 + t312 * t176;
t87 = -t176 * rSges(6,3) + t177 * t284 - t138;
t246 = t176 * (-t260 * t288 + (-t165 * t257 - t166 * t260) * rSges(6,2) + t254) + t87 * t259 + t86 * t257;
t245 = t184 * t286;
t189 = sin(qJ(1));
t243 = t189 * t282;
t242 = pkin(4) * t248;
t154 = t177 * t290;
t82 = -Icges(6,6) * t177 + t176 * t213;
t84 = -Icges(6,5) * t177 + t176 * t215;
t220 = t165 * t82 - t166 * t84;
t20 = -t176 * t220 - t177 * t80;
t21 = -t177 * t81 - t308;
t206 = t220 * t177;
t22 = -t176 * t80 + t206;
t202 = t215 * t184;
t47 = t176 * t202 + t177 * t306;
t229 = -t183 * t83 + t47;
t23 = -t176 * t81 + t177 * t219;
t200 = t213 * t184;
t45 = t176 * t200 + t177 * t305;
t231 = t183 * t85 + t45;
t198 = t211 * t184;
t43 = t176 * t198 + t177 * t304;
t44 = -t176 * t304 + t177 * t198;
t46 = -t176 * t305 + t177 * t200;
t48 = -t176 * t306 + t177 * t202;
t239 = -t176 * ((t176 * t43 + (t22 + t308) * t184) * t176 + (-t23 * t184 + (-t165 * t46 + t166 * t48 - t261 * t82 - t262 * t84 + t280) * t177 + (t279 + t44 + (-t184 * t84 + t229) * t166 + (t184 * t82 - t231) * t165) * t176) * t177) + (-t176 * t21 - t177 * t20) * t259;
t228 = -t183 * t82 + t48;
t230 = t183 * t84 + t46;
t2 = (t177 * t44 + (-t21 + t206) * t184) * t177 + (t20 * t184 + (t165 * t45 - t166 * t47 + t261 * t83 + t262 * t85 - t279) * t176 + (-t280 + t43 + (-t184 * t85 - t228) * t166 + (t184 * t83 + t230) * t165) * t177) * t176;
t5 = -t176 * t23 - t177 * t22;
t238 = -t184 * t5 - t2;
t237 = t259 / 0.2e1;
t236 = -t257 / 0.2e1;
t157 = qJD(3) * t176;
t192 = -t122 * t183 - t181 * t184 - t242;
t18 = t157 + (-t136 - t287) * t259 + t192 * t177 - t255;
t16 = t18 - t243;
t179 = t189 * pkin(1);
t63 = t86 + t256;
t59 = t179 + t63;
t234 = -t184 * t59 - t16;
t190 = cos(qJ(1));
t170 = t190 * t282;
t107 = t136 * t257;
t19 = t107 + (-qJD(3) - t244) * t177 + t192 * t176 + t254;
t17 = t170 + t19;
t180 = t190 * pkin(1);
t112 = t177 * t136;
t64 = -t176 * t181 + t112 - t87;
t60 = t180 + t64;
t233 = t184 * t60 - t17;
t232 = -t184 * t63 - t18;
t227 = t184 * t64 - t19;
t204 = t131 * qJD(4);
t194 = -t184 * t188 - t204;
t30 = t157 + (-t167 - t289) * t259 + t194 * t177 + t253;
t28 = t30 - t243;
t96 = -rSges(5,3) * t177 + t311 * t176;
t67 = t96 - t251;
t65 = t179 + t67;
t226 = -t184 * t65 - t28;
t125 = t167 * t257;
t31 = t125 + (-t184 * t285 - qJD(3)) * t177 + t194 * t176 + t252;
t29 = t170 + t31;
t140 = t177 * t167;
t97 = -t176 * rSges(5,3) + t177 * t285 - t142;
t68 = -t176 * t188 + t140 - t97;
t66 = t180 + t68;
t225 = t184 * t66 - t29;
t224 = -t184 * t67 - t30;
t223 = t184 * t68 - t31;
t133 = t177 * rSges(3,1) - rSges(3,2) * t176;
t110 = rSges(3,1) * t257 - rSges(3,2) * t259;
t92 = -Icges(5,6) * t177 + t176 * t214;
t203 = t216 * t176;
t94 = -Icges(5,5) * t177 + t203;
t218 = t173 * t92 - t174 * t94;
t212 = Icges(5,5) * t174 - Icges(5,6) * t173;
t210 = t120 * t165 - t121 * t166;
t208 = t129 * t173 - t130 * t174;
t196 = -t211 * t183 - t184 * t210;
t207 = (t165 * t229 + t166 * t231 + t196 * t176 + t177 * t299) * t294 + (t165 * t228 + t166 * t230 - t176 * t299 + t196 * t177) * t293 + (-t119 * t177 + t165 * t84 + t166 * t82 - t176 * t210) * t237 + (-t119 * t176 + t165 * t85 + t166 * t83 + t177 * t210) * t236;
t109 = t132 * t184;
t205 = t218 * t177;
t201 = t214 * t184;
t199 = t212 * t184;
t76 = t313 * t176 + t314 * t177 + t154;
t195 = -t212 * qJD(4) - t184 * t208;
t75 = -t313 * t177 + (t290 + t314) * t176;
t61 = t176 * t245 + t157 + (-pkin(2) - t290) * t259 + t313 * t257;
t62 = t184 * t154 + pkin(2) * t257 + (-qJD(3) - t245) * t177 + t313 * t259;
t193 = t166 * t102 + t165 * t103 + t174 * t116 + t173 * t117 - t120 * t262 + t121 * t261 - t129 * t248 + t130 * t247;
t191 = t207 + (t195 * t176 + t177 * t300 - qJD(4) * t217 + t173 * (t177 * t301 + t184 * t203) + t174 * (t176 * t201 + t177 * t302)) * t294 + (-t176 * t300 + t195 * t177 - qJD(4) * t218 + t173 * (-t176 * t301 + t216 * t257) + t174 * (-t176 * t302 + t177 * t201)) * t293 + (-t128 * t177 + t173 * t94 + t174 * t92 - t176 * t208) * t237 + (-t128 * t176 + t173 * t95 + t174 * t93 + t177 * t208) * t236;
t114 = t133 + t180;
t113 = t179 + t132;
t100 = t110 + t170;
t99 = -t109 - t243;
t91 = -Icges(5,3) * t176 - t177 * t212;
t90 = -Icges(5,3) * t177 + t176 * t212;
t79 = t235 * t177;
t77 = t176 * t86;
t72 = t180 + t76;
t71 = t179 + t75;
t70 = -t112 + t140 + t307;
t69 = t251 + t256;
t54 = -t176 * t303 + t177 * t199;
t53 = t176 * t199 + t177 * t303;
t52 = t170 + t62;
t51 = t61 - t243;
t50 = t258 * t283 + (t165 * t258 + t166 * t259) * rSges(6,1) + t255;
t49 = -t177 * t87 + t77;
t39 = -t122 * t257 - t266 + (-t173 * t257 - t176 * t247) * pkin(4);
t27 = -t176 * t91 + t177 * t217;
t26 = -t176 * t90 + t205;
t25 = -t177 * t91 - t309;
t24 = -t176 * t218 - t177 * t90;
t15 = t176 * t69 + t77 + (-t70 - t87) * t177;
t10 = -t177 * t50 + t246;
t3 = t176 * (-t176 * t242 + t107 - t125) + (-t177 * t242 - t50) * t177 + ((-t249 * t177 + t69) * t177 + (t70 - t177 * (t136 - t167) - t307) * t176) * t184 + t246;
t1 = [(t100 * t113 + t114 * t99) * t298 + (t51 * t72 + t52 * t71) * t297 + (t28 * t66 + t29 * t65) * t296 + (t16 * t60 + t17 * t59) * t295 + t193; m(3) * (t100 * t132 - t109 * t114 + t110 * t113 + t133 * t99) + m(4) * (t51 * t76 + t52 * t75 + t61 * t72 + t62 * t71) + m(5) * (t28 * t68 + t29 * t67 + t30 * t66 + t31 * t65) + m(6) * (t16 * t64 + t17 * t63 + t18 * t60 + t19 * t59) + t193; (t18 * t64 + t19 * t63) * t295 + (t30 * t68 + t31 * t67) * t296 + (t61 * t76 + t62 * t75) * t297 + (-t109 * t133 + t110 * t132) * t298 + t193; m(4) * ((-t184 * t71 - t51) * t177 + (t184 * t72 - t52) * t176) + m(5) * (t176 * t225 + t177 * t226) + m(6) * (t176 * t233 + t177 * t234); m(6) * (t176 * t227 + t177 * t232) + m(5) * (t176 * t223 + t177 * t224) + m(4) * ((-t184 * t75 - t61) * t177 + (t184 * t76 - t62) * t176); 0; t191 + m(6) * (-t16 * t78 + t17 * t79 + t38 * t59 + t39 * t60) + (t176 * t226 - t177 * t225) * t291 + (-t176 * t66 + t177 * t65) * t292; t191 + m(6) * (-t18 * t78 + t19 * t79 + t38 * t63 + t39 * t64) + (t176 * t224 - t177 * t223) * t291 + (-t176 * t68 + t177 * t67) * t292; m(6) * ((-t184 * t79 - t39) * t177 + t310 * t176); ((t176 * t96 - t177 * t97) * ((-t177 * t204 + t184 * t96 + t253) * t177 + (-t176 * t204 + (t97 + (-t285 - t289) * t177) * t184 + t252) * t176) + t250 * t131 * t118) * t296 + (-t176 * t25 - t177 * t24) * t259 - t177 * ((t177 * t54 + (-t25 + t205) * t184) * t177 + (t24 * t184 + (t247 * t93 + t248 * t95) * t176 + (t53 + (qJD(4) * t92 - t184 * t95) * t174 + (qJD(4) * t94 + t184 * t93) * t173) * t177) * t176) - t176 * ((t176 * t53 + (t26 + t309) * t184) * t176 + (-t27 * t184 + (-t247 * t92 - t94 * t248) * t177 + (t54 + (-qJD(4) * t93 - t184 * t94) * t174 + (-qJD(4) * t95 + t184 * t92) * t173) * t176) * t177) + (t15 * t3 + t38 * t79 - t39 * t78) * t295 - t177 * t2 + t239 + (t176 * t27 + t177 * t26 - t5) * t257; m(6) * ((-t176 * t60 + t177 * t59) * t104 + (t176 * t234 - t177 * t233) * t122) + t207; m(6) * ((-t176 * t64 + t177 * t63) * t104 + (t176 * t232 - t177 * t227) * t122) + t207; 0; m(6) * (-t122 * t176 * t39 + t10 * t15 - t79 * t241 + t78 * t266 + t3 * t49) + (m(6) * (t104 * t79 - t310 * t122) + t238) * t177 + t239; (t104 * t122 * t250 + t10 * t49) * t295 + t238 * t177 + t239;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
