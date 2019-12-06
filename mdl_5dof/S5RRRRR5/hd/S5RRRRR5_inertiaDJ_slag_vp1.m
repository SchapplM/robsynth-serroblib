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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:58:09
% EndTime: 2019-12-05 18:58:17
% DurationCPUTime: 3.80s
% Computational Cost: add. (12576->404), mult. (7564->560), div. (0->0), fcn. (5570->10), ass. (0->238)
t186 = qJ(1) + qJ(2);
t182 = qJ(3) + t186;
t173 = sin(t182);
t185 = qJ(4) + qJ(5);
t180 = cos(t185);
t289 = rSges(6,1) * t180;
t174 = cos(t182);
t178 = sin(t185);
t287 = rSges(6,2) * t178;
t316 = t174 * rSges(6,3) + t173 * t287;
t86 = -t173 * t289 + t316;
t189 = cos(qJ(4));
t256 = qJD(4) * t189;
t184 = qJD(1) + qJD(2);
t177 = qJD(3) + t184;
t187 = sin(qJ(4));
t267 = t177 * t187;
t320 = t173 * t256 + t174 * t267;
t183 = qJD(4) + qJD(5);
t263 = t180 * t183;
t265 = t178 * t183;
t319 = rSges(6,1) * t265 + rSges(6,2) * t263;
t282 = Icges(5,4) * t189;
t217 = -Icges(5,2) * t187 + t282;
t92 = Icges(5,6) * t174 - t217 * t173;
t283 = Icges(5,4) * t187;
t219 = Icges(5,1) * t189 - t283;
t204 = t219 * t173;
t94 = Icges(5,5) * t174 - t204;
t222 = t187 * t92 - t189 * t94;
t318 = t222 * t174;
t280 = Icges(6,4) * t180;
t216 = -Icges(6,2) * t178 + t280;
t82 = Icges(6,6) * t174 - t216 * t173;
t281 = Icges(6,4) * t178;
t218 = Icges(6,1) * t180 - t281;
t203 = t218 * t173;
t84 = Icges(6,5) * t174 - t203;
t227 = t178 * t82 - t180 * t84;
t317 = t227 * t174;
t135 = Icges(6,1) * t178 + t280;
t271 = t135 * t183;
t315 = -Icges(6,5) * t177 + t271;
t134 = Icges(6,2) * t180 + t281;
t272 = t134 * t183;
t314 = -Icges(6,6) * t177 + t272;
t133 = Icges(6,5) * t178 + Icges(6,6) * t180;
t313 = -Icges(6,3) * t177 + t133 * t183;
t312 = -t319 * t174 + t86 * t177;
t157 = Icges(5,5) * t187 + Icges(5,6) * t189;
t311 = -Icges(5,3) * t177 + t157 * qJD(4);
t158 = Icges(5,2) * t189 + t283;
t310 = -Icges(5,6) * t177 + t158 * qJD(4);
t159 = Icges(5,1) * t187 + t282;
t309 = -Icges(5,5) * t177 + t159 * qJD(4);
t140 = t217 * qJD(4);
t141 = t219 * qJD(4);
t308 = t140 * t187 - t141 * t189 - t157 * t177 + (t158 * t189 + t159 * t187) * qJD(4);
t110 = t216 * t183;
t111 = t218 * t183;
t307 = (t110 + t271) * t178 + (-t111 + t272) * t180 - t133 * t177;
t306 = 2 * m(3);
t305 = 2 * m(4);
t304 = 2 * m(5);
t303 = 2 * m(6);
t302 = t173 / 0.2e1;
t301 = t174 / 0.2e1;
t300 = -rSges(5,3) - pkin(8);
t288 = rSges(5,2) * t187;
t290 = rSges(5,1) * t189;
t145 = (-t288 + t290) * qJD(4);
t299 = m(5) * t145;
t164 = rSges(5,1) * t187 + rSges(5,2) * t189;
t298 = m(5) * t164;
t190 = cos(qJ(1));
t297 = pkin(1) * t190;
t179 = sin(t186);
t296 = pkin(2) * t179;
t181 = cos(t186);
t295 = pkin(2) * t181;
t294 = pkin(8) * t173;
t169 = t174 * pkin(8);
t188 = sin(qJ(1));
t293 = t188 * pkin(1);
t175 = pkin(4) * t189 + pkin(3);
t292 = pkin(3) - t175;
t236 = t173 * t292;
t191 = -pkin(9) - pkin(8);
t268 = t174 * t191;
t291 = t169 - t236 + t268 - t86;
t286 = rSges(6,3) * t173;
t285 = pkin(1) * qJD(1);
t284 = t173 * rSges(5,3);
t168 = t174 * rSges(5,3);
t270 = t173 * t177;
t269 = t174 * t177;
t266 = t177 * t191;
t264 = t179 * t184;
t262 = t181 * t184;
t257 = qJD(4) * t187;
t243 = t173 * t257;
t261 = -pkin(4) * t243 - t173 * t266;
t104 = rSges(4,1) * t270 + rSges(4,2) * t269;
t155 = t173 * t288;
t259 = t155 + t168;
t113 = rSges(3,1) * t264 + rSges(3,2) * t262;
t258 = t173 ^ 2 + t174 ^ 2;
t255 = pkin(2) * t262;
t252 = t190 * t285;
t251 = t173 * t290;
t148 = t174 * t287;
t241 = t174 * t257;
t248 = -pkin(4) * t241 - t174 * t266 - t175 * t270;
t246 = t177 * t148 + t173 * t319;
t240 = t174 * t256;
t245 = -rSges(5,1) * t241 - rSges(5,2) * t240 - t177 * t251;
t244 = rSges(5,1) * t243 + rSges(5,2) * t320;
t163 = pkin(2) * t264;
t78 = t163 + t104;
t239 = -t270 / 0.2e1;
t238 = t269 / 0.2e1;
t237 = -pkin(3) - t290;
t136 = rSges(6,1) * t178 + rSges(6,2) * t180;
t235 = pkin(4) * t187 + t136;
t201 = t216 * t177;
t43 = -t173 * t201 - t314 * t174;
t85 = Icges(6,5) * t173 + t218 * t174;
t234 = t183 * t85 + t43;
t44 = t314 * t173 - t174 * t201;
t233 = t183 * t84 + t44;
t45 = -t315 * t174 - t177 * t203;
t83 = Icges(6,6) * t173 + t216 * t174;
t232 = -t183 * t83 + t45;
t46 = t315 * t173 - t218 * t269;
t231 = t183 * t82 - t46;
t138 = -rSges(3,1) * t181 + t179 * rSges(3,2);
t123 = -rSges(4,1) * t174 + t173 * rSges(4,2);
t230 = -t175 - t289;
t114 = -rSges(3,1) * t262 + rSges(3,2) * t264;
t105 = -rSges(4,1) * t269 + rSges(4,2) * t270;
t137 = -rSges(3,1) * t179 - rSges(3,2) * t181;
t122 = -rSges(4,1) * t173 - rSges(4,2) * t174;
t226 = t178 * t83 - t180 * t85;
t223 = t187 * t94 + t189 * t92;
t93 = Icges(5,6) * t173 + t217 * t174;
t95 = Icges(5,5) * t173 + t219 * t174;
t221 = t187 * t95 + t189 * t93;
t220 = t187 * t93 - t189 * t95;
t215 = Icges(5,5) * t189 - Icges(5,6) * t187;
t214 = Icges(6,5) * t180 - Icges(6,6) * t178;
t213 = t134 * t178 - t135 * t180;
t211 = t158 * t187 - t159 * t189;
t80 = Icges(6,3) * t174 - t214 * t173;
t20 = t227 * t173 + t174 * t80;
t207 = t226 * t173;
t81 = Icges(6,3) * t173 + t214 * t174;
t21 = t174 * t81 + t207;
t22 = t173 * t80 - t317;
t23 = t173 * t81 - t226 * t174;
t199 = t214 * t177;
t41 = -t173 * t199 - t313 * t174;
t42 = t313 * t173 - t174 * t199;
t210 = -(t173 * t21 + t174 * t20) * t270 + t173 * ((t173 * t41 + (-t22 + t207) * t177) * t173 + (t23 * t177 + (-t178 * t44 + t180 * t46 - t82 * t263 - t84 * t265) * t174 + (t227 * t177 - t234 * t178 + t232 * t180 + t42) * t173) * t174) + t174 * ((t174 * t42 + (t21 + t317) * t177) * t174 + (-t20 * t177 + (t178 * t43 - t180 * t45 + t83 * t263 + t85 * t265) * t173 + (t226 * t177 + t233 * t178 + t231 * t180 + t41) * t174) * t173) + (t173 * t23 + t174 * t22) * t269;
t108 = t123 - t295;
t209 = -t174 * t289 - t286;
t206 = t220 * t173;
t196 = t213 * t177 + t214 * t183;
t205 = (t196 * t173 - t307 * t174 + t232 * t178 + t234 * t180) * t302 + (t307 * t173 + t196 * t174 - t231 * t178 + t233 * t180) * t301 + (t133 * t174 + t213 * t173 + t178 * t84 + t180 * t82) * t239 + (t133 * t173 - t213 * t174 + t178 * t85 + t180 * t83) * t238;
t202 = t217 * t177;
t200 = t215 * t177;
t107 = t122 - t296;
t79 = t105 - t255;
t197 = t230 * t174 - t286;
t195 = t215 * qJD(4) + t211 * t177;
t71 = t237 * t173 + t169 + t259;
t194 = t300 * t173 + t237 * t174;
t156 = t174 * t288;
t72 = t156 + t194;
t165 = t173 * t191;
t68 = t148 + t165 + t197;
t67 = t230 * t173 - t268 + t316;
t69 = t71 - t296;
t70 = t72 - t295;
t64 = t68 - t295;
t63 = t67 - t296;
t146 = pkin(3) * t270;
t34 = t146 + (t300 * t174 - t155) * t177 - t245;
t193 = t205 + (-t220 * qJD(4) + t195 * t173 - t308 * t174 + t187 * (-t309 * t174 - t177 * t204) + t189 * (-t173 * t202 - t310 * t174)) * t302 + (-t222 * qJD(4) + t308 * t173 + t195 * t174 + t187 * (t309 * t173 - t219 * t269) + t189 * (t310 * t173 - t174 * t202)) * t301 + (t157 * t174 + t211 * t173 + t223) * t239 + (t157 * t173 - t211 * t174 + t221) * t238;
t192 = t180 * t110 + t178 * t111 - t134 * t265 + t135 * t263 + t189 * t140 + t187 * t141 - t158 * t257 + t159 * t256;
t24 = -t248 - t312;
t35 = t194 * t177 + t244;
t32 = t163 + t34;
t18 = t163 + t24;
t25 = t197 * t177 + t246 - t261;
t33 = t35 - t255;
t19 = t25 - t255;
t176 = t188 * t285;
t116 = t138 - t297;
t115 = t137 - t293;
t112 = (-t287 + t289) * t183;
t101 = t114 - t252;
t100 = t176 + t113;
t99 = t108 - t297;
t98 = t107 - t293;
t97 = t174 * t290 - t156 + t284;
t96 = -t251 + t259;
t91 = Icges(5,3) * t173 + t215 * t174;
t90 = Icges(5,3) * t174 - t215 * t173;
t89 = t235 * t174;
t88 = t235 * t173;
t87 = -t148 - t209;
t77 = -t292 * t174 - t165 - t294;
t75 = t174 * t87;
t74 = t79 - t252;
t73 = t176 + t78;
t66 = t70 - t297;
t65 = t69 - t293;
t62 = t64 - t297;
t61 = t63 - t293;
t56 = t311 * t173 - t174 * t200;
t55 = -t173 * t200 - t311 * t174;
t50 = t209 * t177 + t246;
t49 = -t173 * t86 + t75;
t48 = pkin(4) * t320 + t112 * t173 + t136 * t269;
t47 = t136 * t270 - t112 * t174 + (t173 * t267 - t240) * pkin(4);
t38 = t174 * t312;
t31 = t33 - t252;
t30 = t176 + t32;
t29 = t173 * t91 - t220 * t174;
t28 = t173 * t90 - t318;
t27 = t174 * t91 + t206;
t26 = t222 * t173 + t174 * t90;
t17 = t19 - t252;
t16 = t176 + t18;
t15 = t291 * t173 + t174 * t77 + t75;
t8 = -t86 * t269 + t38 + (-t177 * t87 - t50) * t173;
t3 = t174 * (t146 + t248) + t38 + (-t50 + t261) * t173 + ((-t77 - t87 - t294) * t173 + (-t236 + t291 - t169) * t174) * t177;
t1 = [(t100 * t116 + t101 * t115) * t306 + (t73 * t99 + t74 * t98) * t305 + (t30 * t66 + t31 * t65) * t304 + (t16 * t62 + t17 * t61) * t303 + t192; m(3) * (t100 * t138 + t101 * t137 + t113 * t116 + t114 * t115) + m(4) * (t107 * t74 + t108 * t73 + t78 * t99 + t79 * t98) + m(5) * (t30 * t70 + t31 * t69 + t32 * t66 + t33 * t65) + m(6) * (t16 * t64 + t17 * t63 + t18 * t62 + t19 * t61) + t192; (t18 * t64 + t19 * t63) * t303 + (t32 * t70 + t33 * t69) * t304 + (t107 * t79 + t108 * t78) * t305 + (t113 * t138 + t114 * t137) * t306 + t192; m(4) * (t104 * t99 + t105 * t98 + t122 * t74 + t123 * t73) + m(5) * (t30 * t72 + t71 * t31 + t34 * t66 + t35 * t65) + m(6) * (t16 * t68 + t17 * t67 + t24 * t62 + t25 * t61) + t192; m(6) * (t18 * t68 + t19 * t67 + t24 * t64 + t25 * t63) + m(5) * (t32 * t72 + t71 * t33 + t34 * t70 + t35 * t69) + m(4) * (t104 * t108 + t105 * t107 + t122 * t79 + t123 * t78) + t192; (t24 * t68 + t25 * t67) * t303 + (t34 * t72 + t71 * t35) * t304 + (t104 * t123 + t105 * t122) * t305 + t192; m(6) * (t16 * t88 - t17 * t89 + t47 * t61 + t48 * t62) + t193 + ((t177 * t66 - t31) * t174 + (t177 * t65 + t30) * t173) * t298 + (t173 * t66 - t174 * t65) * t299; m(6) * (t18 * t88 - t19 * t89 + t47 * t63 + t48 * t64) + t193 + ((t177 * t70 - t33) * t174 + (t177 * t69 + t32) * t173) * t298 + (t173 * t70 - t174 * t69) * t299; m(6) * (t24 * t88 - t25 * t89 + t47 * t67 + t48 * t68) + t193 + ((t177 * t72 - t35) * t174 + (t177 * t71 + t34) * t173) * t298 + (t173 * t72 - t174 * t71) * t299; ((-t173 * t96 + t174 * t97) * (t174 * t245 - t173 * t244 + ((-t96 + t168) * t174 + (t284 - t97 + (t288 + t290) * t174) * t173) * t177) + t258 * t164 * t145) * t304 - (t173 * t27 + t174 * t26) * t270 + t174 * ((t174 * t56 + (t27 + t318) * t177) * t174 + (-t26 * t177 + (t256 * t93 + t257 * t95) * t173 + (t223 * qJD(4) + t220 * t177 + t55) * t174) * t173) + (t173 * t29 + t174 * t28) * t269 + t173 * ((t173 * t55 + (-t28 + t206) * t177) * t173 + (t29 * t177 + (-t256 * t92 - t257 * t94) * t174 + (-t221 * qJD(4) + t222 * t177 + t56) * t173) * t174) + (t15 * t3 - t47 * t89 + t48 * t88) * t303 + t210; m(6) * ((t173 * t62 - t174 * t61) * t112 + ((t177 * t62 - t17) * t174 + (t177 * t61 + t16) * t173) * t136) + t205; m(6) * ((t173 * t64 - t174 * t63) * t112 + ((t177 * t64 - t19) * t174 + (t177 * t63 + t18) * t173) * t136) + t205; m(6) * ((t173 * t68 - t174 * t67) * t112 + ((t177 * t68 - t25) * t174 + (t177 * t67 + t24) * t173) * t136) + t205; m(6) * (t8 * t15 + t49 * t3 + (t173 * t88 + t174 * t89) * t112 + ((t177 * t88 - t47) * t174 + (-t177 * t89 + t48) * t173) * t136) + t210; (t112 * t136 * t258 + t49 * t8) * t303 + t210;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
