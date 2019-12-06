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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:46
% DurationCPUTime: 4.37s
% Computational Cost: add. (9454->382), mult. (6258->541), div. (0->0), fcn. (4706->10), ass. (0->234)
t181 = qJ(1) + qJ(2);
t174 = sin(t181);
t175 = cos(t181);
t296 = rSges(4,2) * sin(pkin(9));
t323 = rSges(4,3) + qJ(3);
t326 = -t174 * t296 - t323 * t175;
t183 = cos(pkin(9));
t325 = rSges(4,1) * t183 + pkin(2);
t180 = qJD(1) + qJD(2);
t74 = -t174 * t323 + (-t325 + t296) * t175;
t324 = t180 * t74;
t178 = pkin(9) + qJ(4);
t171 = sin(t178);
t172 = cos(t178);
t258 = qJD(4) * t174;
t269 = t175 * t180;
t322 = t171 * t269 + t172 * t258;
t173 = qJ(5) + t178;
t166 = cos(t173);
t179 = qJD(4) + qJD(5);
t271 = t166 * t179;
t165 = sin(t173);
t272 = t165 * t179;
t321 = rSges(6,1) * t272 + rSges(6,2) * t271;
t284 = Icges(5,4) * t172;
t212 = -Icges(5,2) * t171 + t284;
t90 = Icges(5,6) * t175 - t174 * t212;
t285 = Icges(5,4) * t171;
t214 = Icges(5,1) * t172 - t285;
t92 = Icges(5,5) * t175 - t174 * t214;
t216 = t171 * t90 - t172 * t92;
t320 = t175 * t216;
t282 = Icges(6,4) * t166;
t211 = -Icges(6,2) * t165 + t282;
t197 = t211 * t174;
t80 = Icges(6,6) * t175 - t197;
t283 = Icges(6,4) * t165;
t213 = Icges(6,1) * t166 - t283;
t199 = t213 * t174;
t82 = Icges(6,5) * t175 - t199;
t218 = t165 * t80 - t166 * t82;
t319 = t175 * t218;
t294 = rSges(6,2) * t165;
t263 = t175 * rSges(6,3) + t174 * t294;
t117 = Icges(6,1) * t165 + t282;
t273 = t117 * t179;
t318 = -Icges(6,5) * t180 + t273;
t116 = Icges(6,2) * t166 + t283;
t274 = t116 * t179;
t317 = -Icges(6,6) * t180 + t274;
t115 = Icges(6,5) * t165 + Icges(6,6) * t166;
t316 = -Icges(6,3) * t180 + t115 * t179;
t128 = Icges(5,5) * t171 + Icges(5,6) * t172;
t315 = -Icges(5,3) * t180 + qJD(4) * t128;
t129 = Icges(5,2) * t172 + t285;
t314 = -Icges(5,6) * t180 + qJD(4) * t129;
t130 = Icges(5,1) * t171 + t284;
t313 = -Icges(5,5) * t180 + qJD(4) * t130;
t112 = t212 * qJD(4);
t113 = t214 * qJD(4);
t312 = t112 * t171 - t113 * t172 - t128 * t180 + (t129 * t172 + t130 * t171) * qJD(4);
t100 = t211 * t179;
t101 = t213 * t179;
t311 = t165 * (t100 + t273) + t166 * (-t101 + t274) - t115 * t180;
t310 = 2 * m(3);
t309 = 2 * m(4);
t308 = 2 * m(5);
t307 = 2 * m(6);
t306 = t174 / 0.2e1;
t305 = t175 / 0.2e1;
t295 = rSges(5,2) * t171;
t298 = rSges(5,1) * t172;
t114 = (-t295 + t298) * qJD(4);
t304 = m(5) * t114;
t131 = rSges(5,1) * t171 + rSges(5,2) * t172;
t303 = m(5) * t131;
t185 = sin(qJ(1));
t302 = pkin(1) * t185;
t186 = cos(qJ(1));
t301 = pkin(1) * t186;
t184 = -pkin(7) - qJ(3);
t167 = t183 * pkin(3) + pkin(2);
t177 = -pkin(8) + t184;
t136 = pkin(4) * t172 + t167;
t265 = t136 - t167;
t221 = t174 * t265;
t297 = rSges(6,1) * t166;
t255 = t174 * t297;
t84 = -t255 + t263;
t300 = -(-t177 + t184) * t175 + t221 - t84;
t293 = pkin(1) * qJD(1);
t292 = t174 * rSges(5,3);
t291 = t174 * rSges(6,3);
t163 = t175 * rSges(5,3);
t270 = t174 * t180;
t268 = t177 * t180;
t267 = t180 * t184;
t266 = t167 * t270 + t175 * t267;
t245 = t171 * t258;
t264 = -pkin(4) * t245 - t174 * t268;
t262 = t174 * t295 + t163;
t107 = rSges(3,1) * t270 + rSges(3,2) * t269;
t261 = t174 ^ 2 + t175 ^ 2;
t260 = qJD(4) * t171;
t259 = qJD(4) * t172;
t257 = qJD(4) * t175;
t253 = t186 * t293;
t252 = t174 * t298;
t144 = t175 * t294;
t244 = t171 * t257;
t250 = -pkin(4) * t244 - t136 * t270 - t175 * t268;
t249 = t180 * t144 + t321 * t174;
t248 = -t321 * t175 - t180 * t255;
t242 = t172 * t257;
t247 = -rSges(5,1) * t244 - rSges(5,2) * t242 - t180 * t252;
t246 = rSges(5,1) * t245 + t322 * rSges(5,2);
t241 = -t270 / 0.2e1;
t240 = t269 / 0.2e1;
t118 = rSges(6,1) * t165 + rSges(6,2) * t166;
t238 = pkin(4) * t171 + t118;
t168 = t185 * t293;
t18 = -rSges(6,3) * t269 + (-t180 * t294 - qJD(3)) * t174 - t248 - t250;
t16 = t168 + t18;
t222 = -t136 - t297;
t63 = t174 * t222 - t175 * t177 + t263;
t59 = t63 - t302;
t237 = t180 * t59 + t16;
t159 = qJD(3) * t175;
t192 = t175 * t222 - t291;
t19 = t180 * t192 + t159 + t249 - t264;
t17 = t19 - t253;
t157 = t174 * t177;
t64 = t144 + t157 + t192;
t60 = t64 - t301;
t236 = t180 * t60 - t17;
t45 = -t175 * t317 - t180 * t197;
t83 = Icges(6,5) * t174 + t175 * t213;
t235 = t179 * t83 + t45;
t46 = t174 * t317 - t211 * t269;
t234 = t179 * t82 + t46;
t47 = -t175 * t318 - t180 * t199;
t81 = Icges(6,6) * t174 + t175 * t211;
t233 = -t179 * t81 + t47;
t48 = t174 * t318 - t213 * t269;
t232 = -t179 * t80 + t48;
t231 = t180 * t63 + t18;
t230 = t180 * t64 - t19;
t30 = -rSges(5,3) * t269 + (-t180 * t295 - qJD(3)) * t174 - t247 + t266;
t28 = t168 + t30;
t223 = -t167 - t298;
t67 = t174 * t223 - t175 * t184 + t262;
t65 = t67 - t302;
t229 = t180 * t65 + t28;
t147 = t174 * t267;
t193 = t175 * t223 - t292;
t31 = t180 * t193 + t147 + t159 + t246;
t29 = t31 - t253;
t150 = t175 * t295;
t158 = t174 * t184;
t68 = t150 + t158 + t193;
t66 = t68 - t301;
t228 = t180 * t66 - t29;
t227 = t180 * t67 + t30;
t226 = t180 * t68 - t31;
t102 = (-t294 + t297) * t179;
t38 = t118 * t270 - t102 * t175 + (t171 * t270 - t242) * pkin(4);
t76 = t238 * t174;
t225 = t180 * t76 - t38;
t39 = t322 * pkin(4) + t102 * t174 + t118 * t269;
t77 = t238 * t175;
t224 = -t180 * t77 + t39;
t133 = -t175 * rSges(3,1) + t174 * rSges(3,2);
t108 = -rSges(3,1) * t269 + rSges(3,2) * t270;
t132 = -rSges(3,1) * t174 - rSges(3,2) * t175;
t217 = t165 * t81 - t166 * t83;
t91 = Icges(5,6) * t174 + t175 * t212;
t93 = Icges(5,5) * t174 + t175 * t214;
t215 = t171 * t91 - t172 * t93;
t210 = Icges(5,5) * t172 - Icges(5,6) * t171;
t209 = Icges(6,5) * t166 - Icges(6,6) * t165;
t208 = t116 * t165 - t117 * t166;
t206 = t129 * t171 - t130 * t172;
t78 = Icges(6,3) * t175 - t174 * t209;
t20 = t174 * t218 + t175 * t78;
t202 = t217 * t174;
t79 = Icges(6,3) * t174 + t175 * t209;
t21 = t175 * t79 + t202;
t22 = t174 * t78 - t319;
t23 = t174 * t79 - t175 * t217;
t195 = t209 * t180;
t43 = -t174 * t195 - t175 * t316;
t44 = t174 * t316 - t175 * t195;
t205 = -t270 * (t174 * t21 + t175 * t20) + t174 * ((t174 * t43 + (-t22 + t202) * t180) * t174 + (t23 * t180 + (-t165 * t46 + t166 * t48 - t271 * t80 - t272 * t82) * t175 + (t44 + (-t180 * t82 + t233) * t166 + (t180 * t80 - t235) * t165) * t174) * t175) + t175 * ((t175 * t44 + (t21 + t319) * t180) * t175 + (-t20 * t180 + (t165 * t45 - t166 * t47 + t271 * t81 + t272 * t83) * t174 + (t43 + (-t180 * t83 - t232) * t166 + (t180 * t81 + t234) * t165) * t175) * t174) + (t174 * t23 + t175 * t22) * t269;
t204 = -t175 * t297 - t291;
t191 = t209 * t179 + t180 * t208;
t203 = (t165 * t233 + t166 * t235 + t191 * t174 - t175 * t311) * t306 + (t165 * t232 + t166 * t234 + t174 * t311 + t191 * t175) * t305 + (t115 * t175 + t165 * t82 + t166 * t80 + t174 * t208) * t241 + (t115 * t174 + t165 * t83 + t166 * t81 - t175 * t208) * t240;
t201 = t215 * t174;
t200 = t214 * t180;
t198 = t212 * t180;
t196 = t210 * t180;
t190 = t210 * qJD(4) + t180 * t206;
t73 = -t174 * t325 - t326;
t188 = t166 * t100 + t165 * t101 + t172 * t112 + t171 * t113 - t116 * t272 + t117 * t271 - t129 * t260 + t130 * t259;
t187 = t203 + (t190 * t174 - t175 * t312 - qJD(4) * t215 + t171 * (-t174 * t200 - t175 * t313) + t172 * (-t174 * t198 - t175 * t314)) * t306 + (t174 * t312 + t190 * t175 - qJD(4) * t216 + t171 * (t174 * t313 - t175 * t200) + t172 * (t174 * t314 - t175 * t198)) * t305 + (t128 * t175 + t171 * t92 + t172 * t90 + t174 * t206) * t241 + (t128 * t174 + t171 * t93 + t172 * t91 - t175 * t206) * t240;
t62 = t159 + t324;
t61 = -qJD(3) * t174 + t326 * t180 + t325 * t270;
t110 = t133 - t301;
t109 = t132 - t302;
t98 = t108 - t253;
t97 = t168 + t107;
t95 = t175 * t298 - t150 + t292;
t94 = -t252 + t262;
t89 = Icges(5,3) * t174 + t175 * t210;
t88 = Icges(5,3) * t175 - t174 * t210;
t85 = -t144 - t204;
t75 = t175 * t85;
t72 = t74 - t301;
t71 = t73 - t302;
t70 = t175 * t265 - t157 + t158;
t54 = t174 * t315 - t175 * t196;
t53 = -t174 * t196 - t175 * t315;
t52 = t62 - t253;
t51 = t168 + t61;
t50 = t180 * t204 + t249;
t49 = -t174 * t84 + t75;
t42 = t175 * (t180 * t263 + t248);
t27 = t174 * t89 - t175 * t215;
t26 = t174 * t88 - t320;
t25 = t175 * t89 + t201;
t24 = t174 * t216 + t175 * t88;
t15 = t174 * t300 + t175 * t70 + t75;
t10 = -t84 * t269 + t42 + (-t180 * t85 - t50) * t174;
t3 = t175 * (t250 + t266) + t42 + (t147 - t50 + t264) * t174 + ((-t70 - t85) * t174 + (t221 + t300) * t175) * t180;
t1 = [(t109 * t98 + t110 * t97) * t310 + (t51 * t72 + t52 * t71) * t309 + (t28 * t66 + t29 * t65) * t308 + (t16 * t60 + t17 * t59) * t307 + t188; m(3) * (t107 * t110 + t108 * t109 + t132 * t98 + t133 * t97) + m(4) * (t51 * t74 + t52 * t73 + t61 * t72 + t62 * t71) + m(5) * (t28 * t68 + t29 * t67 + t30 * t66 + t31 * t65) + m(6) * (t16 * t64 + t17 * t63 + t18 * t60 + t19 * t59) + t188; (t18 * t64 + t19 * t63) * t307 + (t30 * t68 + t31 * t67) * t308 + (t61 * t74 + t62 * t73) * t309 + (t107 * t133 + t108 * t132) * t310 + t188; m(4) * ((t180 * t71 + t51) * t175 + (-t180 * t72 + t52) * t174) + m(5) * (-t174 * t228 + t175 * t229) + m(6) * (-t174 * t236 + t175 * t237); m(6) * (-t174 * t230 + t175 * t231) + m(5) * (-t174 * t226 + t175 * t227) + m(4) * ((t180 * t73 + t61) * t175 + (t62 - t324) * t174); 0; t187 + m(6) * (t16 * t76 - t17 * t77 + t38 * t59 + t39 * t60) + (t174 * t229 + t175 * t228) * t303 + (t174 * t66 - t175 * t65) * t304; t187 + m(6) * (t18 * t76 - t19 * t77 + t38 * t63 + t39 * t64) + (t174 * t227 + t175 * t226) * t303 + (t174 * t68 - t175 * t67) * t304; m(6) * (-t174 * t225 + t175 * t224); ((-t174 * t94 + t175 * t95) * (t175 * t247 - t174 * t246 + ((-t94 + t163) * t175 + (t292 - t95 + (t295 + t298) * t175) * t174) * t180) + t261 * t131 * t114) * t308 - (t174 * t25 + t24 * t175) * t270 + t175 * ((t175 * t54 + (t25 + t320) * t180) * t175 + (-t24 * t180 + (t259 * t91 + t260 * t93) * t174 + (t53 + (qJD(4) * t90 - t180 * t93) * t172 + (qJD(4) * t92 + t180 * t91) * t171) * t175) * t174) + (t27 * t174 + t175 * t26) * t269 + t174 * ((t174 * t53 + (-t26 + t201) * t180) * t174 + (t27 * t180 + (-t259 * t90 - t92 * t260) * t175 + (t54 + (-qJD(4) * t91 - t180 * t92) * t172 + (-qJD(4) * t93 + t180 * t90) * t171) * t174) * t175) + (t15 * t3 - t38 * t77 + t39 * t76) * t307 + t205; m(6) * ((t174 * t60 - t175 * t59) * t102 + (t174 * t237 + t175 * t236) * t118) + t203; m(6) * ((t174 * t64 - t175 * t63) * t102 + (t174 * t231 + t175 * t230) * t118) + t203; 0; m(6) * (t10 * t15 + t3 * t49 + (t174 * t76 + t175 * t77) * t102 + (t174 * t224 + t175 * t225) * t118) + t205; (t102 * t118 * t261 + t10 * t49) * t307 + t205;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
