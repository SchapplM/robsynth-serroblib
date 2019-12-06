% Calculate time derivative of joint inertia matrix for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:22
% EndTime: 2019-12-05 15:40:37
% DurationCPUTime: 8.54s
% Computational Cost: add. (4888->439), mult. (14130->663), div. (0->0), fcn. (13422->6), ass. (0->217)
t171 = sin(pkin(7));
t172 = cos(pkin(7));
t173 = sin(qJ(4));
t174 = sin(qJ(2));
t175 = cos(qJ(4));
t249 = t174 * t175;
t157 = t171 * t249 + t172 * t173;
t254 = t173 * t174;
t236 = t171 * t254;
t158 = -t172 * t175 + t236;
t176 = cos(qJ(2));
t257 = t171 * t176;
t83 = Icges(6,5) * t158 + Icges(6,6) * t257 - Icges(6,3) * t157;
t89 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t257;
t334 = t83 - t89;
t85 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t257;
t87 = Icges(6,4) * t158 + Icges(6,2) * t257 - Icges(6,6) * t157;
t325 = t85 + t87;
t91 = Icges(6,1) * t158 + Icges(6,4) * t257 - Icges(6,5) * t157;
t93 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t257;
t331 = t91 + t93;
t156 = t171 * t175 + t172 * t254;
t188 = -t171 * t173 + t172 * t249;
t255 = t172 * t176;
t82 = Icges(6,5) * t156 + Icges(6,6) * t255 - Icges(6,3) * t188;
t88 = Icges(5,4) * t156 + Icges(5,2) * t188 + Icges(5,6) * t255;
t329 = t82 - t88;
t90 = Icges(6,1) * t156 + Icges(6,4) * t255 - Icges(6,5) * t188;
t92 = Icges(5,1) * t156 + Icges(5,4) * t188 + Icges(5,5) * t255;
t332 = t90 + t92;
t330 = -t156 * t331 + t188 * t334 - t255 * t325;
t84 = Icges(5,5) * t156 + Icges(5,6) * t188 + Icges(5,3) * t255;
t86 = Icges(6,4) * t156 + Icges(6,2) * t255 - Icges(6,6) * t188;
t326 = t84 + t86;
t263 = Icges(6,5) * t173;
t205 = -Icges(6,3) * t175 + t263;
t130 = Icges(6,6) * t174 - t176 * t205;
t265 = Icges(5,4) * t173;
t208 = Icges(5,2) * t175 + t265;
t133 = Icges(5,6) * t174 - t176 * t208;
t319 = t130 - t133;
t206 = Icges(5,5) * t173 + Icges(5,6) * t175;
t131 = Icges(5,3) * t174 - t176 * t206;
t207 = Icges(6,4) * t173 - Icges(6,6) * t175;
t132 = Icges(6,2) * t174 - t176 * t207;
t327 = t131 + t132;
t262 = Icges(6,5) * t175;
t210 = Icges(6,1) * t173 - t262;
t134 = Icges(6,4) * t174 - t176 * t210;
t264 = Icges(5,4) * t175;
t211 = Icges(5,1) * t173 + t264;
t135 = Icges(5,5) * t174 - t176 * t211;
t318 = t134 + t135;
t324 = t173 * t331 - t175 * t334;
t323 = t173 * t332 - t175 * t329;
t301 = -t156 * t332 + t188 * t329 - t255 * t326;
t322 = t156 * t318 - t188 * t319 + t255 * t327;
t237 = qJD(4) * t176;
t101 = (-Icges(6,3) * t173 - t262) * t237 + (Icges(6,6) * t176 + t174 * t205) * qJD(2);
t104 = (Icges(5,2) * t173 - t264) * t237 + (Icges(5,6) * t176 + t174 * t208) * qJD(2);
t321 = -t101 + t104;
t105 = (-Icges(6,1) * t175 - t263) * t237 + (Icges(6,4) * t176 + t174 * t210) * qJD(2);
t106 = (-Icges(5,1) * t175 + t265) * t237 + (Icges(5,5) * t176 + t174 * t211) * qJD(2);
t320 = t105 + t106;
t317 = ((Icges(4,5) - Icges(3,6)) * t176 + (Icges(4,4) - Icges(3,5)) * t174) * qJD(2);
t250 = t174 * t132;
t251 = t174 * t131;
t316 = -t250 - t251;
t252 = t174 * ((-Icges(6,4) * t175 - Icges(6,6) * t173) * t237 + (Icges(6,2) * t176 + t174 * t207) * qJD(2));
t253 = t174 * ((-Icges(5,5) * t175 + Icges(5,6) * t173) * t237 + (Icges(5,3) * t176 + t174 * t206) * qJD(2));
t315 = t253 + t252;
t314 = t330 * t171;
t168 = t171 ^ 2;
t169 = t172 ^ 2;
t294 = t168 + t169;
t292 = qJD(2) * t294;
t313 = m(3) * (rSges(3,1) * t174 + rSges(3,2) * t176) * t292;
t312 = -t294 + 0.1e1;
t311 = t174 * t326 - t176 * t323;
t310 = t174 * t325 - t176 * t324;
t309 = t173 * t318 - t175 * t319;
t240 = qJD(2) * t176;
t233 = t173 * t240;
t115 = qJD(4) * t188 + t172 * t233;
t232 = t175 * t240;
t116 = qJD(4) * t156 - t172 * t232;
t241 = qJD(2) * t174;
t234 = t172 * t241;
t67 = Icges(6,4) * t115 - Icges(6,2) * t234 + Icges(6,6) * t116;
t190 = t176 * t67 - t241 * t86;
t63 = Icges(6,5) * t115 - Icges(6,6) * t234 + Icges(6,3) * t116;
t71 = Icges(6,1) * t115 - Icges(6,4) * t234 + Icges(6,5) * t116;
t14 = t115 * t90 + t116 * t82 + t156 * t71 + t172 * t190 - t188 * t63;
t117 = qJD(4) * t157 + t171 * t233;
t118 = -qJD(4) * t236 + (qJD(4) * t172 + t171 * t240) * t175;
t235 = t171 * t241;
t68 = Icges(6,4) * t117 - Icges(6,2) * t235 - Icges(6,6) * t118;
t189 = t176 * t68 - t241 * t87;
t64 = Icges(6,5) * t117 - Icges(6,6) * t235 - Icges(6,3) * t118;
t72 = Icges(6,1) * t117 - Icges(6,4) * t235 - Icges(6,5) * t118;
t15 = t115 * t91 + t116 * t83 + t156 * t72 + t172 * t189 - t188 * t64;
t65 = Icges(5,5) * t115 - Icges(5,6) * t116 - Icges(5,3) * t234;
t192 = t176 * t65 - t241 * t84;
t69 = Icges(5,4) * t115 - Icges(5,2) * t116 - Icges(5,6) * t234;
t73 = Icges(5,1) * t115 - Icges(5,4) * t116 - Icges(5,5) * t234;
t16 = t115 * t92 - t116 * t88 + t156 * t73 + t172 * t192 + t188 * t69;
t66 = Icges(5,5) * t117 + Icges(5,6) * t118 - Icges(5,3) * t235;
t191 = t176 * t66 - t241 * t85;
t70 = Icges(5,4) * t117 + Icges(5,2) * t118 - Icges(5,6) * t235;
t74 = Icges(5,1) * t117 + Icges(5,4) * t118 - Icges(5,5) * t235;
t17 = t115 * t93 - t116 * t89 + t156 * t74 + t172 * t191 + t188 * t70;
t308 = ((t14 + t16 + t315) * t172 + (t15 + t17) * t171) * t176 + (t318 * t115 + t319 * t116 + t320 * t156 + t321 * t188) * t174 + (t322 * t176 + ((t316 + t301) * t172 + t314) * t174) * qJD(2);
t307 = rSges(6,1) + pkin(4);
t306 = ((t63 - t69) * t175 + (-t71 - t73) * t173 + (-t173 * t329 - t175 * t332) * qJD(4) + t326 * qJD(2)) * t176 + (qJD(2) * t323 + t65 + t67) * t174;
t305 = ((t64 - t70) * t175 + (-t72 - t74) * t173 + (-t173 * t334 - t175 * t331) * qJD(4) + t325 * qJD(2)) * t176 + (qJD(2) * t324 + t66 + t68) * t174;
t304 = rSges(6,3) + qJ(5);
t303 = t317 * t171;
t302 = t317 * t172;
t36 = -t157 * t83 + t158 * t91 + t257 * t87;
t38 = t157 * t89 + t158 * t93 + t257 * t85;
t300 = -t38 - t36;
t299 = -t176 * t309 - t316;
t293 = t171 * t310 + t172 * t311;
t268 = rSges(6,2) * t257 - t157 * t304 + t158 * t307;
t269 = rSges(6,2) * t255 + t156 * t307 - t188 * t304;
t291 = t171 * t269 - t172 * t268;
t288 = 2 * m(5);
t287 = 2 * m(6);
t286 = 0.2e1 * qJD(2);
t285 = m(4) / 0.2e1;
t284 = m(5) / 0.2e1;
t283 = m(6) / 0.2e1;
t280 = -t172 / 0.2e1;
t275 = -rSges(6,2) * t235 - qJD(5) * t157 + t117 * t307 - t118 * t304;
t276 = -rSges(6,2) * t234 - qJD(5) * t188 + t115 * t307 + t116 * t304;
t9 = (-t171 * t276 + t172 * t275) * t176 + t291 * t241;
t277 = t176 * t9;
t35 = -t157 * t82 + t158 * t90 + t257 * t86;
t272 = t172 * t35;
t37 = t157 * t88 + t158 * t92 + t257 * t84;
t271 = t172 * t37;
t164 = pkin(2) * t174 - qJ(3) * t176;
t247 = t294 * (-qJD(2) * t164 + qJD(3) * t174);
t178 = -pkin(6) * t241 * t294 + t247;
t22 = t171 * t275 + t172 * t276 + t178;
t270 = t176 * t22;
t258 = t171 * t174;
t256 = t172 * t174;
t223 = rSges(6,1) * t173 - rSges(6,3) * t175;
t248 = (pkin(4) * t241 - qJ(5) * t237) * t173 + (-qJ(5) * t241 + (-pkin(4) * qJD(4) + qJD(5)) * t176) * t175 + (-rSges(6,1) * t175 - rSges(6,3) * t173) * t237 + (rSges(6,2) * t176 + t174 * t223) * qJD(2);
t246 = t174 * rSges(6,2) + (-pkin(4) * t173 + qJ(5) * t175 - t223) * t176;
t220 = pkin(2) * t176 + qJ(3) * t174;
t245 = t294 * t220;
t152 = qJD(2) * t220 - qJD(3) * t176;
t222 = -rSges(4,2) * t176 + rSges(4,3) * t174;
t244 = -qJD(2) * t222 - t152;
t221 = rSges(4,2) * t174 + rSges(4,3) * t176;
t243 = -t164 + t221;
t239 = qJD(4) * t173;
t231 = t173 * t237;
t230 = t176 ^ 2 * t239;
t229 = -pkin(6) * t174 - t164;
t228 = t245 + (t171 * t257 + t172 * t255) * pkin(6);
t224 = rSges(5,1) * t173 + rSges(5,2) * t175;
t137 = rSges(5,3) * t174 - t176 * t224;
t227 = -t137 + t229;
t226 = -pkin(6) * t240 - t152;
t97 = rSges(5,1) * t156 + rSges(5,2) * t188 + rSges(5,3) * t255;
t99 = rSges(5,1) * t158 + rSges(5,2) * t157 + rSges(5,3) * t257;
t217 = t171 * t97 - t172 * t99;
t194 = t229 - t246;
t108 = (-rSges(5,1) * t175 + rSges(5,2) * t173) * t237 + (rSges(5,3) * t176 + t174 * t224) * qJD(2);
t193 = -t108 + t226;
t187 = t226 - t248;
t121 = t243 * t172;
t120 = t243 * t171;
t110 = t244 * t172;
t109 = t244 * t171;
t95 = t227 * t172;
t94 = t227 * t171;
t80 = t194 * t172;
t79 = t194 * t171;
t78 = rSges(5,1) * t117 + rSges(5,2) * t118 - rSges(5,3) * t235;
t76 = rSges(5,1) * t115 - rSges(5,2) * t116 - rSges(5,3) * t234;
t62 = -t137 * t255 + t174 * t97;
t61 = t137 * t257 - t174 * t99;
t60 = t193 * t172;
t59 = t193 * t171;
t54 = t222 * t294 + t245;
t53 = t131 * t257 + t133 * t157 + t135 * t158;
t52 = -t130 * t157 + t132 * t257 + t134 * t158;
t49 = t217 * t176;
t48 = t221 * t292 + t247;
t47 = t187 * t172;
t46 = t187 * t171;
t45 = t174 * t269 - t246 * t255;
t44 = -t174 * t268 + t246 * t257;
t39 = t171 * t99 + t172 * t97 + t228;
t30 = t291 * t176;
t29 = -t108 * t255 + t174 * t76 + (t137 * t256 + t176 * t97) * qJD(2);
t28 = t108 * t257 - t174 * t78 + (-t137 * t258 - t176 * t99) * qJD(2);
t27 = t171 * t78 + t172 * t76 + t178;
t26 = t171 * t268 + t172 * t269 + t228;
t25 = (-t171 * t76 + t172 * t78) * t176 + t217 * t241;
t24 = t276 * t174 - t248 * t255 + (t176 * t269 + t246 * t256) * qJD(2);
t23 = -t275 * t174 + t248 * t257 + (-t176 * t268 - t246 * t258) * qJD(2);
t21 = t117 * t93 + t118 * t89 + t157 * t70 + t158 * t74 + t171 * t191;
t20 = t117 * t92 + t118 * t88 + t157 * t69 + t158 * t73 + t171 * t192;
t19 = t117 * t91 - t118 * t83 - t157 * t64 + t158 * t72 + t171 * t189;
t18 = t117 * t90 - t118 * t82 - t157 * t63 + t158 * t71 + t171 * t190;
t8 = t171 * t20 - t172 * t21;
t7 = t171 * t18 - t172 * t19;
t6 = t16 * t171 - t17 * t172;
t5 = t14 * t171 - t15 * t172;
t4 = (t104 * t157 + t106 * t158 + t117 * t135 + t118 * t133) * t174 + (t20 * t172 + (t21 + t253) * t171) * t176 + (t53 * t176 + (-t271 + (-t38 - t251) * t171) * t174) * qJD(2);
t3 = (-t101 * t157 + t105 * t158 + t117 * t134 - t118 * t130) * t174 + (t18 * t172 + (t19 + t252) * t171) * t176 + (t52 * t176 + (-t272 + (-t36 - t250) * t171) * t174) * qJD(2);
t1 = [0; m(4) * t48 + m(5) * t27 + m(6) * t22 - t313; (t26 * t22 + t46 * t79 + t47 * t80) * t287 + (t39 * t27 + t59 * t94 + t60 * t95) * t288 + 0.2e1 * m(4) * (t109 * t120 + t110 * t121 + t54 * t48) + (-t169 * t303 - t7 - t8) * t172 + (t5 + t6 + t302 * t168 + (-t171 * t303 + t172 * t302) * t172) * t171 + 0.2e1 * t312 * (rSges(3,1) * t176 - rSges(3,2) * t174) * t313; (m(4) + m(5) + m(6)) * t241; m(6) * (t256 * t47 + t258 * t46 - t270) + m(5) * (-t176 * t27 + t256 * t60 + t258 * t59) + m(4) * (t109 * t258 + t110 * t256 - t176 * t48) + ((t174 * t26 + t255 * t80 + t257 * t79) * t283 + (t174 * t39 + t255 * t95 + t257 * t94) * t284 + (t120 * t257 + t121 * t255 + t174 * t54) * t285) * t286; -0.4e1 * (t285 + t284 + t283) * t312 * t174 * t240; m(5) * t25 + m(6) * t9; t4 * t280 + t3 * t280 + m(6) * (-t30 * t22 + t23 * t80 + t24 * t79 + t9 * t26 + t44 * t47 + t45 * t46) + m(5) * (t25 * t39 - t49 * t27 + t28 * t95 + t29 * t94 + t62 * t59 + t61 * t60) + ((t6 / 0.2e1 + t5 / 0.2e1) * t172 + (t7 / 0.2e1 + t8 / 0.2e1) * t171) * t176 + ((-(t300 * t172 + (t35 + t37) * t171) * t171 / 0.2e1 + (-t301 * t171 + t172 * t330) * t280) * t174 + (t311 * t171 - t310 * t172) * t176 / 0.2e1) * qJD(2) + t308 * t171 / 0.2e1 + (t306 * t171 - t305 * t172) * t174 / 0.2e1; m(5) * (-t176 * t25 + t256 * t28 + t258 * t29) + m(6) * (t23 * t256 + t24 * t258 - t277) + ((-t174 * t49 + t255 * t61 + t257 * t62) * t284 + (-t174 * t30 + t255 * t44 + t257 * t45) * t283) * t286; (t23 * t44 + t24 * t45 - t30 * t9) * t287 + (-t25 * t49 + t28 * t61 + t29 * t62) * t288 + (t4 + t3) * t257 + t308 * t255 + (t293 * t240 + (t300 * t171 - t271 - t272) * t235 + (t172 * t301 + t314) * t234) * t176 + (t315 * t174 + t299 * t240 + (-t53 - t52) * t235 - t322 * t234 + ((-t319 * t239 - t320 * t173 + (-qJD(4) * t318 - t321) * t175) * t174 + t306 * t172 + t305 * t171) * t176 + ((t309 * t174 - t293) * t174 + (t174 * t327 + t299) * t176) * qJD(2)) * t174; (-t175 * t241 - t231) * m(6); m(6) * (-t26 * t231 + t116 * t80 - t118 * t79 - t188 * t47 - t157 * t46 + (-t241 * t26 + t270) * t175); m(6) * (t230 + (t116 * t172 - t118 * t171) * t174 + (-t157 * t171 - t172 * t188 + 0.2e1 * t249) * t240); m(6) * (t30 * t231 + t116 * t44 - t118 * t45 - t188 * t23 - t157 * t24 + (t241 * t30 + t277) * t175); (-t188 * t116 + t157 * t118 + (-t174 * t232 - t230) * t175) * t287;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
