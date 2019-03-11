% Calculate time derivative of joint inertia matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:14
% EndTime: 2019-03-09 01:37:23
% DurationCPUTime: 5.84s
% Computational Cost: add. (9812->488), mult. (23977->713), div. (0->0), fcn. (26843->8), ass. (0->228)
t166 = cos(pkin(9));
t171 = cos(qJ(1));
t259 = sin(pkin(9));
t272 = sin(qJ(1));
t221 = t272 * t259;
t138 = t166 * t171 - t221;
t139 = t166 * t272 + t171 * t259;
t168 = sin(qJ(5));
t170 = cos(qJ(5));
t257 = Icges(6,4) * t170;
t202 = -Icges(6,2) * t168 + t257;
t86 = -Icges(6,6) * t139 - t138 * t202;
t258 = Icges(6,4) * t168;
t204 = Icges(6,1) * t170 - t258;
t88 = -Icges(6,5) * t139 - t138 * t204;
t205 = t168 * t86 - t170 * t88;
t189 = t205 * t139;
t85 = -Icges(6,6) * t138 + t139 * t202;
t87 = -Icges(6,5) * t138 + t139 * t204;
t206 = t168 * t85 - t170 * t87;
t290 = -t206 * t138 + t189;
t200 = Icges(6,5) * t170 - Icges(6,6) * t168;
t83 = -Icges(6,3) * t138 + t139 * t200;
t289 = t139 * t83;
t133 = t139 * qJD(1);
t240 = qJD(1) * t171;
t134 = -qJD(1) * t221 + t166 * t240;
t238 = qJD(5) * t168;
t190 = t134 * t170 - t139 * t238;
t188 = t190 * rSges(6,1) + t133 * rSges(6,3);
t237 = qJD(5) * t170;
t251 = t134 * t168;
t191 = t139 * t237 + t251;
t192 = t133 * t170 + t138 * t238;
t253 = t133 * t168;
t193 = t138 * t237 - t253;
t264 = t134 * rSges(6,3);
t247 = t139 * t170;
t248 = t139 * t168;
t89 = rSges(6,1) * t247 - rSges(6,2) * t248 - t138 * rSges(6,3);
t215 = rSges(6,1) * t170 - rSges(6,2) * t168;
t263 = t139 * rSges(6,3);
t90 = -t138 * t215 - t263;
t27 = t133 * t90 - t138 * (t192 * rSges(6,1) - t264) + t134 * t89 + t139 * t188 + (-t138 * t193 - t139 * t191) * rSges(6,2);
t280 = 2 * m(6);
t288 = t27 * t280;
t151 = rSges(6,1) * t168 + rSges(6,2) * t170;
t287 = t151 ^ 2 * t280;
t84 = -Icges(6,3) * t139 - t138 * t200;
t286 = t138 * t84 + t289 + t290;
t167 = sin(qJ(6));
t169 = cos(qJ(6));
t183 = qJD(6) * t138 - t190;
t235 = qJD(6) * t170;
t218 = -t139 * t235 + t133;
t59 = t167 * t183 + t169 * t218;
t60 = t167 * t218 - t169 * t183;
t39 = t60 * rSges(7,1) + t59 * rSges(7,2) + t191 * rSges(7,3);
t285 = t190 * pkin(5) + t191 * pkin(8) + t39;
t268 = pkin(5) * t170;
t219 = pkin(8) * t168 + t268;
t284 = t138 * t219;
t145 = t215 * qJD(5);
t282 = t145 * t151 * t280;
t273 = rSges(7,3) + pkin(8);
t281 = t168 * t273 + pkin(4) + t268;
t279 = 2 * m(7);
t278 = t133 / 0.2e1;
t277 = -t138 / 0.2e1;
t276 = -t139 / 0.2e1;
t275 = t139 / 0.2e1;
t274 = -t170 / 0.2e1;
t271 = m(6) * t145;
t270 = m(6) * t151;
t269 = pkin(5) * t168;
t245 = t167 * t170;
t102 = -t138 * t169 - t139 * t245;
t244 = t169 * t170;
t103 = -t138 * t167 + t139 * t244;
t71 = Icges(7,4) * t103 + Icges(7,2) * t102 + Icges(7,6) * t248;
t73 = Icges(7,1) * t103 + Icges(7,4) * t102 + Icges(7,5) * t248;
t208 = -t167 * t71 + t169 * t73;
t33 = Icges(7,5) * t60 + Icges(7,6) * t59 + Icges(7,3) * t191;
t35 = Icges(7,4) * t60 + Icges(7,2) * t59 + Icges(7,6) * t191;
t37 = Icges(7,1) * t60 + Icges(7,4) * t59 + Icges(7,5) * t191;
t69 = Icges(7,5) * t103 + Icges(7,6) * t102 + Icges(7,3) * t248;
t9 = (qJD(5) * t208 - t33) * t170 + (qJD(5) * t69 - t167 * t35 + t169 * t37 + (-t167 * t73 - t169 * t71) * qJD(6)) * t168;
t267 = t9 * t138;
t266 = -pkin(1) - qJ(3);
t104 = t138 * t245 - t139 * t169;
t105 = -t138 * t244 - t139 * t167;
t249 = t138 * t168;
t72 = Icges(7,4) * t105 + Icges(7,2) * t104 - Icges(7,6) * t249;
t74 = Icges(7,1) * t105 + Icges(7,4) * t104 - Icges(7,5) * t249;
t207 = -t167 * t72 + t169 * t74;
t184 = -qJD(6) * t139 + t192;
t217 = t138 * t235 - t134;
t57 = -t167 * t184 + t169 * t217;
t58 = t167 * t217 + t169 * t184;
t32 = Icges(7,5) * t58 + Icges(7,6) * t57 - Icges(7,3) * t193;
t34 = Icges(7,4) * t58 + Icges(7,2) * t57 - Icges(7,6) * t193;
t36 = Icges(7,1) * t58 + Icges(7,4) * t57 - Icges(7,5) * t193;
t70 = Icges(7,5) * t105 + Icges(7,6) * t104 - Icges(7,3) * t249;
t10 = (qJD(5) * t207 - t32) * t170 + (qJD(5) * t70 - t167 * t34 + t169 * t36 + (-t167 * t74 - t169 * t72) * qJD(6)) * t168;
t265 = t10 * t139;
t75 = t103 * rSges(7,1) + t102 * rSges(7,2) + rSges(7,3) * t248;
t262 = pkin(5) * t247 + pkin(8) * t248 + t75;
t214 = -t105 * rSges(7,1) - t104 * rSges(7,2);
t76 = -rSges(7,3) * t249 - t214;
t261 = -t76 + t284;
t213 = rSges(7,1) * t169 - rSges(7,2) * t167;
t236 = qJD(6) * t168;
t97 = (-rSges(7,1) * t167 - rSges(7,2) * t169) * t236 + (rSges(7,3) * t168 + t170 * t213) * qJD(5);
t260 = t219 * qJD(5) + t97;
t256 = Icges(7,4) * t167;
t255 = Icges(7,4) * t169;
t203 = Icges(7,1) * t169 - t256;
t120 = -Icges(7,5) * t170 + t168 * t203;
t254 = t120 * t169;
t121 = -rSges(7,3) * t170 + t168 * t213;
t152 = -pkin(8) * t170 + t269;
t243 = t121 + t152;
t242 = qJ(2) * t240 + qJD(2) * t272;
t160 = t272 * qJ(2);
t241 = t171 * pkin(1) + t160;
t239 = qJD(5) * t138;
t234 = t272 * pkin(1);
t163 = t272 * pkin(3);
t29 = t168 * t208 - t170 * t69;
t199 = Icges(7,5) * t169 - Icges(7,6) * t167;
t118 = -Icges(7,3) * t170 + t168 * t199;
t201 = -Icges(7,2) * t167 + t255;
t119 = -Icges(7,6) * t170 + t168 * t201;
t49 = t102 * t119 + t103 * t120 + t118 * t248;
t233 = t29 / 0.2e1 + t49 / 0.2e1;
t30 = t168 * t207 - t170 * t70;
t50 = t104 * t119 + t105 * t120 - t118 * t249;
t232 = t50 / 0.2e1 + t30 / 0.2e1;
t231 = t272 * rSges(4,1);
t230 = t272 * rSges(3,3);
t229 = qJD(3) * t171 + t242;
t228 = t171 * qJ(3) + t241;
t226 = t119 * t237;
t222 = t163 + t228;
t220 = t58 * rSges(7,1) + t57 * rSges(7,2);
t159 = qJD(2) * t171;
t216 = -qJD(3) * t272 + t159;
t23 = t102 * t71 + t103 * t73 + t248 * t69;
t24 = t102 * t72 + t103 * t74 + t248 * t70;
t212 = -t138 * t24 + t139 * t23;
t25 = t104 * t71 + t105 * t73 - t249 * t69;
t26 = t104 * t72 + t105 * t74 - t249 * t70;
t211 = -t138 * t26 + t139 * t25;
t210 = t30 * t138 - t29 * t139;
t209 = t138 * t75 + t139 * t76;
t198 = pkin(4) + t215;
t197 = t168 * t32 + t237 * t70;
t196 = t168 * t33 + t237 * t69;
t195 = t138 * t272 - t139 * t171;
t194 = t138 * t171 + t139 * t272;
t187 = -qJ(3) * t272 - t234;
t186 = t139 * pkin(4) - pkin(7) * t138 + t222;
t162 = t171 * qJ(2);
t182 = t171 * pkin(3) + t162 + t187;
t181 = rSges(3,2) * t272 + t171 * rSges(3,3) - t234;
t94 = (-Icges(7,5) * t167 - Icges(7,6) * t169) * t236 + (Icges(7,3) * t168 + t170 * t199) * qJD(5);
t96 = (-Icges(7,1) * t167 - t255) * t236 + (Icges(7,5) * t168 + t170 * t203) * qJD(5);
t180 = t118 * t238 - t170 * t94 + t237 * t254 + (-t119 * t236 + t168 * t96) * t169;
t179 = t139 * pkin(7) + t182;
t176 = t171 * rSges(4,1) - rSges(4,3) * t272 + t187;
t175 = pkin(3) * t240 + qJD(1) * t187 + t229;
t174 = t134 * pkin(4) + t133 * pkin(7) + t175;
t173 = (t171 * t266 - t160 - t163) * qJD(1) + t216;
t172 = t134 * pkin(7) + t173;
t125 = -t171 * rSges(3,2) + t230 + t241;
t124 = t162 + t181;
t113 = t171 * rSges(4,3) + t228 + t231;
t112 = t162 + t176;
t109 = t159 + (-t230 - t160 + (rSges(3,2) - pkin(1)) * t171) * qJD(1);
t108 = qJD(1) * t181 + t242;
t99 = (-t231 - t160 + (-rSges(4,3) + t266) * t171) * qJD(1) + t216;
t98 = qJD(1) * t176 + t229;
t95 = (-Icges(7,2) * t169 - t256) * t236 + (Icges(7,6) * t168 + t170 * t201) * qJD(5);
t92 = rSges(5,1) * t139 + rSges(5,2) * t138 + t222;
t91 = t138 * rSges(5,1) - t139 * rSges(5,2) + t182;
t82 = t243 * t138;
t81 = t243 * t139;
t80 = -t133 * rSges(5,1) - t134 * rSges(5,2) + t173;
t79 = t134 * rSges(5,1) - t133 * rSges(5,2) + t175;
t78 = -t118 * t170 + (-t119 * t167 + t254) * t168;
t77 = t78 * t238;
t68 = t186 + t89;
t67 = t138 * t198 + t179 + t263;
t62 = Icges(6,5) * t190 - Icges(6,6) * t191 + Icges(6,3) * t133;
t61 = Icges(6,5) * t192 + Icges(6,6) * t193 - Icges(6,3) * t134;
t54 = -t121 * t249 + t170 * t76;
t53 = -t121 * t248 - t170 * t75;
t52 = -t134 * t243 - t139 * t260;
t51 = -t133 * t243 + t138 * t260;
t48 = -t133 * t198 - t151 * t239 + t172 + t264;
t47 = -rSges(6,2) * t191 + t174 + t188;
t46 = t186 + t262;
t45 = t281 * t138 + t179 + t214;
t40 = t209 * t168;
t38 = -rSges(7,3) * t193 + t220;
t31 = (-t226 + (-qJD(6) * t120 - t95) * t168) * t167 + t180;
t28 = t138 * t261 + t139 * t262;
t22 = (t170 * t273 - t269) * t239 - t281 * t133 + t172 - t220;
t21 = t174 + t285;
t20 = (-qJD(5) * t121 * t139 - t39) * t170 + (qJD(5) * t75 - t121 * t134 - t139 * t97) * t168;
t19 = (-t121 * t239 + t38) * t170 + (-qJD(5) * t76 + t121 * t133 - t138 * t97) * t168;
t18 = t102 * t95 + t103 * t96 + t118 * t191 + t119 * t59 + t120 * t60 + t248 * t94;
t17 = t104 * t95 + t105 * t96 - t118 * t193 + t119 * t57 + t120 * t58 - t249 * t94;
t16 = -t138 * t25 - t139 * t26;
t15 = -t138 * t23 - t139 * t24;
t14 = t168 * t211 - t170 * t50;
t13 = t168 * t212 - t170 * t49;
t12 = t209 * t237 + (-t133 * t75 + t134 * t76 + t138 * t39 + t139 * t38) * t168;
t11 = t285 * t139 + t262 * t134 + (-t152 * t239 - t38) * t138 + (-t261 - t284) * t133;
t8 = t102 * t34 + t103 * t36 + t139 * t197 + t251 * t70 + t59 * t72 + t60 * t74;
t7 = t102 * t35 + t103 * t37 + t139 * t196 + t251 * t69 + t59 * t71 + t60 * t73;
t6 = t104 * t34 + t105 * t36 - t138 * t197 + t253 * t70 + t57 * t72 + t58 * t74;
t5 = t104 * t35 + t105 * t37 - t138 * t196 + t253 * t69 + t57 * t71 + t58 * t73;
t4 = t133 * t23 - t134 * t24 - t138 * t7 - t139 * t8;
t3 = t133 * t25 - t134 * t26 - t138 * t5 - t139 * t6;
t2 = (qJD(5) * t212 - t18) * t170 + (qJD(5) * t49 + t133 * t24 + t134 * t23 - t138 * t8 + t139 * t7) * t168;
t1 = (qJD(5) * t211 - t17) * t170 + (qJD(5) * t50 + t133 * t26 + t134 * t25 - t138 * t6 + t139 * t5) * t168;
t41 = [(t21 * t46 + t22 * t45) * t279 + (t47 * t68 + t48 * t67) * t280 + 0.2e1 * m(4) * (t112 * t99 + t113 * t98) + 0.2e1 * m(5) * (t79 * t92 + t80 * t91) + 0.2e1 * m(3) * (t108 * t125 + t109 * t124) + t180 + (-Icges(6,2) * t170 + t204 - t258) * t238 + (Icges(6,1) * t168 + t202 + t257) * t237 + (-t120 * t236 - t168 * t95 - t226) * t167; m(7) * (t272 * t22 - t171 * t21 + (t171 * t45 + t272 * t46) * qJD(1)) + m(6) * (t272 * t48 - t171 * t47 + (t171 * t67 + t272 * t68) * qJD(1)) + m(4) * (t272 * t99 - t171 * t98 + (t112 * t171 + t113 * t272) * qJD(1)) + m(5) * (t272 * t80 - t171 * t79 + (t171 * t91 + t272 * t92) * qJD(1)) + m(3) * (t272 * t109 - t171 * t108 + (t124 * t171 + t125 * t272) * qJD(1)); 0; m(7) * (t272 * t21 + t171 * t22 + (t171 * t46 - t272 * t45) * qJD(1)) + m(6) * (t272 * t47 + t171 * t48 + (t171 * t68 - t272 * t67) * qJD(1)) + m(4) * (t272 * t98 + t171 * t99 + (-t112 * t272 + t113 * t171) * qJD(1)) + m(5) * (t272 * t79 + t171 * t80 + (t171 * t92 - t272 * t91) * qJD(1)); 0; 0; 0; 0; 0; 0; -t267 / 0.2e1 - t265 / 0.2e1 + (t138 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1) * t200 * qJD(5) + (-t168 * t88 / 0.2e1 + t86 * t274 - t232) * t134 + (t168 * t87 / 0.2e1 + t170 * t85 / 0.2e1 + t233) * t133 + m(7) * (t21 * t82 - t22 * t81 + t52 * t45 + t51 * t46) + (t138 * t68 - t139 * t67) * t271 + (-t133 * t68 - t134 * t67 + t138 * t47 - t139 * t48) * t270 + (-t133 * t138 + t134 * t139) * (Icges(6,5) * t168 + Icges(6,6) * t170) + (-qJD(5) * t206 + t168 * (Icges(6,1) * t190 - Icges(6,4) * t191 + Icges(6,5) * t133) + t170 * (Icges(6,4) * t190 - Icges(6,2) * t191 + Icges(6,6) * t133) + t18) * t277 + (-qJD(5) * t205 + t168 * (Icges(6,1) * t192 + Icges(6,4) * t193 - Icges(6,5) * t134) + t170 * (Icges(6,4) * t192 + Icges(6,2) * t193 - Icges(6,6) * t134) + t17) * t276; m(7) * (t52 * t272 - t51 * t171 + (-t171 * t81 + t272 * t82) * qJD(1)) - t194 * t271 + (qJD(1) * t195 + t133 * t171 - t134 * t272) * t270; m(7) * (t51 * t272 + t52 * t171 + (t171 * t82 + t272 * t81) * qJD(1)) + t195 * t271 + (qJD(1) * t194 - t133 * t272 - t134 * t171) * t270; m(6) * t27 + m(7) * t11; (t28 * t11 + t51 * t82 - t52 * t81) * t279 + (t89 * t288 - t3 + (-t139 * t61 + t282) * t139) * t139 + (-t90 * t288 - t4 + (-t138 * t62 + t282) * t138 + (-t138 * t61 - t139 * t62) * t139) * t138 + (-t16 + (-0.3e1 * t139 * t84 + t287) * t139 + (-t286 - t289 + t290) * t138) * t134 + (t15 + (0.3e1 * t138 * t83 - t287) * t138 + (-t189 + t286 + (t206 + t84) * t138) * t139) * t133; t77 + m(7) * (t19 * t45 + t20 * t46 + t21 * t53 + t22 * t54) + (-t31 + (-t138 * t232 + t139 * t233) * qJD(5)) * t170 + ((t9 / 0.2e1 + t18 / 0.2e1) * t139 + (-t17 / 0.2e1 - t10 / 0.2e1) * t138 + t233 * t134 + t232 * t133) * t168; m(7) * (t19 * t272 - t20 * t171 + (t171 * t54 + t272 * t53) * qJD(1)); m(7) * (t20 * t272 + t19 * t171 + (t171 * t53 - t272 * t54) * qJD(1)); m(7) * t12; t13 * t278 + t2 * t277 - t134 * t14 / 0.2e1 + t1 * t276 + (t29 * t133 - t30 * t134 - t265 - t267) * t274 + m(7) * (t40 * t11 + t12 * t28 - t19 * t81 + t20 * t82 + t53 * t51 + t54 * t52) + (t15 * t275 + t16 * t277) * t237 + (t134 * t15 / 0.2e1 + t4 * t275 + t16 * t278 + t3 * t277 + qJD(5) * (-t138 * t29 - t139 * t30) / 0.2e1) * t168; (t12 * t40 + t19 * t54 + t20 * t53) * t279 + (t31 * t170 - t77 + (t139 * t13 - t138 * t14 + t170 * t210) * qJD(5)) * t170 + (t134 * t13 + t139 * t2 + t133 * t14 - t138 * t1 - t170 * (-t10 * t138 + t30 * t133 + t29 * t134 + t9 * t139) + (-t168 * t210 - t78 * t170) * qJD(5)) * t168;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t41(1) t41(2) t41(4) t41(7) t41(11) t41(16); t41(2) t41(3) t41(5) t41(8) t41(12) t41(17); t41(4) t41(5) t41(6) t41(9) t41(13) t41(18); t41(7) t41(8) t41(9) t41(10) t41(14) t41(19); t41(11) t41(12) t41(13) t41(14) t41(15) t41(20); t41(16) t41(17) t41(18) t41(19) t41(20) t41(21);];
Mq  = res;
