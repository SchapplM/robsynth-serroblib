% Calculate time derivative of joint inertia matrix for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:26
% DurationCPUTime: 4.38s
% Computational Cost: add. (11570->294), mult. (11705->475), div. (0->0), fcn. (10719->10), ass. (0->190)
t177 = qJ(2) + qJ(3);
t169 = pkin(9) + t177;
t166 = sin(t169);
t167 = cos(t169);
t170 = sin(t177);
t171 = cos(t177);
t321 = Icges(4,5) * t170 + Icges(5,5) * t166 + Icges(4,6) * t171 + Icges(5,6) * t167;
t176 = qJD(2) + qJD(3);
t320 = t321 * t176;
t178 = sin(pkin(8));
t174 = t178 ^ 2;
t179 = cos(pkin(8));
t175 = t179 ^ 2;
t309 = t174 + t175;
t308 = t176 * t309;
t319 = t320 * t178;
t318 = t320 * t179;
t317 = t321 * t308;
t162 = rSges(4,1) * t170 + rSges(4,2) * t171;
t310 = t162 * t309;
t182 = cos(qJ(5));
t260 = t179 * t182;
t180 = sin(qJ(5));
t263 = t178 * t180;
t149 = -t167 * t263 - t260;
t261 = t179 * t180;
t262 = t178 * t182;
t150 = t167 * t262 - t261;
t266 = t176 * t167;
t264 = t178 * t176;
t249 = t167 * t264;
t267 = t176 * t166;
t251 = t180 * t267;
t90 = -qJD(5) * t150 + t178 * t251;
t250 = t182 * t267;
t91 = qJD(5) * t149 - t178 * t250;
t54 = Icges(6,5) * t91 + Icges(6,6) * t90 + Icges(6,3) * t249;
t269 = t166 * t178;
t75 = Icges(6,5) * t150 + Icges(6,6) * t149 + Icges(6,3) * t269;
t206 = t166 * t54 + t266 * t75;
t56 = Icges(6,4) * t91 + Icges(6,2) * t90 + Icges(6,6) * t249;
t58 = Icges(6,1) * t91 + Icges(6,4) * t90 + Icges(6,5) * t249;
t77 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t269;
t79 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t269;
t13 = t149 * t56 + t150 * t58 + t178 * t206 + t77 * t90 + t79 * t91;
t265 = t176 * t179;
t248 = t167 * t265;
t152 = t167 * t260 + t263;
t92 = -qJD(5) * t152 + t179 * t251;
t151 = -t167 * t261 + t262;
t93 = qJD(5) * t151 - t179 * t250;
t55 = Icges(6,5) * t93 + Icges(6,6) * t92 + Icges(6,3) * t248;
t268 = t166 * t179;
t76 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t268;
t205 = t166 * t55 + t266 * t76;
t57 = Icges(6,4) * t93 + Icges(6,2) * t92 + Icges(6,6) * t248;
t59 = Icges(6,1) * t93 + Icges(6,4) * t92 + Icges(6,5) * t248;
t78 = Icges(6,4) * t152 + Icges(6,2) * t151 + Icges(6,6) * t268;
t80 = Icges(6,1) * t152 + Icges(6,4) * t151 + Icges(6,5) * t268;
t14 = t149 * t57 + t150 * t59 + t178 * t205 + t78 * t90 + t80 * t91;
t9 = -t13 * t179 + t14 * t178;
t306 = -t9 + (-t318 * t179 + t317) * t178 + t319 * t175;
t235 = pkin(4) * t167 + pkin(7) * t166;
t231 = rSges(6,1) * t182 - rSges(6,2) * t180;
t68 = t231 * t266 + (rSges(6,3) * t176 + (-rSges(6,1) * t180 - rSges(6,2) * t182) * qJD(5)) * t166;
t304 = -t235 * t176 - t68;
t159 = rSges(5,1) * t166 + rSges(5,2) * t167;
t301 = t159 * t308;
t119 = -rSges(6,3) * t167 + t166 * t231;
t160 = pkin(4) * t166 - pkin(7) * t167;
t300 = -t119 - t160;
t181 = sin(qJ(2));
t183 = cos(qJ(2));
t299 = qJD(2) * (rSges(3,1) * t181 + rSges(3,2) * t183);
t60 = rSges(6,1) * t91 + rSges(6,2) * t90 + rSges(6,3) * t249;
t61 = rSges(6,1) * t93 + rSges(6,2) * t92 + rSges(6,3) * t248;
t297 = -t160 * t308 + t178 * t60 + t179 * t61;
t293 = 2 * m(4);
t292 = 2 * m(5);
t291 = 2 * m(6);
t290 = pkin(2) * t181;
t289 = pkin(3) * t170;
t288 = pkin(3) * t176;
t285 = t309 * pkin(3) * t171;
t284 = t309 * t170 * t288;
t283 = pkin(2) * qJD(2);
t214 = Icges(6,5) * t182 - Icges(6,6) * t180;
t282 = t167 * (t214 * t266 + (Icges(6,3) * t176 + (-Icges(6,5) * t180 - Icges(6,6) * t182) * qJD(5)) * t166);
t33 = t151 * t77 + t152 * t79 + t268 * t75;
t281 = t178 * t33;
t32 = t149 * t78 + t150 * t80 + t269 * t76;
t280 = t179 * t32;
t273 = Icges(6,4) * t180;
t272 = Icges(6,4) * t182;
t106 = -Icges(6,3) * t167 + t166 * t214;
t271 = t106 * t167;
t71 = t176 * t310;
t259 = t309 * t183 * pkin(2);
t233 = rSges(4,1) * t171 - rSges(4,2) * t170;
t74 = t309 * t233;
t15 = t151 * t56 + t152 * t58 + t179 * t206 + t77 * t92 + t79 * t93;
t16 = t151 * t57 + t152 * t59 + t179 * t205 + t78 * t92 + t80 * t93;
t10 = -t15 * t179 + t16 * t178;
t255 = (t10 - t318 * t174 + (t319 * t178 - t317) * t179) * t178;
t254 = t171 * t288;
t252 = t183 * t283;
t246 = -t162 - t290;
t245 = -t159 - t289;
t232 = rSges(5,1) * t167 - rSges(5,2) * t166;
t39 = t232 * t309 + t285;
t142 = t232 * t176;
t240 = -t142 - t254;
t239 = -t289 + t300;
t148 = t233 * t176;
t238 = -t148 - t252;
t236 = -t289 - t290;
t228 = -t180 * t77 + t182 * t79;
t36 = t166 * t228 - t167 * t75;
t227 = -t180 * t78 + t182 * t80;
t37 = t166 * t227 - t167 * t76;
t230 = t36 * t178 + t37 * t179;
t81 = rSges(6,1) * t150 + rSges(6,2) * t149 + rSges(6,3) * t269;
t82 = rSges(6,1) * t152 + rSges(6,2) * t151 + rSges(6,3) * t268;
t229 = -t178 * t82 + t179 * t81;
t219 = Icges(6,1) * t182 - t273;
t215 = -Icges(6,2) * t180 + t272;
t107 = -Icges(6,6) * t167 + t166 * t215;
t108 = -Icges(6,5) * t167 + t166 * t219;
t213 = t107 * t180 - t108 * t182;
t209 = -t254 + t304;
t208 = t309 * t181 * t283;
t30 = t178 * t81 + t179 * t82 + t235 * t309 + t285;
t207 = -t159 + t236;
t11 = (t176 * t228 - t54) * t167 + (t176 * t75 - t180 * t56 + t182 * t58 + (-t180 * t79 - t182 * t77) * qJD(5)) * t166;
t12 = (t176 * t227 - t55) * t167 + (t176 * t76 - t180 * t57 + t182 * t59 + (-t180 * t80 - t182 * t78) * qJD(5)) * t166;
t31 = t149 * t77 + t150 * t79 + t269 * t75;
t41 = t106 * t269 + t107 * t149 + t108 * t150;
t66 = t215 * t266 + (Icges(6,6) * t176 + (-Icges(6,2) * t182 - t273) * qJD(5)) * t166;
t67 = t219 * t266 + (Icges(6,5) * t176 + (-Icges(6,1) * t180 - t272) * qJD(5)) * t166;
t3 = (-t107 * t90 - t108 * t91 - t149 * t66 - t150 * t67 + (t280 + (t31 - t271) * t178) * t176) * t167 + (t14 * t179 + t41 * t176 + (t13 - t282) * t178) * t166;
t34 = t151 * t78 + t152 * t80 + t268 * t76;
t42 = t106 * t268 + t107 * t151 + t108 * t152;
t4 = (-t107 * t92 - t108 * t93 - t151 * t66 - t152 * t67 + (t281 + (t34 - t271) * t179) * t176) * t167 + (t15 * t178 + t42 * t176 + (t16 - t282) * t179) * t166;
t203 = t178 * t4 / 0.2e1 - t179 * t3 / 0.2e1 + (t178 * t37 - t179 * t36) * t267 / 0.2e1 - t167 * (-t11 * t179 + t12 * t178) / 0.2e1 + t9 * t269 / 0.2e1 + t10 * t268 / 0.2e1 + (t178 * (t178 * t32 - t179 * t31) + t179 * (t178 * t34 - t179 * t33)) * t266 / 0.2e1;
t200 = t236 + t300;
t192 = -t252 - t254;
t189 = qJD(2) * (-Icges(3,5) * t181 - Icges(3,6) * t183);
t188 = t306 * t179 + t255;
t187 = -t208 - t284;
t186 = -t142 + t192;
t185 = t192 + t304;
t154 = t179 * t189;
t153 = t178 * t189;
t137 = t246 * t179;
t136 = t246 * t178;
t110 = t245 * t179;
t109 = t245 * t178;
t99 = t238 * t179;
t98 = t238 * t178;
t97 = t207 * t179;
t96 = t207 * t178;
t89 = t240 * t179;
t88 = t240 * t178;
t85 = t309 * t299;
t84 = t186 * t179;
t83 = t186 * t178;
t70 = t239 * t179;
t69 = t239 * t178;
t64 = t200 * t179;
t63 = t200 * t178;
t62 = -t208 - t71;
t51 = -t119 * t268 - t167 * t82;
t50 = t119 * t269 + t167 * t81;
t49 = t209 * t179;
t48 = t209 * t178;
t47 = t74 + t259;
t46 = t185 * t179;
t45 = t185 * t178;
t44 = -t166 * t213 - t271;
t43 = t229 * t166;
t40 = -t284 - t301;
t38 = t187 - t301;
t35 = t39 + t259;
t29 = (-t119 * t265 - t61) * t167 + (t176 * t82 - t179 * t68) * t166;
t28 = (t119 * t264 + t60) * t167 + (-t176 * t81 + t178 * t68) * t166;
t27 = t30 + t259;
t26 = -t284 + t297;
t25 = t229 * t266 + (-t178 * t61 + t179 * t60) * t166;
t23 = t187 + t297;
t1 = [0; -m(3) * t85 + m(4) * t62 + m(5) * t38 + m(6) * t23; (t23 * t27 + t45 * t63 + t46 * t64) * t291 + (t35 * t38 + t83 * t96 + t84 * t97) * t292 + (t136 * t98 + t137 * t99 + t47 * t62) * t293 + t178 * t174 * t154 + t255 + 0.2e1 * m(3) * (t299 - t85) * t309 * (rSges(3,1) * t183 - rSges(3,2) * t181) + (-t175 * t153 + (-t178 * t153 + t179 * t154) * t178 + t306) * t179; -m(4) * t71 + m(5) * t40 + m(6) * t26; m(6) * (t23 * t30 + t26 * t27 + t45 * t69 + t46 * t70 + t48 * t63 + t49 * t64) + m(5) * (t109 * t83 + t110 * t84 + t40 * t35 + t39 * t38 + t88 * t96 + t89 * t97) + m(4) * (-t71 * t47 + t74 * t62 + (-t178 * t98 - t179 * t99) * t162 + (-t136 * t178 - t137 * t179) * t148) + t188; (t26 * t30 + t48 * t69 + t49 * t70) * t291 + (t109 * t88 + t110 * t89 + t39 * t40) * t292 + (t148 * t310 - t74 * t71) * t293 + t188; 0; m(6) * (t178 * t46 - t179 * t45) + m(5) * (t178 * t84 - t179 * t83); m(6) * (t178 * t49 - t179 * t48) + m(5) * (t178 * t89 - t179 * t88); 0; m(6) * t25; m(6) * (t23 * t43 + t25 * t27 + t28 * t64 + t29 * t63 + t45 * t51 + t46 * t50) + t203; m(6) * (t25 * t30 + t26 * t43 + t28 * t70 + t29 * t69 + t48 * t51 + t49 * t50) + t203; m(6) * (t178 * t28 - t179 * t29); (t25 * t43 + t28 * t50 + t29 * t51) * t291 + (-t167 * t42 + (t179 * t34 + t281) * t166) * t248 + t4 * t268 + (-t167 * t41 + (t178 * t31 + t280) * t166) * t249 + t3 * t269 + (t166 * t230 - t167 * t44) * t267 - t167 * ((t282 + (t167 * t213 + t230) * t176) * t167 + (t12 * t179 + t11 * t178 - (t106 * t176 - t180 * t66 + t182 * t67 + (-t107 * t182 - t108 * t180) * qJD(5)) * t167 + t44 * t176) * t166);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
