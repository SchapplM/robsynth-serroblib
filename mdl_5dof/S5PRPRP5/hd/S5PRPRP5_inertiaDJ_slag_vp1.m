% Calculate time derivative of joint inertia matrix for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:37:56
% DurationCPUTime: 9.51s
% Computational Cost: add. (9419->466), mult. (14854->717), div. (0->0), fcn. (14080->8), ass. (0->235)
t184 = pkin(8) + qJ(4);
t180 = sin(t184);
t181 = cos(t184);
t187 = sin(pkin(7));
t189 = cos(pkin(7));
t192 = cos(qJ(2));
t268 = t189 * t192;
t155 = t180 * t268 - t181 * t187;
t156 = t180 * t187 + t181 * t268;
t191 = sin(qJ(2));
t269 = t189 * t191;
t85 = Icges(6,5) * t156 + Icges(6,6) * t269 + Icges(6,3) * t155;
t91 = Icges(5,4) * t156 - Icges(5,2) * t155 + Icges(5,6) * t269;
t346 = t85 - t91;
t87 = Icges(5,5) * t156 - Icges(5,6) * t155 + Icges(5,3) * t269;
t89 = Icges(6,4) * t156 + Icges(6,2) * t269 + Icges(6,6) * t155;
t347 = t87 + t89;
t93 = Icges(6,1) * t156 + Icges(6,4) * t269 + Icges(6,5) * t155;
t95 = Icges(5,1) * t156 - Icges(5,4) * t155 + Icges(5,5) * t269;
t345 = t93 + t95;
t271 = t187 * t192;
t153 = t180 * t271 + t181 * t189;
t154 = -t180 * t189 + t181 * t271;
t272 = t187 * t191;
t84 = Icges(6,5) * t154 + Icges(6,6) * t272 + Icges(6,3) * t153;
t90 = Icges(5,4) * t154 - Icges(5,2) * t153 + Icges(5,6) * t272;
t343 = t84 - t90;
t92 = Icges(6,1) * t154 + Icges(6,4) * t272 + Icges(6,5) * t153;
t94 = Icges(5,1) * t154 - Icges(5,4) * t153 + Icges(5,5) * t272;
t341 = t92 + t94;
t344 = t153 * t346 + t154 * t345 + t272 * t347;
t86 = Icges(5,5) * t154 - Icges(5,6) * t153 + Icges(5,3) * t272;
t88 = Icges(6,4) * t154 + Icges(6,2) * t272 + Icges(6,6) * t153;
t342 = t88 + t86;
t277 = Icges(6,5) * t181;
t218 = Icges(6,3) * t180 + t277;
t140 = -Icges(6,6) * t192 + t191 * t218;
t279 = Icges(5,4) * t181;
t221 = -Icges(5,2) * t180 + t279;
t143 = -Icges(5,6) * t192 + t191 * t221;
t331 = t140 - t143;
t219 = Icges(5,5) * t181 - Icges(5,6) * t180;
t141 = -Icges(5,3) * t192 + t191 * t219;
t220 = Icges(6,4) * t181 + Icges(6,6) * t180;
t142 = -Icges(6,2) * t192 + t191 * t220;
t340 = t141 + t142;
t278 = Icges(6,5) * t180;
t223 = Icges(6,1) * t181 + t278;
t144 = -Icges(6,4) * t192 + t191 * t223;
t280 = Icges(5,4) * t180;
t224 = Icges(5,1) * t181 - t280;
t145 = -Icges(5,5) * t192 + t191 * t224;
t339 = t144 + t145;
t336 = t180 * t346 + t181 * t345;
t335 = t180 * t343 + t181 * t341;
t316 = t153 * t343 + t154 * t341 + t272 * t342;
t334 = t153 * t331 + t154 * t339 + t272 * t340;
t253 = qJD(4) * t191;
t100 = (Icges(6,3) * t181 - t278) * t253 + (Icges(6,6) * t191 + t192 * t218) * qJD(2);
t103 = (-Icges(5,2) * t181 - t280) * t253 + (Icges(5,6) * t191 + t192 * t221) * qJD(2);
t333 = t100 - t103;
t104 = (-Icges(6,1) * t180 + t277) * t253 + (Icges(6,4) * t191 + t192 * t223) * qJD(2);
t105 = (-Icges(5,1) * t180 - t279) * t253 + (Icges(5,5) * t191 + t192 * t224) * qJD(2);
t332 = t104 + t105;
t266 = t192 * ((-Icges(6,4) * t180 + Icges(6,6) * t181) * t253 + (Icges(6,2) * t191 + t192 * t220) * qJD(2));
t267 = t192 * ((-Icges(5,5) * t180 - Icges(5,6) * t181) * t253 + (Icges(5,3) * t191 + t192 * t219) * qJD(2));
t329 = -t267 - t266;
t275 = t142 * t192;
t276 = t141 * t192;
t328 = -t275 - t276;
t327 = t344 * t189;
t182 = t187 ^ 2;
t183 = t189 ^ 2;
t311 = t182 + t183;
t326 = 0.1e1 - t311;
t308 = qJD(2) * t311;
t325 = m(3) * (rSges(3,1) * t191 + rSges(3,2) * t192) * t308;
t324 = t191 * t335 - t192 * t342;
t323 = t191 * t336 - t192 * t347;
t322 = -t180 * t331 - t181 * t339;
t254 = qJD(4) * t181;
t244 = t192 * t254;
t257 = qJD(2) * t191;
t250 = t187 * t257;
t125 = -t187 * t244 + (qJD(4) * t189 + t250) * t180;
t126 = -qJD(4) * t153 - t181 * t250;
t256 = qJD(2) * t192;
t249 = t187 * t256;
t69 = Icges(6,4) * t126 + Icges(6,2) * t249 - Icges(6,6) * t125;
t209 = t191 * t69 + t256 * t88;
t65 = Icges(6,5) * t126 + Icges(6,6) * t249 - Icges(6,3) * t125;
t73 = Icges(6,1) * t126 + Icges(6,4) * t249 - Icges(6,5) * t125;
t14 = -t125 * t84 + t126 * t92 + t153 * t65 + t154 * t73 + t187 * t209;
t248 = t189 * t257;
t255 = qJD(4) * t180;
t127 = t180 * t248 - t187 * t255 - t189 * t244;
t128 = -qJD(4) * t155 - t181 * t248;
t247 = t189 * t256;
t70 = Icges(6,4) * t128 + Icges(6,2) * t247 - Icges(6,6) * t127;
t208 = t191 * t70 + t256 * t89;
t66 = Icges(6,5) * t128 + Icges(6,6) * t247 - Icges(6,3) * t127;
t74 = Icges(6,1) * t128 + Icges(6,4) * t247 - Icges(6,5) * t127;
t15 = -t125 * t85 + t126 * t93 + t153 * t66 + t154 * t74 + t187 * t208;
t67 = Icges(5,5) * t126 + Icges(5,6) * t125 + Icges(5,3) * t249;
t211 = t191 * t67 + t256 * t86;
t71 = Icges(5,4) * t126 + Icges(5,2) * t125 + Icges(5,6) * t249;
t75 = Icges(5,1) * t126 + Icges(5,4) * t125 + Icges(5,5) * t249;
t16 = t125 * t90 + t126 * t94 - t153 * t71 + t154 * t75 + t187 * t211;
t68 = Icges(5,5) * t128 + Icges(5,6) * t127 + Icges(5,3) * t247;
t210 = t191 * t68 + t256 * t87;
t72 = Icges(5,4) * t128 + Icges(5,2) * t127 + Icges(5,6) * t247;
t76 = Icges(5,1) * t128 + Icges(5,4) * t127 + Icges(5,5) * t247;
t17 = t125 * t91 + t126 * t95 - t153 * t72 + t154 * t76 + t187 * t210;
t321 = (t331 * t125 - t126 * t339 - t333 * t153 - t332 * t154) * t192 + ((t15 + t17) * t189 + (t14 + t16 + t329) * t187) * t191 + (((t328 + t316) * t187 + t327) * t192 + t334 * t191) * qJD(2);
t320 = rSges(6,1) + pkin(4);
t319 = (-qJD(2) * t335 + t67 + t69) * t192 + ((-t73 - t75) * t181 + (-t65 + t71) * t180 + (t180 * t341 - t181 * t343) * qJD(4) - t342 * qJD(2)) * t191;
t318 = (-qJD(2) * t336 + t68 + t70) * t192 + ((-t74 - t76) * t181 + (-t66 + t72) * t180 + (t180 * t345 - t181 * t346) * qJD(4) - t347 * qJD(2)) * t191;
t317 = rSges(6,3) + qJ(5);
t37 = t155 * t85 + t156 * t93 + t269 * t89;
t39 = -t155 * t91 + t156 * t95 + t269 * t87;
t315 = t37 + t39;
t314 = t191 * t322 - t328;
t310 = t187 * t324 + t189 * t323;
t186 = sin(pkin(8));
t188 = cos(pkin(8));
t236 = rSges(4,1) * t188 - rSges(4,2) * t186;
t309 = rSges(4,3) * t192 - t191 * t236;
t293 = pkin(3) * t188;
t307 = pkin(6) * t192 - t191 * t293;
t304 = 2 * m(5);
t303 = 2 * m(6);
t302 = 0.2e1 * qJD(2);
t301 = m(4) / 0.2e1;
t300 = m(5) / 0.2e1;
t299 = m(6) / 0.2e1;
t298 = t187 / 0.2e1;
t292 = rSges(6,2) * t249 + qJD(5) * t153 - t125 * t317 + t126 * t320;
t291 = rSges(6,2) * t247 + qJD(5) * t155 - t127 * t317 + t128 * t320;
t36 = t155 * t84 + t156 * t92 + t269 * t88;
t290 = t187 * t36;
t38 = -t155 * t90 + t156 * t94 + t269 * t86;
t289 = t187 * t38;
t234 = rSges(6,1) * t181 + rSges(6,3) * t180;
t285 = (-rSges(6,1) * t180 + rSges(6,3) * t181) * t253 + (rSges(6,2) * t191 + t192 * t234) * qJD(2) + (pkin(4) * t256 + qJ(5) * t253) * t181 + (qJ(5) * t256 + (-pkin(4) * qJD(4) + qJD(5)) * t191) * t180;
t284 = rSges(6,2) * t272 + t153 * t317 + t154 * t320;
t283 = rSges(6,2) * t269 + t155 * t317 + t156 * t320;
t273 = t187 * t186;
t270 = t189 * t186;
t233 = pkin(2) * t192 + qJ(3) * t191;
t167 = qJD(2) * t233 - qJD(3) * t192;
t200 = pkin(6) * t191 + t192 * t293;
t264 = -qJD(2) * t200 - t167;
t176 = pkin(2) * t191 - qJ(3) * t192;
t263 = t311 * (-qJD(2) * t176 + qJD(3) * t191);
t262 = -t176 + t307;
t261 = -(rSges(4,3) * t191 + t192 * t236) * qJD(2) - t167;
t260 = -rSges(6,2) * t192 + (pkin(4) * t181 + qJ(5) * t180 + t234) * t191;
t259 = -t176 + t309;
t258 = t311 * t233;
t235 = rSges(5,1) * t181 - rSges(5,2) * t180;
t107 = (-rSges(5,1) * t180 - rSges(5,2) * t181) * t253 + (rSges(5,3) * t191 + t192 * t235) * qJD(2);
t252 = -t107 + t264;
t148 = -rSges(5,3) * t192 + t191 * t235;
t251 = -t148 + t262;
t246 = t191 * t256;
t245 = t181 * t253;
t243 = t187 * t260;
t242 = t189 * t260;
t241 = t264 - t285;
t240 = t187 * (-pkin(3) * t270 + t187 * t200) + t189 * (pkin(3) * t273 + t189 * t200) + t258;
t239 = t307 * t308 + t263;
t238 = -t260 + t262;
t97 = rSges(5,1) * t154 - rSges(5,2) * t153 + rSges(5,3) * t272;
t99 = rSges(5,1) * t156 - rSges(5,2) * t155 + rSges(5,3) * t269;
t226 = -t187 * t99 + t189 * t97;
t201 = qJD(2) * (-Icges(3,5) * t191 - Icges(3,6) * t192);
t199 = -t187 * t283 + t189 * t284;
t195 = qJD(2) * (Icges(4,5) * t192 + (-Icges(4,1) * t188 + Icges(4,4) * t186) * t191);
t194 = qJD(2) * (Icges(4,6) * t192 + (-Icges(4,4) * t188 + Icges(4,2) * t186) * t191);
t185 = t191 ^ 2;
t171 = t188 * t268 + t273;
t170 = -t186 * t268 + t187 * t188;
t169 = t188 * t271 - t270;
t168 = -t186 * t271 - t188 * t189;
t162 = t189 * t201;
t161 = t187 * t201;
t134 = t189 * t195;
t133 = t187 * t195;
t132 = t189 * t194;
t131 = t187 * t194;
t124 = t259 * t189;
t123 = t259 * t187;
t109 = t261 * t189;
t108 = t261 * t187;
t82 = t251 * t189;
t81 = t251 * t187;
t80 = rSges(5,1) * t128 + rSges(5,2) * t127 + rSges(5,3) * t247;
t78 = rSges(5,1) * t126 + rSges(5,2) * t125 + rSges(5,3) * t249;
t64 = -t148 * t269 - t192 * t99;
t63 = t148 * t272 + t192 * t97;
t62 = t238 * t189;
t61 = t238 * t187;
t56 = t252 * t189;
t55 = t252 * t187;
t54 = t226 * t191;
t53 = t141 * t269 - t143 * t155 + t145 * t156;
t52 = t140 * t155 + t142 * t269 + t144 * t156;
t49 = t308 * t309 + t263;
t48 = t187 * (rSges(4,1) * t169 + rSges(4,2) * t168 + rSges(4,3) * t272) + t189 * (rSges(4,1) * t171 + rSges(4,2) * t170 + rSges(4,3) * t269) + t258;
t47 = -t191 * t242 - t192 * t283;
t46 = t191 * t243 + t192 * t284;
t45 = t241 * t189;
t44 = t241 * t187;
t31 = t199 * t191;
t30 = -t107 * t269 - t192 * t80 + (-t148 * t268 + t191 * t99) * qJD(2);
t29 = t107 * t272 + t192 * t78 + (t148 * t271 - t191 * t97) * qJD(2);
t28 = t187 * t97 + t189 * t99 + t240;
t27 = (-t187 * t80 + t189 * t78) * t191 + t226 * t256;
t26 = t187 * t78 + t189 * t80 + t239;
t25 = t187 * t284 + t189 * t283 + t240;
t24 = -t291 * t192 - t285 * t269 + (t191 * t283 - t192 * t242) * qJD(2);
t23 = t292 * t192 + t285 * t272 + (-t191 * t284 + t192 * t243) * qJD(2);
t22 = t187 * t292 + t189 * t291 + t239;
t21 = t127 * t91 + t128 * t95 - t155 * t72 + t156 * t76 + t189 * t210;
t20 = t127 * t90 + t128 * t94 - t155 * t71 + t156 * t75 + t189 * t211;
t19 = -t127 * t85 + t128 * t93 + t155 * t66 + t156 * t74 + t189 * t208;
t18 = -t127 * t84 + t128 * t92 + t155 * t65 + t156 * t73 + t189 * t209;
t9 = (-t187 * t291 + t189 * t292) * t191 + t199 * t256;
t8 = t187 * t21 - t189 * t20;
t7 = -t18 * t189 + t187 * t19;
t6 = -t16 * t189 + t17 * t187;
t5 = -t14 * t189 + t15 * t187;
t4 = -(-t103 * t155 + t105 * t156 + t127 * t143 + t128 * t145) * t192 + (t20 * t187 + (t21 - t267) * t189) * t191 + (t53 * t191 + (t289 + (t39 - t276) * t189) * t192) * qJD(2);
t3 = -(t155 * t100 + t156 * t104 - t127 * t140 + t128 * t144) * t192 + (t18 * t187 + (t19 - t266) * t189) * t191 + (t52 * t191 + (t290 + (t37 - t275) * t189) * t192) * qJD(2);
t1 = [0; m(4) * t49 + m(5) * t26 + m(6) * t22 - t325; (t22 * t25 + t44 * t61 + t45 * t62) * t303 + (t26 * t28 + t55 * t81 + t56 * t82) * t304 + 0.2e1 * m(4) * (t108 * t123 + t109 * t124 + t48 * t49) + (-t183 * t161 - t5 - t6 + (t131 * t168 + t133 * t169) * t189) * t189 + (t8 + t7 + t182 * t162 + (t170 * t132 + t171 * t134) * t187 + (-t170 * t131 - t168 * t132 - t171 * t133 - t169 * t134 - t161 * t187 + t162 * t189) * t189) * t187 + 0.2e1 * t326 * (rSges(3,1) * t192 - rSges(3,2) * t191) * t325; (m(4) + m(5) + m(6)) * t257; m(6) * (-t192 * t22 + t269 * t45 + t272 * t44) + m(5) * (-t192 * t26 + t269 * t56 + t272 * t55) + m(4) * (t108 * t272 + t109 * t269 - t192 * t49) + ((t191 * t25 + t268 * t62 + t271 * t61) * t299 + (t191 * t28 + t268 * t82 + t271 * t81) * t300 + (t123 * t271 + t124 * t268 + t191 * t48) * t301) * t302; -0.4e1 * (t301 + t300 + t299) * t326 * t246; m(5) * t27 + m(6) * t9; t4 * t298 + t3 * t298 + m(6) * (t22 * t31 + t23 * t62 + t24 * t61 + t25 * t9 + t44 * t47 + t45 * t46) + m(5) * (t26 * t54 + t27 * t28 + t29 * t82 + t30 * t81 + t55 * t64 + t56 * t63) + ((t8 / 0.2e1 + t7 / 0.2e1) * t189 + (t6 / 0.2e1 + t5 / 0.2e1) * t187) * t191 + (((t187 * t344 - t316 * t189) * t298 + ((-t36 - t38) * t189 + t315 * t187) * t189 / 0.2e1) * t192 + (t323 * t187 - t324 * t189) * t191 / 0.2e1) * qJD(2) - t321 * t189 / 0.2e1 - (-t318 * t187 + t319 * t189) * t192 / 0.2e1; m(5) * (-t192 * t27 + t269 * t29 + t272 * t30) + m(6) * (-t192 * t9 + t23 * t269 + t24 * t272) + ((t191 * t54 + t268 * t63 + t271 * t64) * t300 + (t191 * t31 + t268 * t46 + t271 * t47) * t299) * t302; (t23 * t46 + t24 * t47 + t31 * t9) * t303 + (t27 * t54 + t29 * t63 + t30 * t64) * t304 + t321 * t272 + (t3 + t4) * t269 + (t310 * t257 + (t187 * t316 + t327) * t249 + (t315 * t189 + t289 + t290) * t247) * t191 + (t329 * t192 + t314 * t257 - t334 * t249 + (-t52 - t53) * t247 + ((t180 * t333 + t181 * t332 + t254 * t331 - t255 * t339) * t192 + t318 * t189 + t319 * t187) * t191 + ((-t192 * t322 - t310) * t192 + (t192 * t340 + t314) * t191) * qJD(2)) * t192; (t180 * t256 + t245) * m(6); m(6) * (t25 * t245 - t125 * t61 - t127 * t62 + t153 * t44 + t155 * t45 + (t191 * t22 + t25 * t256) * t180); m(6) * ((-t125 * t187 - t127 * t189 - t244) * t191 + (t180 * t185 + (t153 * t187 + t155 * t189 - t180 * t192) * t192) * qJD(2)); m(6) * (t31 * t245 - t125 * t47 - t127 * t46 + t153 * t24 + t155 * t23 + (t191 * t9 + t256 * t31) * t180); (-t153 * t125 - t155 * t127 + (t180 * t246 + t185 * t254) * t180) * t303;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
