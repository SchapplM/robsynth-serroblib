% Calculate time derivative of joint inertia matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR8_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:49
% EndTime: 2019-12-31 17:08:00
% DurationCPUTime: 6.85s
% Computational Cost: add. (4066->451), mult. (11375->658), div. (0->0), fcn. (10700->6), ass. (0->227)
t169 = sin(qJ(2));
t172 = cos(qJ(2));
t265 = Icges(4,5) * t172;
t269 = Icges(3,4) * t172;
t304 = -t265 + t269 + (Icges(3,1) + Icges(4,1)) * t169;
t303 = t304 * qJD(2);
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t200 = t168 * t172 - t169 * t171;
t293 = qJD(2) - qJD(4);
t296 = t293 * t200;
t302 = t296 / 0.2e1;
t199 = t168 * t169 + t171 * t172;
t87 = t293 * t199;
t301 = -t87 / 0.2e1;
t170 = sin(qJ(1));
t300 = t170 / 0.2e1;
t173 = cos(qJ(1));
t299 = -t173 / 0.2e1;
t298 = -qJD(1) / 0.2e1;
t297 = qJD(1) / 0.2e1;
t295 = Icges(5,5) * t301 + Icges(5,6) * t302;
t208 = -Icges(3,2) * t169 + t269;
t104 = Icges(3,6) * t170 + t173 * t208;
t270 = Icges(3,4) * t169;
t212 = Icges(3,1) * t172 - t270;
t108 = Icges(3,5) * t170 + t173 * t212;
t201 = t104 * t169 - t108 * t172;
t191 = t201 * t170;
t103 = -Icges(3,6) * t173 + t170 * t208;
t107 = -Icges(3,5) * t173 + t170 * t212;
t202 = t103 * t169 - t107 * t172;
t192 = t202 * t173;
t266 = Icges(4,5) * t169;
t210 = Icges(4,1) * t172 + t266;
t106 = Icges(4,4) * t170 + t173 * t210;
t204 = Icges(4,3) * t169 + t265;
t98 = Icges(4,6) * t170 + t173 * t204;
t213 = t106 * t172 + t169 * t98;
t193 = t213 * t170;
t105 = -Icges(4,4) * t173 + t170 * t210;
t97 = -Icges(4,6) * t173 + t170 * t204;
t214 = t105 * t172 + t169 * t97;
t194 = t214 * t173;
t257 = t169 * t173;
t294 = -rSges(3,2) * t257 + t170 * rSges(3,3);
t205 = Icges(3,5) * t172 - Icges(3,6) * t169;
t99 = -Icges(3,3) * t173 + t170 * t205;
t206 = Icges(4,4) * t172 + Icges(4,6) * t169;
t101 = -Icges(4,2) * t173 + t170 * t206;
t292 = 2 * m(3);
t291 = 2 * m(4);
t290 = 2 * m(5);
t166 = t170 ^ 2;
t167 = t173 ^ 2;
t289 = m(4) / 0.2e1;
t288 = m(5) / 0.2e1;
t287 = -pkin(2) - pkin(3);
t286 = -rSges(4,1) - pkin(2);
t285 = -rSges(5,3) - pkin(6);
t146 = rSges(3,1) * t169 + rSges(3,2) * t172;
t284 = m(3) * t146;
t283 = pkin(3) * t169;
t282 = pkin(6) * t170;
t281 = pkin(6) * t173;
t244 = qJD(1) * t170;
t54 = t173 * t87 + t200 * t244;
t55 = t173 * t296 - t199 * t244;
t279 = t55 * rSges(5,1) + t54 * rSges(5,2);
t278 = rSges(4,1) * t169;
t277 = rSges(4,2) * t173;
t276 = rSges(3,3) * t173;
t162 = t170 * rSges(4,2);
t120 = t200 * t173;
t121 = t199 * t173;
t65 = Icges(5,5) * t121 - Icges(5,6) * t120 - Icges(5,3) * t170;
t275 = t170 * t65;
t274 = t173 * t65;
t273 = -rSges(4,3) - qJ(3);
t118 = t200 * t170;
t119 = t199 * t170;
t220 = -t119 * rSges(5,1) + t118 * rSges(5,2);
t70 = rSges(5,3) * t173 - t220;
t272 = t170 * t172 * pkin(3) + t281 + t70;
t256 = t172 * t173;
t156 = pkin(3) * t256;
t255 = t121 * rSges(5,1) - t120 * rSges(5,2);
t71 = -rSges(5,3) * t170 + t255;
t271 = t156 - t282 + t71;
t260 = qJ(3) * t169;
t259 = qJ(3) * t172;
t91 = -rSges(5,1) * t200 - rSges(5,2) * t199;
t258 = qJD(1) * t91;
t219 = pkin(2) * t172 + t260;
t123 = t219 * t170;
t157 = pkin(2) * t256;
t124 = qJ(3) * t257 + t157;
t254 = t170 * t123 + t173 * t124;
t122 = qJD(2) * t219 - qJD(3) * t172;
t221 = rSges(4,1) * t172 + rSges(4,3) * t169;
t253 = -t221 * qJD(2) - t122;
t144 = pkin(2) * t169 - t259;
t145 = -rSges(4,3) * t172 + t278;
t252 = -t144 - t145;
t239 = qJD(2) * t173;
t230 = t172 * t239;
t238 = qJD(3) * t169;
t251 = qJ(3) * t230 + t173 * t238;
t243 = qJD(1) * t173;
t250 = rSges(4,2) * t243 + rSges(4,3) * t230;
t229 = t169 * t244;
t249 = rSges(3,2) * t229 + rSges(3,3) * t243;
t248 = t173 * pkin(1) + t170 * pkin(5);
t247 = t166 + t167;
t100 = Icges(3,3) * t170 + t173 * t205;
t246 = qJD(1) * t100;
t102 = Icges(4,2) * t170 + t173 * t206;
t245 = qJD(1) * t102;
t242 = qJD(2) * t169;
t241 = qJD(2) * t170;
t240 = qJD(2) * t172;
t24 = Icges(5,4) * t55 + Icges(5,2) * t54 - Icges(5,6) * t243;
t25 = Icges(5,1) * t55 + Icges(5,4) * t54 - Icges(5,5) * t243;
t49 = Icges(5,4) * t87 - Icges(5,2) * t296;
t50 = Icges(5,1) * t87 - Icges(5,4) * t296;
t67 = Icges(5,4) * t121 - Icges(5,2) * t120 - Icges(5,6) * t170;
t69 = Icges(5,1) * t121 - Icges(5,4) * t120 - Icges(5,5) * t170;
t88 = -Icges(5,5) * t200 - Icges(5,6) * t199;
t89 = -Icges(5,4) * t200 - Icges(5,2) * t199;
t90 = -Icges(5,1) * t200 - Icges(5,4) * t199;
t237 = -t120 * t49 / 0.2e1 + t121 * t50 / 0.2e1 + t170 * t295 - t243 * t88 / 0.2e1 + t54 * t89 / 0.2e1 + t55 * t90 / 0.2e1 - t199 * t24 / 0.2e1 - t200 * t25 / 0.2e1 - t67 * t296 / 0.2e1 + t69 * t87 / 0.2e1;
t56 = -qJD(1) * t120 + t170 * t87;
t57 = qJD(1) * t121 + t170 * t296;
t178 = Icges(5,4) * t57 + Icges(5,2) * t56 - Icges(5,6) * t244;
t179 = Icges(5,1) * t57 + Icges(5,4) * t56 - Icges(5,5) * t244;
t66 = Icges(5,4) * t119 - Icges(5,2) * t118 + Icges(5,6) * t173;
t68 = Icges(5,1) * t119 - Icges(5,4) * t118 + Icges(5,5) * t173;
t236 = t118 * t49 / 0.2e1 - t119 * t50 / 0.2e1 + t173 * t295 + t244 * t88 / 0.2e1 - t56 * t89 / 0.2e1 - t57 * t90 / 0.2e1 + t199 * t178 / 0.2e1 + t200 * t179 / 0.2e1 + t66 * t302 + t68 * t301;
t232 = t169 * t241;
t150 = pkin(2) * t232;
t231 = t169 * t239;
t181 = -t172 * t244 - t231;
t235 = t123 * t243 + t170 * (qJD(1) * t157 + t170 * t238 - t150 + (t169 * t243 + t240 * t170) * qJ(3)) + t173 * (pkin(2) * t181 - qJ(3) * t229 + t251);
t160 = pkin(5) * t243;
t233 = t160 + t251;
t114 = rSges(4,1) * t256 + rSges(4,3) * t257 + t162;
t228 = t91 + t283;
t51 = rSges(5,1) * t87 - rSges(5,2) * t296;
t227 = t247 * t51;
t96 = t252 * t173;
t226 = t248 + t124;
t225 = -t144 - t228;
t224 = -t57 * rSges(5,1) - t56 * rSges(5,2);
t177 = Icges(5,5) * t57 + Icges(5,6) * t56 - Icges(5,3) * t244;
t64 = Icges(5,5) * t119 - Icges(5,6) * t118 + Icges(5,3) * t173;
t182 = -t120 * t66 + t121 * t68 - t170 * t64;
t19 = -t120 * t67 + t121 * t69 - t275;
t23 = Icges(5,5) * t55 + Icges(5,6) * t54 - Icges(5,3) * t243;
t1 = (-t120 * t24 + t121 * t25 - t170 * t23 + t54 * t67 + t55 * t69) * t170 - (-t120 * t178 + t121 * t179 - t170 * t177 - t243 * t64 + t54 * t66 + t55 * t68) * t173 + (t19 * t173 + (t182 - t274) * t170) * qJD(1);
t18 = -t118 * t67 + t119 * t69 + t274;
t183 = -t118 * t66 + t119 * t68 + t173 * t64;
t2 = (-t118 * t24 + t119 * t25 + t173 * t23 + t56 * t67 + t57 * t69) * t170 - (-t118 * t178 + t119 * t179 + t173 * t177 - t244 * t64 + t56 * t66 + t57 * t68) * t173 + (t18 * t173 + (t183 - t275) * t170) * qJD(1);
t223 = t170 * t1 - t173 * t2;
t222 = rSges(3,1) * t172 - rSges(3,2) * t169;
t207 = Icges(3,2) * t172 + t270;
t198 = -pkin(3) * t240 - t122 - t51;
t115 = rSges(3,1) * t256 + t294;
t197 = -pkin(1) - t222;
t59 = t225 * t173;
t195 = qJD(2) * t146;
t188 = qJD(2) * t207;
t187 = qJD(2) * (-Icges(4,4) * t169 + Icges(4,6) * t172);
t186 = qJD(2) * (-Icges(3,5) * t169 - Icges(3,6) * t172);
t184 = t172 * t287 - pkin(1) - t260;
t180 = t169 * t273 + t172 * t286 - pkin(1);
t176 = t180 * t170;
t175 = t170 * t184 + t173 * t285;
t164 = t173 * pkin(5);
t135 = t222 * qJD(2);
t125 = t144 * t244;
t113 = t170 * t222 - t276;
t112 = t170 * t221 - t277;
t95 = t252 * t170;
t94 = t115 + t248;
t93 = t170 * t197 + t164 + t276;
t79 = t170 * t187 + t245;
t78 = -qJD(1) * t101 + t173 * t187;
t77 = t170 * t186 + t246;
t76 = -qJD(1) * t99 + t173 * t186;
t73 = t226 + t114;
t72 = t164 + t176 + t277;
t61 = t146 * t241 + ((-rSges(3,3) - pkin(5)) * t170 + t197 * t173) * qJD(1);
t60 = rSges(3,1) * t181 - rSges(3,2) * t230 - pkin(1) * t244 + t160 + t249;
t58 = t225 * t170;
t47 = qJD(1) * t96 + t170 * t253;
t46 = t145 * t244 + t173 * t253 + t125;
t45 = t170 * t100 - t201 * t173;
t44 = t170 * t99 - t192;
t43 = t170 * t102 + t213 * t173;
t42 = t170 * t101 + t194;
t41 = -t100 * t173 - t191;
t40 = -t170 * t202 - t173 * t99;
t39 = -t102 * t173 + t193;
t38 = -t101 * t173 + t170 * t214;
t37 = t170 * t285 + t156 + t226 + t255;
t36 = t164 + t175 + t220;
t35 = t170 * t112 + t114 * t173 + t254;
t34 = t150 + (-t238 + (t172 * t273 + t278) * qJD(2)) * t170 + ((-rSges(4,2) - pkin(5)) * t170 + t180 * t173) * qJD(1);
t33 = qJD(1) * t176 + t231 * t286 + t233 + t250;
t32 = -t170 * t70 - t173 * t71;
t31 = -t199 * t67 - t200 * t69;
t30 = -t199 * t66 - t200 * t68;
t29 = -t120 * t89 + t121 * t90 - t170 * t88;
t28 = -t118 * t89 + t119 * t90 + t173 * t88;
t27 = -rSges(5,3) * t244 - t224;
t26 = -rSges(5,3) * t243 + t279;
t22 = qJD(1) * t59 + t170 * t198;
t21 = t173 * t198 + t228 * t244 + t125;
t20 = t170 * t272 + t173 * t271 + t254;
t15 = t150 + (-t238 + (-t259 + t283) * qJD(2)) * t170 + ((-pkin(5) - t285) * t170 + t184 * t173) * qJD(1) + t224;
t14 = qJD(1) * t175 + t231 * t287 + t233 + t279;
t13 = t173 * t250 + (-t145 * t166 - t167 * t278) * qJD(2) + (t173 * t112 + (-t114 - t124 + t162) * t170) * qJD(1) + t235;
t12 = t170 * t19 - t173 * t182;
t11 = t170 * t18 - t173 * t183;
t10 = -t170 * t27 - t173 * t26 + (t170 * t71 - t173 * t70) * qJD(1);
t5 = (-pkin(3) * t231 + t26 + (t272 - t281) * qJD(1)) * t173 + (-pkin(3) * t232 + t27 + (-t124 - t271 - t282) * qJD(1)) * t170 + t235;
t3 = [t87 * t90 - t200 * t50 - t296 * t89 - t199 * t49 + (t14 * t37 + t15 * t36) * t290 + (t33 * t73 + t34 * t72) * t291 + (t60 * t94 + t61 * t93) * t292 + (-Icges(4,3) * t172 - t207 + t210 + t212 + t266) * t242 + (t208 - t204 + t304) * t240; m(5) * (t14 * t58 + t15 * t59 + t21 * t36 + t22 * t37) + m(4) * (t33 * t95 + t34 * t96 + t46 * t72 + t47 * t73) + ((t104 * t298 + t188 * t300 + t98 * t297) * t172 + (t303 * t300 + (t106 + t108) * t298) * t169 + t236 + m(3) * (-t135 * t93 - t146 * t61) + (t31 / 0.2e1 + t29 / 0.2e1 - t94 * t284 + (-t98 / 0.2e1 + t104 / 0.2e1) * t172 + (t106 / 0.2e1 + t108 / 0.2e1) * t169) * qJD(1)) * t173 + ((t103 * t298 + t188 * t299 + t97 * t297) * t172 + (t303 * t299 + (t105 + t107) * t298) * t169 + t237 + m(3) * (-t135 * t94 - t146 * t60) + (t30 / 0.2e1 + t28 / 0.2e1 + t93 * t284 + (-t97 / 0.2e1 + t103 / 0.2e1) * t172 + (t105 / 0.2e1 + t107 / 0.2e1) * t169) * qJD(1)) * t170 + (t193 / 0.2e1 - t191 / 0.2e1 - t194 / 0.2e1 + t192 / 0.2e1 + (t205 + t206) * (t166 / 0.2e1 + t167 / 0.2e1)) * qJD(2); (t20 * t5 + t21 * t59 + t22 * t58) * t290 + (t35 * t13 + t46 * t96 + t47 * t95) * t291 + t170 * ((t170 * t76 + (t44 + t191) * qJD(1)) * t170 + (t45 * qJD(1) + (t103 * t240 + t107 * t242) * t173 + (-t77 + (-t104 * t172 - t108 * t169) * qJD(2) + (t100 - t202) * qJD(1)) * t170) * t173) - t173 * ((t173 * t79 + (t39 - t194) * qJD(1)) * t173 + (t38 * qJD(1) + (-t106 * t242 + t240 * t98 + t245) * t170 + (-t78 + (t105 * t169 - t172 * t97) * qJD(2) + t213 * qJD(1)) * t173) * t170) - t173 * ((t173 * t77 + (t41 + t192) * qJD(1)) * t173 + (t40 * qJD(1) + (-t104 * t240 - t108 * t242 + t246) * t170 + (-t76 + (t103 * t172 + t107 * t169) * qJD(2) - t201 * qJD(1)) * t173) * t170) + t170 * ((t170 * t78 + (t42 - t193) * qJD(1)) * t170 + (t43 * qJD(1) + (t105 * t242 - t240 * t97) * t173 + (-t79 + (-t106 * t169 + t172 * t98) * qJD(2) + (t102 + t214) * qJD(1)) * t170) * t173) + ((t170 * t113 + t115 * t173) * ((qJD(1) * t113 - t173 * t195 + t249) * t173 + (-t170 * t195 + (-t115 + t294) * qJD(1)) * t170) + t247 * t146 * t135) * t292 + t223 + (t11 + (-t38 - t40) * t173 + (t39 + t41) * t170) * t244 + (t12 + (-t42 - t44) * t173 + (t43 + t45) * t170) * t243; 0.2e1 * ((t170 * t37 + t173 * t36) * t288 + (t170 * t73 + t173 * t72) * t289) * t240 + 0.2e1 * ((t14 * t170 + t15 * t173 + t243 * t37 - t244 * t36) * t288 + (t170 * t33 + t173 * t34 + t243 * t73 - t244 * t72) * t289) * t169; 0.2e1 * ((t239 * t59 + t241 * t58 - t5) * t288 + (t239 * t96 + t241 * t95 - t13) * t289) * t172 + 0.2e1 * ((qJD(2) * t20 + t170 * t22 + t173 * t21 + t243 * t58 - t244 * t59) * t288 + (qJD(2) * t35 + t170 * t47 + t173 * t46 + t243 * t95 - t244 * t96) * t289) * t169; 0.4e1 * (t289 + t288) * (-0.1e1 + t247) * t169 * t240; (m(5) * (t15 * t91 + t258 * t37 + t36 * t51) - t236) * t173 + (m(5) * (t14 * t91 - t258 * t36 + t37 * t51) - t237) * t170 + ((t31 + t29) * t173 + (t30 + t28) * t170) * t298; m(5) * (t10 * t20 + t32 * t5) + (m(5) * (t21 * t91 + t258 * t58 + t51 * t59) + t2 - qJD(1) * t12) * t173 + (m(5) * (t22 * t91 - t258 * t59 + t51 * t58) - qJD(1) * t11 - t1) * t170; m(5) * (-t10 * t172 + t169 * t227 + (t172 * t247 * t91 + t169 * t32) * qJD(2)); (t32 * t10 + t227 * t91) * t290 + (t170 * t11 + t173 * t12) * qJD(1) + t223;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
