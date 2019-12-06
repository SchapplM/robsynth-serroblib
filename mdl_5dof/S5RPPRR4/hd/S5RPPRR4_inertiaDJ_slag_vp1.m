% Calculate time derivative of joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:05
% EndTime: 2019-12-05 17:44:25
% DurationCPUTime: 6.70s
% Computational Cost: add. (11955->501), mult. (12820->746), div. (0->0), fcn. (12113->10), ass. (0->236)
t187 = -pkin(6) - qJ(3);
t179 = -pkin(7) + t187;
t184 = sin(pkin(8));
t283 = (rSges(6,3) - t179) * t184;
t186 = cos(pkin(8));
t185 = cos(pkin(9));
t171 = t185 * pkin(3) + pkin(2);
t181 = pkin(9) + qJ(4);
t174 = cos(t181);
t163 = pkin(4) * t174 + t171;
t239 = t163 - t171;
t282 = t186 * t239;
t188 = sin(qJ(1));
t173 = sin(t181);
t183 = sin(pkin(9));
t274 = pkin(3) * t183;
t164 = pkin(4) * t173 + t274;
t220 = -t164 + t274;
t281 = t220 * t188;
t262 = Icges(5,4) * t174;
t141 = -Icges(5,6) * t186 + (-Icges(5,2) * t173 + t262) * t184;
t263 = Icges(5,4) * t173;
t142 = -Icges(5,5) * t186 + (Icges(5,1) * t174 - t263) * t184;
t280 = -t141 * t174 - t142 * t173;
t279 = t184 * (rSges(4,3) + qJ(3)) + t186 * (rSges(4,1) * t185 - rSges(4,2) * t183 + pkin(2));
t278 = 2 * m(5);
t277 = 2 * m(6);
t276 = m(5) / 0.2e1;
t275 = m(6) / 0.2e1;
t189 = cos(qJ(1));
t195 = -t184 * t187 - t282;
t248 = t188 * t173;
t270 = pkin(4) * qJD(4);
t215 = t248 * t270;
t236 = qJD(1) * t189;
t227 = t184 * t236;
t232 = t189 * t270;
t230 = t174 * t232 + t179 * t227 + t186 * t215;
t182 = qJD(4) + qJD(5);
t175 = qJ(5) + t181;
t169 = sin(t175);
t170 = cos(t175);
t249 = t188 * t170;
t251 = t186 * t189;
t199 = t169 * t251 - t249;
t200 = -t169 * t189 + t186 * t249;
t104 = qJD(1) * t199 + t182 * t200;
t250 = t188 * t169;
t150 = t170 * t189 + t186 * t250;
t153 = t170 * t251 + t250;
t105 = -qJD(1) * t153 + t150 * t182;
t245 = t105 * rSges(6,1) + t104 * rSges(6,2);
t65 = -rSges(6,3) * t227 + t245;
t273 = -t65 - (t189 * t195 + t281) * qJD(1) - t230;
t233 = t189 * t274;
t253 = t184 * t188;
t240 = t189 * t164 + t179 * t253;
t86 = t188 * t195 - t233 + t240;
t242 = -rSges(6,1) * t200 + t150 * rSges(6,2);
t98 = -rSges(6,3) * t253 + t242;
t272 = -t86 - t98;
t252 = t184 * t189;
t168 = t187 * t252;
t87 = t168 - t281 + (-t179 * t184 + t282) * t189;
t205 = -t153 * rSges(6,1) + rSges(6,2) * t199;
t99 = rSges(6,3) * t252 - t205;
t271 = -t87 - t99;
t254 = t182 * t184;
t261 = Icges(6,4) * t169;
t131 = (-Icges(6,2) * t170 - t261) * t254;
t138 = -Icges(6,5) * t186 + (Icges(6,1) * t170 - t261) * t184;
t130 = (-Icges(6,5) * t169 - Icges(6,6) * t170) * t254;
t260 = Icges(6,4) * t170;
t132 = (-Icges(6,1) * t169 - t260) * t254;
t218 = t184 * t170 * t132 - t186 * t130;
t137 = -Icges(6,6) * t186 + (-Icges(6,2) * t169 + t260) * t184;
t231 = t182 * t170 * t137;
t42 = (-t231 + (-t138 * t182 - t131) * t169) * t184 + t218;
t269 = t42 * t186;
t268 = -rSges(3,3) - qJ(2);
t266 = rSges(5,3) - t187;
t133 = (-rSges(6,1) * t169 - rSges(6,2) * t170) * t254;
t102 = qJD(1) * t150 - t153 * t182;
t103 = -qJD(1) * t200 - t182 * t199;
t206 = -t103 * rSges(6,1) - t102 * rSges(6,2);
t237 = qJD(1) * t188;
t228 = t184 * t237;
t64 = -rSges(6,3) * t228 - t206;
t264 = t133 * t252 + t186 * t64;
t139 = -t186 * rSges(6,3) + (rSges(6,1) * t170 - rSges(6,2) * t169) * t184;
t79 = t139 * t252 + t186 * t99;
t257 = t163 * t186;
t256 = t171 * t186;
t234 = qJD(4) * t184;
t147 = (-Icges(5,2) * t174 - t263) * t234;
t255 = t173 * t147;
t247 = t188 * t174;
t246 = -qJ(2) - t164;
t197 = t173 * t251 - t247;
t198 = -t173 * t189 + t186 * t247;
t118 = qJD(1) * t197 + qJD(4) * t198;
t156 = t174 * t189 + t186 * t248;
t159 = t174 * t251 + t248;
t119 = -qJD(1) * t159 + qJD(4) * t156;
t244 = t119 * rSges(5,1) + t118 * rSges(5,2);
t243 = t133 * t253 + t139 * t227;
t241 = -rSges(5,1) * t198 + t156 * rSges(5,2);
t238 = t179 - t187;
t235 = qJD(3) * t184;
t58 = Icges(6,5) * t103 + Icges(6,6) * t102 - Icges(6,3) * t228;
t60 = Icges(6,4) * t103 + Icges(6,2) * t102 - Icges(6,6) * t228;
t62 = Icges(6,1) * t103 + Icges(6,4) * t102 - Icges(6,5) * t228;
t95 = Icges(6,4) * t153 - Icges(6,2) * t199 + Icges(6,6) * t252;
t97 = Icges(6,1) * t153 - Icges(6,4) * t199 + Icges(6,5) * t252;
t10 = -t186 * t58 + ((-t182 * t95 + t62) * t170 + (-t182 * t97 - t60) * t169) * t184;
t136 = -Icges(6,3) * t186 + (Icges(6,5) * t170 - Icges(6,6) * t169) * t184;
t18 = t102 * t137 + t103 * t138 - t199 * t131 + t153 * t132 + (t130 * t189 - t136 * t237) * t184;
t92 = -Icges(6,5) * t200 + Icges(6,6) * t150 - Icges(6,3) * t253;
t94 = -Icges(6,4) * t200 + Icges(6,2) * t150 - Icges(6,6) * t253;
t96 = -Icges(6,1) * t200 + Icges(6,4) * t150 - Icges(6,5) * t253;
t26 = t153 * t96 - t199 * t94 + t252 * t92;
t93 = Icges(6,5) * t153 - Icges(6,6) * t199 + Icges(6,3) * t252;
t27 = t153 * t97 - t199 * t95 + t252 * t93;
t40 = -t186 * t92 + (-t169 * t94 + t170 * t96) * t184;
t41 = -t186 * t93 + (-t169 * t95 + t170 * t97) * t184;
t59 = Icges(6,5) * t105 + Icges(6,6) * t104 - Icges(6,3) * t227;
t61 = Icges(6,4) * t105 + Icges(6,2) * t104 - Icges(6,6) * t227;
t63 = Icges(6,1) * t105 + Icges(6,4) * t104 - Icges(6,5) * t227;
t9 = -t186 * t59 + ((-t182 * t94 + t63) * t170 + (-t182 * t96 - t61) * t169) * t184;
t229 = -t186 * (-t269 + (t10 * t189 - t188 * t9 + (-t188 * t41 - t189 * t40) * qJD(1)) * t184) + (-t18 * t186 + (-(t102 * t94 + t103 * t96 + t153 * t63 - t199 * t61) * t188 - t26 * t236 + (t102 * t95 + t103 * t97 + t153 * t62 - t199 * t60) * t189 - t27 * t237 + (-(t189 * t59 - t237 * t92) * t188 + (t189 * t58 - t237 * t93) * t189) * t184) * t184) * t252;
t223 = -pkin(1) - t257;
t222 = -pkin(1) - t256;
t221 = -qJ(2) - t274;
t219 = t246 * t188;
t146 = (-Icges(5,5) * t173 - Icges(5,6) * t174) * t234;
t148 = (-Icges(5,1) * t173 - t262) * t234;
t217 = t184 * t174 * t148 - t186 * t146;
t216 = pkin(1) * t237 - qJD(2) * t188;
t176 = qJD(2) * t189;
t212 = -t188 * t235 + t176;
t211 = t221 * t188;
t210 = rSges(3,1) * t186 - rSges(3,2) * t184;
t209 = rSges(4,1) * t183 + rSges(4,2) * t185;
t116 = qJD(1) * t156 - qJD(4) * t159;
t192 = t197 * qJD(4);
t117 = -qJD(1) * t198 - t192;
t208 = -t117 * rSges(5,1) - t116 * rSges(5,2);
t207 = -t159 * rSges(5,1) + rSges(5,2) * t197;
t204 = -pkin(1) - t210;
t202 = -rSges(6,3) * t184 + t223;
t201 = -qJ(2) - t209;
t19 = t104 * t137 + t105 * t138 + t150 * t131 - t200 * t132 + (-t130 * t188 - t136 * t236) * t184;
t52 = -t136 * t253 + t137 * t150 - t138 * t200;
t53 = t136 * t252 - t137 * t199 + t153 * t138;
t196 = -(t19 + t9) * t253 / 0.2e1 + (t10 + t18) * t252 / 0.2e1 - ((t41 + t53) * t188 + (t40 + t52) * t189) * qJD(1) * t184 / 0.2e1;
t194 = -t184 * t266 + t222;
t193 = -t189 * t235 + t216;
t135 = t188 * t268 + t189 * t204;
t191 = -pkin(1) - t279;
t24 = t150 * t94 - t200 * t96 - t253 * t92;
t25 = t150 * t95 - t200 * t97 - t253 * t93;
t2 = -t19 * t186 + (-(t104 * t94 + t105 * t96 + t150 * t61 - t200 * t63) * t188 - t24 * t236 + (t104 * t95 + t105 * t97 + t150 * t60 - t200 * t62) * t189 - t25 * t237 + (-(-t188 * t59 - t236 * t92) * t188 + (-t188 * t58 - t236 * t93) * t189) * t184) * t184;
t5 = -t52 * t186 + (-t188 * t24 + t189 * t25) * t184;
t6 = -t53 * t186 + (-t188 * t26 + t189 * t27) * t184;
t190 = (-t188 * t2 + (-t188 * t6 - t189 * t5) * qJD(1)) * t184 + t229;
t89 = t188 * t201 + t189 * t191;
t180 = t184 ^ 2;
t178 = t189 * qJ(2);
t161 = t237 * t256;
t149 = (-rSges(5,1) * t173 - rSges(5,2) * t174) * t234;
t143 = -t186 * rSges(5,3) + (rSges(5,1) * t174 - rSges(5,2) * t173) * t184;
t140 = -Icges(5,3) * t186 + (Icges(5,5) * t174 - Icges(5,6) * t173) * t184;
t134 = t189 * rSges(3,3) + t188 * t204 + t178;
t127 = t139 * t253;
t122 = qJD(1) * t135 + t176;
t121 = (t188 * t210 + t189 * t268) * qJD(1) + t216;
t120 = t184 * t239 + t186 * t238;
t113 = rSges(5,3) * t252 - t207;
t112 = -rSges(5,3) * t253 + t241;
t111 = Icges(5,1) * t159 - Icges(5,4) * t197 + Icges(5,5) * t252;
t110 = -Icges(5,1) * t198 + Icges(5,4) * t156 - Icges(5,5) * t253;
t109 = Icges(5,4) * t159 - Icges(5,2) * t197 + Icges(5,6) * t252;
t108 = -Icges(5,4) * t198 + Icges(5,2) * t156 - Icges(5,6) * t253;
t107 = Icges(5,5) * t159 - Icges(5,6) * t197 + Icges(5,3) * t252;
t106 = -Icges(5,5) * t198 + Icges(5,6) * t156 - Icges(5,3) * t253;
t90 = t98 * t228;
t88 = t188 * t191 + t189 * t209 + t178;
t85 = t186 * t113 + t143 * t252;
t84 = -t112 * t186 + t143 * t253;
t83 = t168 + t211 + (-rSges(5,3) * t184 + t222) * t189 + t207;
t82 = t188 * t194 + t178 + t233 + t241;
t81 = qJD(1) * t89 + t212;
t80 = (t279 * t188 + t201 * t189) * qJD(1) + t193;
t78 = -t186 * t98 + t127;
t77 = -rSges(5,3) * t227 + t244;
t76 = -rSges(5,3) * t228 - t208;
t75 = Icges(5,1) * t119 + Icges(5,4) * t118 - Icges(5,5) * t227;
t74 = Icges(5,1) * t117 + Icges(5,4) * t116 - Icges(5,5) * t228;
t73 = Icges(5,4) * t119 + Icges(5,2) * t118 - Icges(5,6) * t227;
t72 = Icges(5,4) * t117 + Icges(5,2) * t116 - Icges(5,6) * t228;
t71 = Icges(5,5) * t119 + Icges(5,6) * t118 - Icges(5,3) * t227;
t70 = Icges(5,5) * t117 + Icges(5,6) * t116 - Icges(5,3) * t228;
t68 = t161 - pkin(4) * t192 + (-t220 * t189 + (t184 * t238 - t257) * t188) * qJD(1);
t67 = t219 + (t223 - t283) * t189 + t205;
t66 = t188 * t202 + t178 + t240 + t242;
t57 = t140 * t252 - t141 * t197 + t159 * t142;
t56 = -t140 * t253 + t141 * t156 - t142 * t198;
t54 = (-t188 * t99 - t189 * t98) * t184;
t49 = (t189 * t194 + t211) * qJD(1) + t212 + t244;
t48 = t161 + (t189 * t221 + t253 * t266) * qJD(1) + t193 + t208;
t47 = (t280 * qJD(4) - t255) * t184 + t217;
t46 = -t186 * t77 + (t143 * t236 + t149 * t188) * t184;
t45 = t186 * t76 + (-t143 * t237 + t149 * t189) * t184;
t44 = -t186 * t107 + (-t109 * t173 + t111 * t174) * t184;
t43 = -t186 * t106 + (-t108 * t173 + t110 * t174) * t184;
t39 = -t186 * t65 + t243;
t38 = -t139 * t228 + t264;
t35 = t107 * t252 - t109 * t197 + t159 * t111;
t34 = t106 * t252 - t108 * t197 + t159 * t110;
t33 = -t107 * t253 + t109 * t156 - t111 * t198;
t32 = -t106 * t253 + t108 * t156 - t110 * t198;
t31 = (t189 * t202 + t219) * qJD(1) + t212 + t230 + t245;
t30 = t197 * t270 + (t246 * t189 + (t257 + t283) * t188) * qJD(1) + t193 + t206;
t29 = t120 * t252 + t186 * t87 + t79;
t28 = t120 * t253 + t186 * t272 + t127;
t23 = (t188 * t271 + t189 * t272) * t184;
t22 = t118 * t141 + t119 * t142 + t156 * t147 - t198 * t148 + (-t140 * t236 - t146 * t188) * t184;
t21 = t116 * t141 + t117 * t142 - t197 * t147 + t159 * t148 + (-t140 * t237 + t146 * t189) * t184;
t20 = (-t188 * t76 - t189 * t77 + (t112 * t188 - t113 * t189) * qJD(1)) * t184;
t15 = t120 * t227 - t180 * t215 + t186 * t273 + t243;
t14 = -t180 * t173 * t232 + t186 * t68 + (-t120 - t139) * t228 + t264;
t13 = t90 + (-t188 * t64 + (-qJD(1) * t99 - t65) * t189) * t184;
t12 = -t186 * t70 + (-t173 * t72 + t174 * t74 + (-t109 * t174 - t111 * t173) * qJD(4)) * t184;
t11 = -t186 * t71 + (-t173 * t73 + t174 * t75 + (-t108 * t174 - t110 * t173) * qJD(4)) * t184;
t4 = t90 + ((qJD(1) * t86 - t64 - t68) * t188 + (qJD(1) * t271 + t273) * t189) * t184;
t1 = [0.2e1 * m(3) * (t121 * t135 + t122 * t134) + 0.2e1 * m(4) * (t80 * t89 + t81 * t88) + (t48 * t83 + t49 * t82) * t278 - t169 * t138 * t254 + (t30 * t67 + t31 * t66) * t277 + t217 + t218 + t280 * t234 + (-t169 * t131 - t231 - t255) * t184; m(3) * (t189 * t121 + t188 * t122 + (t134 * t189 - t135 * t188) * qJD(1)) + m(4) * (t188 * t81 + t189 * t80 + (-t188 * t89 + t189 * t88) * qJD(1)) + m(5) * (t188 * t49 + t189 * t48 + (-t188 * t83 + t189 * t82) * qJD(1)) + m(6) * (t188 * t31 + t189 * t30 + (-t188 * t67 + t189 * t66) * qJD(1)); 0; 0.2e1 * (m(4) * (-t188 * t80 + t189 * t81 - t236 * t89 - t237 * t88) / 0.2e1 + (-t188 * t48 + t189 * t49 - t236 * t83 - t237 * t82) * t276 + (-t188 * t30 + t189 * t31 - t236 * t67 - t237 * t66) * t275) * t184; 0; 0; (-t47 - t42) * t186 + m(5) * (t45 * t83 + t46 * t82 + t48 * t85 + t49 * t84) + m(6) * (t14 * t67 + t15 * t66 + t28 * t31 + t29 * t30) + ((t12 / 0.2e1 + t21 / 0.2e1) * t189 + (-t11 / 0.2e1 - t22 / 0.2e1) * t188 + ((-t43 / 0.2e1 - t56 / 0.2e1) * t189 + (-t44 / 0.2e1 - t57 / 0.2e1) * t188) * qJD(1)) * t184 + t196; m(5) * (t46 * t188 + t45 * t189 + (-t188 * t85 + t189 * t84) * qJD(1)) + m(6) * (t14 * t189 + t15 * t188 + (-t188 * t29 + t189 * t28) * qJD(1)); 0.2e1 * (-m(5) * t20 / 0.2e1 - m(6) * t4 / 0.2e1) * t186 + 0.2e1 * ((-t188 * t45 + t189 * t46 - t236 * t85 - t237 * t84) * t276 + (-t14 * t188 + t15 * t189 - t236 * t29 - t237 * t28) * t275) * t184; (t85 * t45 + t84 * t46 + (-t112 * t189 - t113 * t188) * t20 * t184) * t278 - t186 * (-t47 * t186 + (-t11 * t188 + t12 * t189 + (-t188 * t44 - t189 * t43) * qJD(1)) * t184) + (-t21 * t186 + (-(t116 * t108 + t117 * t110 + t159 * t75 - t197 * t73) * t188 - t34 * t236 + (t116 * t109 + t117 * t111 + t159 * t74 - t197 * t72) * t189 - t35 * t237 + (-(-t106 * t237 + t189 * t71) * t188 + (-t107 * t237 + t189 * t70) * t189) * t184) * t184) * t252 + (t14 * t29 + t15 * t28 + t23 * t4) * t277 + t229 + (t22 * t186 - (-(t118 * t108 + t119 * t110 + t156 * t73 - t198 * t75) * t188 - t32 * t236 + (t118 * t109 + t119 * t111 + t156 * t72 - t198 * t74) * t189 - t33 * t237 + (-(-t106 * t236 - t188 * t71) * t188 + (-t107 * t236 - t188 * t70) * t189) * t184) * t184 - t2) * t253 + (t57 * t186 - (-t188 * t34 + t189 * t35) * t184 - t6) * t228 + (t56 * t186 - (-t188 * t32 + t189 * t33) * t184 - t5) * t227; -t269 + m(6) * (t30 * t79 + t31 * t78 + t38 * t67 + t39 * t66) + t196; m(6) * (t39 * t188 + t38 * t189 + (-t188 * t79 + t189 * t78) * qJD(1)); m(6) * (-t13 * t186 + (-t188 * t38 + t189 * t39 + (-t188 * t78 - t189 * t79) * qJD(1)) * t184); m(6) * (t13 * t23 + t14 * t79 + t15 * t78 + t28 * t39 + t29 * t38 + t4 * t54) + t190; (t13 * t54 + t38 * t79 + t39 * t78) * t277 + t190;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
