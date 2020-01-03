% Calculate time derivative of joint inertia matrix for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:22
% DurationCPUTime: 5.94s
% Computational Cost: add. (4216->382), mult. (12156->581), div. (0->0), fcn. (11746->6), ass. (0->197)
t148 = sin(pkin(6));
t149 = cos(pkin(6));
t152 = cos(qJ(3));
t150 = sin(qJ(3));
t153 = cos(qJ(2));
t210 = t150 * t153;
t137 = -t148 * t152 + t149 * t210;
t209 = t152 * t153;
t138 = t148 * t150 + t149 * t209;
t151 = sin(qJ(2));
t211 = t149 * t151;
t81 = Icges(5,5) * t138 + Icges(5,6) * t211 + Icges(5,3) * t137;
t87 = Icges(4,4) * t138 - Icges(4,2) * t137 + Icges(4,6) * t211;
t281 = t81 - t87;
t83 = Icges(4,5) * t138 - Icges(4,6) * t137 + Icges(4,3) * t211;
t85 = Icges(5,4) * t138 + Icges(5,2) * t211 + Icges(5,6) * t137;
t282 = t83 + t85;
t89 = Icges(5,1) * t138 + Icges(5,4) * t211 + Icges(5,5) * t137;
t91 = Icges(4,1) * t138 - Icges(4,4) * t137 + Icges(4,5) * t211;
t280 = t89 + t91;
t135 = t148 * t210 + t149 * t152;
t212 = t149 * t150;
t136 = t148 * t209 - t212;
t213 = t148 * t151;
t80 = Icges(5,5) * t136 + Icges(5,6) * t213 + Icges(5,3) * t135;
t86 = Icges(4,4) * t136 - Icges(4,2) * t135 + Icges(4,6) * t213;
t278 = t80 - t86;
t88 = Icges(5,1) * t136 + Icges(5,4) * t213 + Icges(5,5) * t135;
t90 = Icges(4,1) * t136 - Icges(4,4) * t135 + Icges(4,5) * t213;
t276 = t88 + t90;
t279 = t281 * t135 + t280 * t136 + t282 * t213;
t82 = Icges(4,5) * t136 - Icges(4,6) * t135 + Icges(4,3) * t213;
t84 = Icges(5,4) * t136 + Icges(5,2) * t213 + Icges(5,6) * t135;
t277 = t82 + t84;
t218 = Icges(5,5) * t152;
t170 = Icges(5,3) * t150 + t218;
t119 = -Icges(5,6) * t153 + t151 * t170;
t220 = Icges(4,4) * t152;
t173 = -Icges(4,2) * t150 + t220;
t122 = -Icges(4,6) * t153 + t151 * t173;
t266 = t119 - t122;
t171 = Icges(4,5) * t152 - Icges(4,6) * t150;
t120 = -Icges(4,3) * t153 + t151 * t171;
t172 = Icges(5,4) * t152 + Icges(5,6) * t150;
t121 = -Icges(5,2) * t153 + t151 * t172;
t275 = t120 + t121;
t219 = Icges(5,5) * t150;
t175 = Icges(5,1) * t152 + t219;
t123 = -Icges(5,4) * t153 + t151 * t175;
t221 = Icges(4,4) * t150;
t176 = Icges(4,1) * t152 - t221;
t124 = -Icges(4,5) * t153 + t151 * t176;
t274 = t123 + t124;
t271 = t281 * t150 + t280 * t152;
t270 = t278 * t150 + t276 * t152;
t253 = t278 * t135 + t276 * t136 + t277 * t213;
t269 = t266 * t135 + t274 * t136 + t275 * t213;
t200 = qJD(3) * t151;
t100 = (-Icges(4,2) * t152 - t221) * t200 + (Icges(4,6) * t151 + t153 * t173) * qJD(2);
t97 = (Icges(5,3) * t152 - t219) * t200 + (Icges(5,6) * t151 + t153 * t170) * qJD(2);
t268 = t100 - t97;
t101 = (-Icges(5,1) * t150 + t218) * t200 + (Icges(5,4) * t151 + t153 * t175) * qJD(2);
t102 = (-Icges(4,1) * t150 - t220) * t200 + (Icges(4,5) * t151 + t153 * t176) * qJD(2);
t267 = t101 + t102;
t216 = t121 * t153;
t217 = t120 * t153;
t264 = -t216 - t217;
t227 = t153 * ((-Icges(5,4) * t150 + Icges(5,6) * t152) * t200 + (Icges(5,2) * t151 + t153 * t172) * qJD(2));
t228 = t153 * ((-Icges(4,5) * t150 - Icges(4,6) * t152) * t200 + (Icges(4,3) * t151 + t153 * t171) * qJD(2));
t263 = -t228 - t227;
t262 = t279 * t149;
t261 = t270 * t151 - t153 * t277;
t260 = t271 * t151 - t153 * t282;
t259 = -t266 * t150 - t152 * t274;
t240 = t149 ^ 2;
t247 = t148 ^ 2 + t240;
t199 = qJD(3) * t152;
t191 = t153 * t199;
t198 = t148 * qJD(2);
t109 = -t148 * t191 + (qJD(3) * t149 + t151 * t198) * t150;
t203 = qJD(2) * t151;
t193 = t152 * t203;
t110 = -qJD(3) * t135 - t148 * t193;
t202 = qJD(2) * t153;
t195 = t153 * t198;
t63 = Icges(5,4) * t110 + Icges(5,2) * t195 - Icges(5,6) * t109;
t162 = t151 * t63 + t202 * t84;
t59 = Icges(5,5) * t110 + Icges(5,6) * t195 - Icges(5,3) * t109;
t67 = Icges(5,1) * t110 + Icges(5,4) * t195 - Icges(5,5) * t109;
t14 = -t109 * t80 + t110 * t88 + t135 * t59 + t136 * t67 + t148 * t162;
t201 = qJD(3) * t150;
t111 = -t148 * t201 - t149 * t191 + t203 * t212;
t112 = -qJD(3) * t137 - t149 * t193;
t194 = t149 * t202;
t64 = Icges(5,4) * t112 + Icges(5,2) * t194 - Icges(5,6) * t111;
t161 = t151 * t64 + t202 * t85;
t60 = Icges(5,5) * t112 + Icges(5,6) * t194 - Icges(5,3) * t111;
t68 = Icges(5,1) * t112 + Icges(5,4) * t194 - Icges(5,5) * t111;
t15 = -t109 * t81 + t110 * t89 + t135 * t60 + t136 * t68 + t148 * t161;
t61 = Icges(4,5) * t110 + Icges(4,6) * t109 + Icges(4,3) * t195;
t164 = t151 * t61 + t202 * t82;
t65 = Icges(4,4) * t110 + Icges(4,2) * t109 + Icges(4,6) * t195;
t69 = Icges(4,1) * t110 + Icges(4,4) * t109 + Icges(4,5) * t195;
t16 = t109 * t86 + t110 * t90 - t135 * t65 + t136 * t69 + t148 * t164;
t62 = Icges(4,5) * t112 + Icges(4,6) * t111 + Icges(4,3) * t194;
t163 = t151 * t62 + t202 * t83;
t66 = Icges(4,4) * t112 + Icges(4,2) * t111 + Icges(4,6) * t194;
t70 = Icges(4,1) * t112 + Icges(4,4) * t111 + Icges(4,5) * t194;
t17 = t109 * t87 + t110 * t91 - t135 * t66 + t136 * t70 + t148 * t163;
t258 = (t266 * t109 - t110 * t274 + t268 * t135 - t267 * t136) * t153 + ((t15 + t17) * t149 + (t14 + t16 + t263) * t148) * t151 + (((t264 + t253) * t148 + t262) * t153 + t269 * t151) * qJD(2);
t257 = rSges(5,1) + pkin(3);
t256 = (-t270 * qJD(2) + t61 + t63) * t153 + ((-t67 - t69) * t152 + (-t59 + t65) * t150 + (t276 * t150 - t278 * t152) * qJD(3) - t277 * qJD(2)) * t151;
t255 = (-t271 * qJD(2) + t62 + t64) * t153 + ((-t68 - t70) * t152 + (-t60 + t66) * t150 + (t280 * t150 - t281 * t152) * qJD(3) - t282 * qJD(2)) * t151;
t254 = rSges(5,3) + qJ(4);
t36 = t137 * t81 + t138 * t89 + t211 * t85;
t38 = -t137 * t87 + t138 * t91 + t211 * t83;
t252 = t36 + t38;
t251 = t259 * t151 - t264;
t248 = t261 * t148 + t260 * t149;
t246 = qJD(2) * (t151 * rSges(3,1) + rSges(3,2) * t153);
t243 = 2 * m(4);
t242 = 2 * m(5);
t239 = t148 / 0.2e1;
t234 = rSges(5,2) * t195 + qJD(4) * t135 - t254 * t109 + t257 * t110;
t233 = rSges(5,2) * t194 + qJD(4) * t137 - t254 * t111 + t257 * t112;
t35 = t137 * t80 + t138 * t88 + t211 * t84;
t232 = t148 * t35;
t37 = -t137 * t86 + t138 * t90 + t211 * t82;
t231 = t148 * t37;
t185 = rSges(5,1) * t152 + rSges(5,3) * t150;
t226 = (-rSges(5,1) * t150 + rSges(5,3) * t152) * t200 + (rSges(5,2) * t151 + t153 * t185) * qJD(2) + (pkin(3) * t202 + qJ(4) * t200) * t152 + (qJ(4) * t202 + (-pkin(3) * qJD(3) + qJD(4)) * t151) * t150;
t225 = rSges(5,2) * t213 + t254 * t135 + t257 * t136;
t224 = rSges(5,2) * t211 + t254 * t137 + t257 * t138;
t186 = rSges(4,1) * t152 - rSges(4,2) * t150;
t126 = -rSges(4,3) * t153 + t151 * t186;
t215 = t126 * t153;
t104 = (-rSges(4,1) * t150 - rSges(4,2) * t152) * t200 + (rSges(4,3) * t151 + t153 * t186) * qJD(2);
t188 = pkin(2) * t153 + pkin(5) * t151;
t141 = t188 * qJD(2);
t208 = -t104 - t141;
t145 = t151 * pkin(2) - pkin(5) * t153;
t207 = t247 * qJD(2) * t145;
t206 = -rSges(5,2) * t153 + (pkin(3) * t152 + qJ(4) * t150 + t185) * t151;
t205 = -t126 - t145;
t204 = t247 * t188;
t197 = -t141 - t226;
t196 = -t145 - t206;
t192 = t151 * t199;
t190 = t148 * t206;
t189 = t149 * t206;
t93 = rSges(4,1) * t136 - rSges(4,2) * t135 + rSges(4,3) * t213;
t95 = rSges(4,1) * t138 - rSges(4,2) * t137 + rSges(4,3) * t211;
t182 = -t148 * t95 + t149 * t93;
t156 = qJD(2) * (-Icges(3,5) * t151 - Icges(3,6) * t153);
t155 = t150 * t202 + t192;
t154 = -t224 * t148 + t225 * t149;
t129 = t148 * t156;
t106 = t205 * t149;
t105 = t205 * t148;
t79 = t247 * t246;
t78 = t196 * t149;
t77 = t196 * t148;
t76 = t208 * t149;
t75 = t208 * t148;
t74 = t112 * rSges(4,1) + t111 * rSges(4,2) + rSges(4,3) * t194;
t72 = t110 * rSges(4,1) + t109 * rSges(4,2) + rSges(4,3) * t195;
t58 = -t126 * t211 - t153 * t95;
t57 = t126 * t213 + t153 * t93;
t52 = t120 * t211 - t122 * t137 + t124 * t138;
t51 = t119 * t137 + t121 * t211 + t123 * t138;
t48 = t182 * t151;
t47 = t197 * t149;
t46 = t197 * t148;
t45 = t148 * t93 + t149 * t95 + t204;
t44 = -t151 * t189 - t224 * t153;
t43 = t151 * t190 + t225 * t153;
t30 = t154 * t151;
t29 = t148 * t72 + t149 * t74 - t207;
t28 = -t104 * t211 - t153 * t74 + (-t149 * t215 + t151 * t95) * qJD(2);
t27 = t104 * t213 + t153 * t72 + (t148 * t215 - t151 * t93) * qJD(2);
t26 = t225 * t148 + t224 * t149 + t204;
t25 = (-t148 * t74 + t149 * t72) * t151 + t182 * t202;
t24 = t234 * t148 + t233 * t149 - t207;
t23 = -t233 * t153 - t226 * t211 + (t224 * t151 - t153 * t189) * qJD(2);
t22 = t234 * t153 + t226 * t213 + (-t225 * t151 + t153 * t190) * qJD(2);
t21 = t111 * t87 + t112 * t91 - t137 * t66 + t138 * t70 + t149 * t163;
t20 = t111 * t86 + t112 * t90 - t137 * t65 + t138 * t69 + t149 * t164;
t19 = -t111 * t81 + t112 * t89 + t137 * t60 + t138 * t68 + t149 * t161;
t18 = -t111 * t80 + t112 * t88 + t137 * t59 + t138 * t67 + t149 * t162;
t9 = (-t233 * t148 + t234 * t149) * t151 + t154 * t202;
t8 = t148 * t21 - t149 * t20;
t7 = t148 * t19 - t149 * t18;
t6 = t148 * t17 - t149 * t16;
t5 = -t14 * t149 + t148 * t15;
t4 = -(-t137 * t100 + t138 * t102 + t111 * t122 + t112 * t124) * t153 + (t20 * t148 + (t21 - t228) * t149) * t151 + (t52 * t151 + (t231 + (t38 - t217) * t149) * t153) * qJD(2);
t3 = -(t138 * t101 - t111 * t119 + t112 * t123 + t137 * t97) * t153 + (t18 * t148 + (t19 - t227) * t149) * t151 + (t51 * t151 + (t232 + (t36 - t216) * t149) * t153) * qJD(2);
t1 = [0; -m(3) * t79 + m(4) * t29 + m(5) * t24; (t26 * t24 + t46 * t77 + t47 * t78) * t242 + (t105 * t75 + t106 * t76 + t45 * t29) * t243 + 0.2e1 * m(3) * (t246 - t79) * t247 * (rSges(3,1) * t153 - rSges(3,2) * t151) + (-t240 * t129 - t5 - t6) * t149 + (t7 + t8 + (-t148 * t129 + t247 * t156) * t149) * t148; m(4) * t25 + m(5) * t9; t4 * t239 + t3 * t239 + m(5) * (t22 * t78 + t23 * t77 + t30 * t24 + t9 * t26 + t43 * t47 + t44 * t46) + m(4) * (t105 * t28 + t106 * t27 + t25 * t45 + t48 * t29 + t57 * t76 + t58 * t75) + ((t8 / 0.2e1 + t7 / 0.2e1) * t149 + (t6 / 0.2e1 + t5 / 0.2e1) * t148) * t151 + (((t279 * t148 - t253 * t149) * t239 + ((-t35 - t37) * t149 + t252 * t148) * t149 / 0.2e1) * t153 + (t260 * t148 - t261 * t149) * t151 / 0.2e1) * qJD(2) - t258 * t149 / 0.2e1 - (-t255 * t148 + t256 * t149) * t153 / 0.2e1; (t48 * t25 + t27 * t57 + t28 * t58) * t243 + (t22 * t43 + t23 * t44 + t30 * t9) * t242 + t258 * t213 + (t3 + t4) * t211 + (t248 * t203 + (t253 * t148 + t262) * t195 + (t252 * t149 + t231 + t232) * t194) * t151 + (t263 * t153 + t251 * t203 - t269 * t195 + (-t52 - t51) * t194 + ((-t268 * t150 + t267 * t152 + t266 * t199 - t201 * t274) * t153 + t255 * t149 + t256 * t148) * t151 + ((-t259 * t153 - t248) * t153 + (t275 * t153 + t251) * t151) * qJD(2)) * t153; t155 * m(5); m(5) * (t26 * t192 - t109 * t77 - t111 * t78 + t135 * t46 + t137 * t47 + (t151 * t24 + t202 * t26) * t150); m(5) * (t30 * t192 - t109 * t44 - t111 * t43 + t135 * t23 + t137 * t22 + (t151 * t9 + t202 * t30) * t150); (t150 * t151 * t155 - t135 * t109 - t137 * t111) * t242;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
