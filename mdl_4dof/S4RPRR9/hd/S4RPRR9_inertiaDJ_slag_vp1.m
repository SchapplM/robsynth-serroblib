% Calculate time derivative of joint inertia matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR9_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:09
% EndTime: 2019-12-31 16:56:15
% DurationCPUTime: 3.70s
% Computational Cost: add. (4515->416), mult. (12275->630), div. (0->0), fcn. (11725->6), ass. (0->212)
t143 = sin(qJ(3));
t144 = sin(qJ(1));
t146 = cos(qJ(3));
t203 = qJD(3) * t146;
t188 = t144 * t203;
t147 = cos(qJ(1));
t206 = qJD(1) * t147;
t265 = t143 * t206 + t188;
t247 = rSges(5,3) + pkin(6);
t192 = t247 * t146;
t245 = pkin(3) * t143;
t264 = -t192 + t245;
t142 = sin(qJ(4));
t145 = cos(qJ(4));
t223 = Icges(5,4) * t145;
t160 = -Icges(5,2) * t142 + t223;
t93 = Icges(5,6) * t143 + t146 * t160;
t224 = Icges(5,4) * t142;
t163 = Icges(5,1) * t145 - t224;
t96 = Icges(5,5) * t143 + t146 * t163;
t263 = -t142 * t93 + t145 * t96;
t226 = Icges(4,4) * t143;
t161 = Icges(4,2) * t146 + t226;
t94 = Icges(4,6) * t147 + t144 * t161;
t225 = Icges(4,4) * t146;
t164 = Icges(4,1) * t143 + t225;
t97 = Icges(4,5) * t147 + t144 * t164;
t172 = t143 * t97 + t146 * t94;
t262 = t147 * t172;
t179 = rSges(4,1) * t143 + rSges(4,2) * t146;
t154 = t147 * t179;
t261 = -t142 * t96 - t145 * t93;
t191 = t146 * t206;
t195 = rSges(4,1) * t265 + rSges(4,2) * t191;
t253 = -pkin(1) - pkin(5);
t200 = -rSges(4,3) + t253;
t205 = qJD(3) * t143;
t209 = qJ(2) * t206 + qJD(2) * t144;
t54 = (-rSges(4,2) * t205 + qJD(1) * t200) * t144 + t195 + t209;
t242 = rSges(4,2) * t143;
t120 = rSges(4,1) * t146 - t242;
t134 = qJD(2) * t147;
t202 = qJD(3) * t147;
t55 = t134 + t120 * t202 + (t200 * t147 + (-qJ(2) - t179) * t144) * qJD(1);
t260 = t144 * t55 - t147 * t54;
t182 = qJD(1) * t143 + qJD(4);
t259 = t144 * t182 - t146 * t202;
t158 = Icges(4,5) * t143 + Icges(4,6) * t146;
t258 = -Icges(4,3) * t144 + t147 * t158;
t257 = -Icges(4,6) * t144 + t147 * t161;
t256 = -Icges(4,5) * t144 + t147 * t164;
t255 = 2 * m(4);
t254 = 2 * m(5);
t140 = t144 ^ 2;
t141 = t147 ^ 2;
t252 = t143 / 0.2e1;
t251 = -t144 / 0.2e1;
t249 = t147 / 0.2e1;
t248 = rSges(3,2) - pkin(1);
t246 = m(4) * t120;
t244 = pkin(3) * t146;
t157 = Icges(5,5) * t145 - Icges(5,6) * t142;
t201 = qJD(4) * t146;
t72 = (-Icges(5,5) * t142 - Icges(5,6) * t145) * t201 + (Icges(5,3) * t146 - t143 * t157) * qJD(3);
t78 = (-Icges(5,1) * t142 - t223) * t201 + (Icges(5,5) * t146 - t143 * t163) * qJD(3);
t90 = Icges(5,3) * t143 + t146 * t157;
t148 = t146 * t145 * t78 + t143 * t72 + t90 * t203 - t205 * t263;
t75 = (-Icges(5,2) * t145 - t224) * t201 + (Icges(5,6) * t146 - t143 * t160) * qJD(3);
t241 = t142 * t75;
t47 = t143 * t90 + t146 * t263;
t243 = ((t261 * qJD(4) - t241) * t146 + t148) * t143 + t47 * t203;
t238 = t143 * t94;
t237 = t143 * t257;
t236 = t144 * rSges(4,3);
t232 = t146 * t97;
t231 = t146 * t256;
t137 = t147 * rSges(4,3);
t217 = t143 * t144;
t131 = pkin(3) * t217;
t214 = t144 * t146;
t213 = t145 * t147;
t216 = t144 * t142;
t106 = -t143 * t216 + t213;
t215 = t144 * t145;
t218 = t142 * t147;
t107 = t143 * t215 + t218;
t210 = t107 * rSges(5,1) + t106 * rSges(5,2);
t70 = -rSges(5,3) * t214 + t210;
t229 = pkin(6) * t214 - t131 - t70;
t180 = pkin(6) * t146 - t245;
t108 = t143 * t218 + t215;
t109 = -t143 * t213 + t216;
t178 = -t109 * rSges(5,1) - t108 * rSges(5,2);
t212 = t146 * t147;
t71 = rSges(5,3) * t212 - t178;
t228 = t180 * t147 + t71;
t177 = rSges(5,1) * t145 - rSges(5,2) * t142;
t81 = (-rSges(5,1) * t142 - rSges(5,2) * t145) * t201 + (rSges(5,3) * t146 - t143 * t177) * qJD(3);
t227 = t180 * qJD(3) + t81;
t91 = Icges(4,3) * t147 + t144 * t158;
t219 = qJD(1) * t91;
t100 = rSges(5,3) * t143 + t146 * t177;
t121 = pkin(6) * t143 + t244;
t211 = t100 + t121;
t208 = t147 * pkin(1) + t144 * qJ(2);
t207 = qJD(1) * t144;
t204 = qJD(3) * t144;
t190 = t143 * t204;
t183 = qJD(4) * t143 + qJD(1);
t61 = -t183 * t215 + (-t147 * t182 - t188) * t142;
t62 = t182 * t213 + (-t142 * t183 + t145 * t203) * t144;
t199 = t62 * rSges(5,1) + t61 * rSges(5,2) + rSges(5,3) * t190;
t66 = Icges(5,4) * t107 + Icges(5,2) * t106 - Icges(5,6) * t214;
t68 = Icges(5,1) * t107 + Icges(5,4) * t106 - Icges(5,5) * t214;
t176 = t142 * t66 - t145 * t68;
t64 = Icges(5,5) * t107 + Icges(5,6) * t106 - Icges(5,3) * t214;
t28 = t143 * t64 - t146 * t176;
t39 = t106 * t93 + t107 * t96 - t214 * t90;
t198 = t28 / 0.2e1 + t39 / 0.2e1;
t67 = Icges(5,4) * t109 + Icges(5,2) * t108 + Icges(5,6) * t212;
t69 = Icges(5,1) * t109 + Icges(5,4) * t108 + Icges(5,5) * t212;
t175 = t142 * t67 - t145 * t69;
t65 = Icges(5,5) * t109 + Icges(5,6) * t108 + Icges(5,3) * t212;
t29 = t143 * t65 - t146 * t175;
t40 = t108 * t93 + t109 * t96 + t212 * t90;
t197 = -t29 / 0.2e1 - t40 / 0.2e1;
t196 = t253 * t144;
t194 = -pkin(3) * t265 - pkin(6) * t190;
t99 = rSges(4,1) * t217 + rSges(4,2) * t214 + t137;
t193 = t147 * pkin(5) + t208;
t189 = t143 * t202;
t115 = t179 * qJD(3);
t185 = (t140 + t141) * t115;
t184 = qJD(1) * t211;
t156 = t183 * t147;
t59 = -t142 * t259 + t145 * t156;
t60 = t142 * t156 + t145 * t259;
t181 = t60 * rSges(5,1) + t59 * rSges(5,2);
t171 = -t143 * t256 - t146 * t257;
t24 = t106 * t66 + t107 * t68 - t214 * t64;
t25 = t106 * t67 + t107 * t69 - t214 * t65;
t17 = t25 * t144 + t147 * t24;
t170 = t144 * t24 - t147 * t25;
t26 = t108 * t66 + t109 * t68 + t212 * t64;
t27 = t108 * t67 + t109 * t69 + t212 * t65;
t18 = t27 * t144 + t147 * t26;
t169 = t144 * t26 - t147 * t27;
t168 = t29 * t144 + t147 * t28;
t167 = t144 * t28 - t147 * t29;
t166 = t144 * t71 + t147 * t70;
t165 = Icges(4,1) * t146 - t226;
t162 = -Icges(4,2) * t143 + t225;
t159 = Icges(4,5) * t146 - Icges(4,6) * t143;
t155 = rSges(3,3) * t147 + t144 * t248;
t153 = t171 * t144;
t152 = qJD(3) * t165;
t151 = qJD(3) * t162;
t150 = t190 - t191;
t149 = -t146 * t207 - t189;
t136 = t147 * qJ(2);
t103 = -rSges(3,2) * t147 + t144 * rSges(3,3) + t208;
t102 = t136 + t155;
t101 = t236 - t154;
t87 = t134 + (t248 * t147 + (-rSges(3,3) - qJ(2)) * t144) * qJD(1);
t86 = qJD(1) * t155 + t209;
t85 = t193 + t99;
t84 = t144 * t200 + t136 + t154;
t83 = t211 * t147;
t82 = t211 * t144;
t74 = qJD(1) * t258 + t159 * t204;
t73 = -t159 * t202 + t219;
t53 = -t144 * t192 + t131 + t193 + t210;
t52 = t147 * t264 + t136 + t178 + t196;
t51 = t100 * t212 - t143 * t71;
t50 = t100 * t214 + t143 * t70;
t49 = -t144 * t258 - t147 * t171;
t48 = t144 * t91 - t262;
t46 = -t147 * t258 + t153;
t45 = t172 * t144 + t147 * t91;
t43 = t166 * t146;
t42 = t144 * t227 + t147 * t184;
t41 = t144 * t184 - t147 * t227;
t38 = -rSges(5,3) * t191 + t199;
t37 = rSges(5,3) * t149 + t181;
t36 = Icges(5,1) * t62 + Icges(5,4) * t61 + Icges(5,5) * t150;
t35 = Icges(5,1) * t60 + Icges(5,4) * t59 + Icges(5,5) * t149;
t34 = Icges(5,4) * t62 + Icges(5,2) * t61 + Icges(5,6) * t150;
t33 = Icges(5,4) * t60 + Icges(5,2) * t59 + Icges(5,6) * t149;
t32 = Icges(5,5) * t62 + Icges(5,6) * t61 + Icges(5,3) * t150;
t31 = Icges(5,5) * t60 + Icges(5,6) * t59 + Icges(5,3) * t149;
t30 = t144 * t229 + t147 * t228;
t23 = t134 + (t143 * t247 + t244) * t202 + (t253 * t147 + (-qJ(2) - t264) * t144) * qJD(1) - t181;
t22 = (-t147 * t192 + t196) * qJD(1) - t194 + t199 + t209;
t21 = (-t100 * t204 + t38) * t143 + (qJD(3) * t70 + t100 * t206 + t144 * t81) * t146;
t20 = (-t100 * t202 - t37) * t143 + (-qJD(3) * t71 - t100 * t207 + t147 * t81) * t146;
t16 = -t90 * t191 + t106 * t75 + t107 * t78 + t61 * t93 + t62 * t96 + (-t146 * t72 + t205 * t90) * t144;
t15 = -t90 * t189 + t108 * t75 + t109 * t78 + t59 * t93 + t60 * t96 + (t147 * t72 - t207 * t90) * t146;
t14 = t166 * t205 + (-t144 * t37 - t147 * t38 + (t144 * t70 - t147 * t71) * qJD(1)) * t146;
t13 = (-qJD(1) * t228 + t194 - t38) * t144 + (t37 - t121 * t202 + (t131 + t229) * qJD(1)) * t147;
t12 = t40 * t143 - t146 * t169;
t11 = t39 * t143 - t146 * t170;
t10 = (qJD(3) * t175 + t31) * t143 + (qJD(3) * t65 - t142 * t33 + t145 * t35 + (-t142 * t69 - t145 * t67) * qJD(4)) * t146;
t9 = (qJD(3) * t176 + t32) * t143 + (qJD(3) * t64 - t142 * t34 + t145 * t36 + (-t142 * t68 - t145 * t66) * qJD(4)) * t146;
t8 = -t65 * t191 + t106 * t33 + t107 * t35 + t61 * t67 + t62 * t69 + (-t146 * t31 + t205 * t65) * t144;
t7 = -t64 * t191 + t106 * t34 + t107 * t36 + t61 * t66 + t62 * t68 + (-t146 * t32 + t205 * t64) * t144;
t6 = -t65 * t189 + t108 * t33 + t109 * t35 + t59 * t67 + t60 * t69 + (t147 * t31 - t207 * t65) * t146;
t5 = -t64 * t189 + t108 * t34 + t109 * t36 + t59 * t66 + t60 * t68 + (t147 * t32 - t207 * t64) * t146;
t4 = -qJD(1) * t170 + t8 * t144 + t147 * t7;
t3 = -qJD(1) * t169 + t6 * t144 + t147 * t5;
t2 = (qJD(3) * t170 + t16) * t143 + (-qJD(1) * t17 + qJD(3) * t39 - t144 * t7 + t147 * t8) * t146;
t1 = (qJD(3) * t169 + t15) * t143 + (-qJD(1) * t18 + qJD(3) * t40 - t144 * t5 + t147 * t6) * t146;
t19 = [0.2e1 * m(3) * (t102 * t87 + t103 * t86) - t143 * t152 - t164 * t203 + t161 * t205 + (t54 * t85 + t55 * t84) * t255 + (t22 * t53 + t23 * t52) * t254 + t148 + t261 * t201 + (-t151 - t241) * t146; m(3) * (t144 * t87 - t147 * t86 + (t102 * t147 + t103 * t144) * qJD(1)) + m(4) * ((t144 * t85 + t147 * t84) * qJD(1) + t260) + m(5) * (t144 * t23 - t147 * t22 + (t144 * t53 + t147 * t52) * qJD(1)); 0; m(4) * (t260 * t120 - (t144 * t84 - t147 * t85) * t115) + m(5) * (-t22 * t83 + t23 * t82 + t41 * t53 + t42 * t52) - (t141 / 0.2e1 + t140 / 0.2e1) * t158 * qJD(3) + ((t238 / 0.2e1 - t232 / 0.2e1 + t85 * t246 - t198) * t144 + (t237 / 0.2e1 - t231 / 0.2e1 + t84 * t246 - t197) * t147) * qJD(1) + (-qJD(3) * t171 - t143 * (qJD(1) * t94 - t162 * t202) + t146 * (qJD(1) * t97 - t165 * t202) + t10 + t15) * t144 / 0.2e1 + (-qJD(3) * t172 - t143 * (qJD(1) * t257 + t144 * t151) + t146 * (qJD(1) * t256 + t144 * t152) + t16 + t9) * t249; m(5) * (t42 * t144 - t41 * t147 + (-t144 * t83 + t147 * t82) * qJD(1)) - m(4) * t185; ((t101 * t147 - t144 * t99) * (-t144 * t195 + (-t120 * t141 + t140 * t242) * qJD(3) + ((-t99 + t137) * t147 + (-t101 + t154 + t236) * t144) * qJD(1)) - t120 * t185) * t255 + t147 * ((t147 * t74 + (t46 + t262) * qJD(1)) * t147 + (-t45 * qJD(1) + (-t203 * t256 + t205 * t257) * t144 + (t73 + (t232 - t238) * qJD(3) + (t171 - t91) * qJD(1)) * t147) * t144) + t144 * ((t144 * t73 + (-t48 + t153) * qJD(1)) * t144 + (t49 * qJD(1) + (-t203 * t97 + t205 * t94 + t219) * t147 + (t74 + (t231 - t237) * qJD(3) + t172 * qJD(1)) * t144) * t147) + (t30 * t13 - t41 * t83 + t42 * t82) * t254 + t147 * t4 + t144 * t3 + (-t46 * t144 - t147 * t45 - t17) * t207 + (t49 * t144 + t147 * t48 + t18) * t206; m(5) * (t20 * t52 + t21 * t53 + t22 * t50 + t23 * t51) + (t144 * t198 + t147 * t197) * t205 + ((t10 / 0.2e1 + t15 / 0.2e1) * t147 + (-t9 / 0.2e1 - t16 / 0.2e1) * t144 + (t144 * t197 - t147 * t198) * qJD(1)) * t146 + t243; m(5) * (t20 * t144 - t147 * t21 + (t144 * t50 + t147 * t51) * qJD(1)); m(5) * (-t43 * t13 + t14 * t30 + t20 * t82 - t21 * t83 + t50 * t41 + t51 * t42) + (t2 / 0.2e1 + qJD(1) * t12 / 0.2e1 - t18 * t205 / 0.2e1 + (qJD(1) * t29 + t9) * t252) * t147 + (-qJD(1) * t11 / 0.2e1 + t17 * t205 / 0.2e1 + t1 / 0.2e1 + (-qJD(1) * t28 + t10) * t252) * t144 + (t4 * t251 + t3 * t249 + qJD(3) * t168 / 0.2e1 + (-t147 * t17 / 0.2e1 + t18 * t251) * qJD(1)) * t146; (-t14 * t43 + t20 * t51 + t21 * t50) * t254 + ((t144 * t11 - t147 * t12 + t143 * t167) * qJD(3) + t243) * t143 + (-t144 * t2 + t147 * t1 + t143 * (t10 * t147 - t9 * t144) + (t47 * t143 - t146 * t167) * qJD(3) + (-t147 * t11 - t144 * t12 - t143 * t168) * qJD(1)) * t146;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1), t19(2), t19(4), t19(7); t19(2), t19(3), t19(5), t19(8); t19(4), t19(5), t19(6), t19(9); t19(7), t19(8), t19(9), t19(10);];
Mq = res;
