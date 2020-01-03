% Calculate time derivative of joint inertia matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:59:02
% DurationCPUTime: 6.39s
% Computational Cost: add. (1701->365), mult. (4502->514), div. (0->0), fcn. (3448->4), ass. (0->186)
t140 = sin(qJ(2));
t142 = cos(qJ(2));
t239 = Icges(4,5) * t142;
t242 = Icges(5,4) * t142;
t245 = Icges(3,4) * t142;
t278 = -t239 - t242 + t245 + (Icges(3,1) + Icges(4,1) + Icges(5,1)) * t140;
t260 = rSges(5,1) + pkin(3);
t208 = t260 * t140;
t204 = -rSges(5,2) * t142 + t208;
t141 = sin(qJ(1));
t249 = -qJ(4) - rSges(5,3);
t203 = t249 * t141;
t143 = cos(qJ(1));
t227 = t142 * t143;
t228 = t140 * t143;
t247 = rSges(5,2) * t228 + t260 * t227 + t203;
t202 = t249 * t143;
t274 = t141 / 0.2e1;
t273 = -t143 / 0.2e1;
t271 = -qJD(1) / 0.2e1;
t270 = qJD(1) / 0.2e1;
t175 = -Icges(3,2) * t140 + t245;
t69 = Icges(3,6) * t141 + t143 * t175;
t246 = Icges(3,4) * t140;
t181 = Icges(3,1) * t142 - t246;
t75 = Icges(3,5) * t141 + t143 * t181;
t183 = t140 * t69 - t142 * t75;
t159 = t183 * t141;
t68 = -Icges(3,6) * t143 + t141 * t175;
t74 = -Icges(3,5) * t143 + t141 * t181;
t184 = t140 * t68 - t142 * t74;
t160 = t184 * t143;
t172 = Icges(5,2) * t140 + t242;
t65 = -Icges(5,6) * t141 + t143 * t172;
t243 = Icges(5,4) * t140;
t177 = Icges(5,1) * t142 + t243;
t71 = -Icges(5,5) * t141 + t143 * t177;
t185 = t140 * t65 + t142 * t71;
t161 = t185 * t141;
t64 = Icges(5,6) * t143 + t141 * t172;
t70 = Icges(5,5) * t143 + t141 * t177;
t186 = t140 * t64 + t142 * t70;
t162 = t186 * t143;
t169 = Icges(4,3) * t140 + t239;
t61 = Icges(4,6) * t141 + t143 * t169;
t240 = Icges(4,5) * t140;
t179 = Icges(4,1) * t142 + t240;
t73 = Icges(4,4) * t141 + t143 * t179;
t187 = t140 * t61 + t142 * t73;
t163 = t187 * t141;
t60 = -Icges(4,6) * t143 + t141 * t169;
t72 = -Icges(4,4) * t143 + t141 * t179;
t188 = t140 * t60 + t142 * t72;
t164 = t188 * t143;
t267 = -rSges(3,2) * t228 + t141 * rSges(3,3);
t167 = Icges(5,5) * t142 + Icges(5,6) * t140;
t58 = Icges(5,3) * t143 + t141 * t167;
t170 = Icges(3,5) * t142 - Icges(3,6) * t140;
t62 = -Icges(3,3) * t143 + t141 * t170;
t173 = Icges(4,4) * t142 + Icges(4,6) * t140;
t66 = -Icges(4,2) * t143 + t141 * t173;
t266 = 2 * m(3);
t265 = 2 * m(4);
t264 = 2 * m(5);
t138 = t141 ^ 2;
t139 = t143 ^ 2;
t263 = m(4) / 0.2e1;
t262 = m(5) / 0.2e1;
t261 = -rSges(4,1) - pkin(2);
t115 = rSges(3,1) * t140 + rSges(3,2) * t142;
t259 = m(3) * t115;
t258 = pkin(2) * t140;
t195 = pkin(2) * t142 + qJ(3) * t140;
t86 = t195 * t141;
t129 = pkin(2) * t227;
t87 = qJ(3) * t228 + t129;
t257 = t141 * t86 + t143 * t87;
t197 = rSges(4,1) * t142 + rSges(4,3) * t140;
t85 = qJD(2) * t195 - qJD(3) * t142;
t256 = -t197 * qJD(2) - t85;
t255 = rSges(4,1) * t140;
t254 = rSges(4,2) * t143;
t252 = rSges(3,3) * t143;
t134 = t141 * rSges(4,2);
t251 = -rSges(5,2) - qJ(3);
t250 = -rSges(4,3) - qJ(3);
t196 = rSges(5,1) * t142 + rSges(5,2) * t140;
t248 = -t202 + (pkin(3) * t142 + t196) * t141;
t59 = -Icges(5,3) * t141 + t143 * t167;
t231 = qJD(1) * t59;
t63 = Icges(3,3) * t141 + t143 * t170;
t230 = qJD(1) * t63;
t67 = Icges(4,2) * t141 + t143 * t173;
t229 = qJD(1) * t67;
t112 = -qJ(3) * t142 + t258;
t114 = -rSges(4,3) * t142 + t255;
t226 = -t112 - t114;
t214 = qJD(2) * t143;
t205 = t142 * t214;
t213 = qJD(3) * t140;
t225 = qJ(3) * t205 + t143 * t213;
t218 = qJD(1) * t143;
t224 = rSges(4,2) * t218 + rSges(4,3) * t205;
t219 = qJD(1) * t141;
t207 = t140 * t219;
t223 = rSges(3,2) * t207 + rSges(3,3) * t218;
t221 = t143 * pkin(1) + t141 * pkin(5);
t220 = t138 + t139;
t217 = qJD(2) * t140;
t216 = qJD(2) * t141;
t215 = qJD(2) * t142;
t212 = -pkin(2) - t260;
t120 = t216 * t258;
t206 = t140 * t214;
t148 = -t142 * t219 - t206;
t211 = t141 * (qJD(1) * t129 + t141 * t213 - t120 + (t140 * t218 + t215 * t141) * qJ(3)) + t143 * (pkin(2) * t148 - qJ(3) * t207 + t225) + t86 * t218;
t132 = pkin(5) * t218;
t209 = t132 + t225;
t81 = rSges(4,1) * t227 + rSges(4,3) * t228 + t134;
t57 = t226 * t143;
t201 = t221 + t87;
t200 = (-t167 / 0.2e1 + t173 / 0.2e1 + t170 / 0.2e1) * qJD(2);
t199 = -t112 - t204;
t198 = rSges(3,1) * t142 - rSges(3,2) * t140;
t136 = t143 * pkin(5);
t146 = t140 * t251 + t142 * t212 - pkin(1);
t144 = t141 * t146 + t202;
t25 = t136 + t144;
t26 = t201 + t247;
t182 = t141 * t26 + t143 * t25;
t174 = Icges(3,2) * t142 + t246;
t166 = -pkin(3) * t215 - t196 * qJD(2) - t85;
t82 = rSges(3,1) * t227 + t267;
t165 = -pkin(1) - t198;
t52 = t199 * t143;
t158 = qJD(2) * t115;
t154 = qJD(2) * t174;
t153 = qJD(2) * (-Icges(4,4) * t140 + Icges(4,6) * t142);
t151 = qJD(2) * (-Icges(3,5) * t140 - Icges(3,6) * t142);
t149 = qJD(2) * (-Icges(5,5) * t140 + Icges(5,6) * t142);
t147 = t140 * t250 + t142 * t261 - pkin(1);
t145 = t147 * t141;
t118 = rSges(5,2) * t205;
t100 = t198 * qJD(2);
t88 = t112 * t219;
t79 = t141 * t198 - t252;
t78 = t141 * t197 - t254;
t56 = t226 * t141;
t55 = t82 + t221;
t54 = t141 * t165 + t136 + t252;
t51 = t199 * t141;
t42 = t141 * t153 + t229;
t41 = -qJD(1) * t66 + t143 * t153;
t38 = t141 * t151 + t230;
t37 = -qJD(1) * t62 + t143 * t151;
t34 = t141 * t149 + t231;
t33 = -qJD(1) * t58 + t143 * t149;
t32 = t201 + t81;
t31 = t136 + t145 + t254;
t28 = t115 * t216 + ((-rSges(3,3) - pkin(5)) * t141 + t165 * t143) * qJD(1);
t27 = rSges(3,1) * t148 - rSges(3,2) * t205 - pkin(1) * t219 + t132 + t223;
t24 = qJD(1) * t57 + t141 * t256;
t23 = t114 * t219 + t143 * t256 + t88;
t22 = t141 * t63 - t183 * t143;
t21 = t141 * t62 - t160;
t20 = t141 * t67 + t187 * t143;
t19 = t141 * t66 + t164;
t18 = -t141 * t59 + t185 * t143;
t17 = -t141 * t58 + t162;
t16 = -t143 * t63 - t159;
t15 = -t141 * t184 - t143 * t62;
t14 = -t143 * t67 + t163;
t13 = t141 * t188 - t143 * t66;
t12 = t143 * t59 + t161;
t11 = t141 * t186 + t143 * t58;
t10 = qJD(1) * t52 + t141 * t166;
t9 = t143 * t166 + t204 * t219 + t88;
t8 = t141 * t78 + t143 * t81 + t257;
t7 = t120 + (-t213 + (t142 * t250 + t255) * qJD(2)) * t141 + ((-rSges(4,2) - pkin(5)) * t141 + t147 * t143) * qJD(1);
t6 = qJD(1) * t145 + t206 * t261 + t209 + t224;
t5 = t141 * t248 + t143 * t247 + t257;
t4 = -qJD(4) * t143 + t120 + (-t213 + (t142 * t251 + t208) * qJD(2)) * t141 + ((-pkin(5) - t249) * t141 + t146 * t143) * qJD(1);
t3 = qJD(1) * t144 - qJD(4) * t141 + t206 * t212 + t118 + t209;
t2 = t143 * t224 + (-t114 * t138 - t139 * t255) * qJD(2) + (t143 * t78 + (-t81 - t87 + t134) * t141) * qJD(1) + t211;
t1 = t143 * t118 + (-t138 * t204 - t139 * t208) * qJD(2) + ((t202 + t248) * t143 + (-t87 + t203 - t247) * t141) * qJD(1) + t211;
t29 = [(t27 * t55 + t28 * t54) * t266 + (t31 * t7 + t32 * t6) * t265 + (t25 * t4 + t26 * t3) * t264 + (-t174 + t177 + t179 + t181 + t240 + t243 + (-Icges(5,2) - Icges(4,3)) * t142) * t217 + (t175 - t169 - t172 + t278) * t215; m(5) * (t10 * t26 + t25 * t9 + t3 * t51 + t4 * t52) + m(4) * (t23 * t31 + t24 * t32 + t56 * t6 + t57 * t7) + (m(3) * (-t100 * t54 - t115 * t28) + t200 * t143 + (t154 * t274 + t69 * t271 + (t61 + t65) * t270) * t142) * t143 + (m(3) * (-t100 * t55 - t115 * t27) + t200 * t141 + (t154 * t273 + t68 * t271 + (t60 + t64) * t270) * t142) * t141 + (((t71 + t73 + t75) * t143 + (t70 + t72 + t74) * t141) * t271 + (t141 * t273 + t143 * t274) * t278 * qJD(2)) * t140 + (-t162 / 0.2e1 - t164 / 0.2e1 + t160 / 0.2e1 - t159 / 0.2e1 + t161 / 0.2e1 + t163 / 0.2e1) * qJD(2) + ((-t55 * t259 + (-t65 / 0.2e1 - t61 / 0.2e1 + t69 / 0.2e1) * t142 + (t71 / 0.2e1 + t73 / 0.2e1 + t75 / 0.2e1) * t140) * t143 + (t54 * t259 + (-t64 / 0.2e1 - t60 / 0.2e1 + t68 / 0.2e1) * t142 + (t70 / 0.2e1 + t72 / 0.2e1 + t74 / 0.2e1) * t140) * t141) * qJD(1); (t2 * t8 + t23 * t57 + t24 * t56) * t265 + (t1 * t5 + t10 * t51 + t52 * t9) * t264 - t143 * ((-t143 * t34 + (t12 - t162) * qJD(1)) * t143 + (t11 * qJD(1) + (t215 * t65 - t217 * t71 - t231) * t141 + (t33 + (t140 * t70 - t142 * t64) * qJD(2) + t185 * qJD(1)) * t143) * t141) - t143 * ((t143 * t42 + (t14 - t164) * qJD(1)) * t143 + (t13 * qJD(1) + (t215 * t61 - t217 * t73 + t229) * t141 + (-t41 + (t140 * t72 - t142 * t60) * qJD(2) + t187 * qJD(1)) * t143) * t141) + t141 * ((t141 * t37 + (t21 + t159) * qJD(1)) * t141 + (t22 * qJD(1) + (t215 * t68 + t217 * t74) * t143 + (-t38 + (-t140 * t75 - t142 * t69) * qJD(2) + (-t184 + t63) * qJD(1)) * t141) * t143) - t143 * ((t143 * t38 + (t16 + t160) * qJD(1)) * t143 + (t15 * qJD(1) + (-t215 * t69 - t217 * t75 + t230) * t141 + (-t37 + (t140 * t74 + t142 * t68) * qJD(2) - t183 * qJD(1)) * t143) * t141) + ((t141 * t79 + t143 * t82) * ((qJD(1) * t79 - t143 * t158 + t223) * t143 + (-t141 * t158 + (-t82 + t267) * qJD(1)) * t141) + t220 * t115 * t100) * t266 + t141 * ((t141 * t41 + (t19 - t163) * qJD(1)) * t141 + (t20 * qJD(1) + (-t215 * t60 + t217 * t72) * t143 + (-t42 + (-t140 * t73 + t142 * t61) * qJD(2) + (t188 + t67) * qJD(1)) * t141) * t143) + t141 * ((-t141 * t33 + (t17 - t161) * qJD(1)) * t141 + (t18 * qJD(1) + (-t215 * t64 + t217 * t70) * t143 + (t34 + (-t140 * t71 + t142 * t65) * qJD(2) + (t186 - t59) * qJD(1)) * t141) * t143) + ((-t11 - t13 - t15) * t143 + (t12 + t14 + t16) * t141) * t219 + ((-t17 - t19 - t21) * t143 + (t18 + t20 + t22) * t141) * t218; 0.2e1 * ((t141 * t32 + t143 * t31) * t263 + t182 * t262) * t215 + 0.2e1 * ((t141 * t6 + t143 * t7 + t218 * t32 - t219 * t31) * t263 + (t141 * t3 + t143 * t4 + t218 * t26 - t219 * t25) * t262) * t140; 0.2e1 * ((t214 * t57 + t216 * t56 - t2) * t263 + (t214 * t52 + t216 * t51 - t1) * t262) * t142 + 0.2e1 * ((qJD(2) * t8 + t141 * t24 + t143 * t23 + t218 * t56 - t219 * t57) * t263 + (qJD(2) * t5 + t10 * t141 + t143 * t9 + t218 * t51 - t219 * t52) * t262) * t140; 0.4e1 * (t263 + t262) * (-0.1e1 + t220) * t140 * t215; m(5) * (-qJD(1) * t182 - t141 * t4 + t143 * t3); m(5) * (t10 * t143 - t141 * t9 + (-t141 * t51 - t143 * t52) * qJD(1)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t29(1), t29(2), t29(4), t29(7); t29(2), t29(3), t29(5), t29(8); t29(4), t29(5), t29(6), t29(9); t29(7), t29(8), t29(9), t29(10);];
Mq = res;
