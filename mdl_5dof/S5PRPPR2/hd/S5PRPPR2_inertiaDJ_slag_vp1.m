% Calculate time derivative of joint inertia matrix for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:59
% EndTime: 2019-12-05 15:24:10
% DurationCPUTime: 4.18s
% Computational Cost: add. (6587->294), mult. (7475->495), div. (0->0), fcn. (6903->10), ass. (0->171)
t145 = sin(qJ(2));
t146 = cos(qJ(2));
t140 = sin(pkin(7));
t135 = t140 ^ 2;
t142 = cos(pkin(7));
t136 = t142 ^ 2;
t263 = t135 + t136;
t259 = t263 * qJD(2);
t266 = m(3) * (t145 * rSges(3,1) + rSges(3,2) * t146) * t259;
t265 = 0.1e1 - t263;
t138 = qJ(2) + pkin(8);
t132 = sin(t138);
t134 = cos(t138);
t264 = (-Icges(3,5) * t145 - Icges(4,5) * t132 - Icges(3,6) * t146 - Icges(4,6) * t134) * qJD(2);
t262 = 2 * m(4);
t261 = t264 * t140;
t260 = t264 * t142;
t139 = sin(pkin(9));
t141 = cos(pkin(9));
t192 = rSges(5,1) * t141 - rSges(5,2) * t139;
t252 = t134 * rSges(5,3) - t192 * t132;
t144 = -pkin(6) - qJ(4);
t214 = qJ(4) + t144;
t129 = pkin(4) * t141 + pkin(3);
t238 = pkin(3) - t129;
t249 = t238 * t132 - t214 * t134;
t245 = 2 * m(6);
t244 = m(5) / 0.2e1;
t243 = m(6) / 0.2e1;
t242 = t140 / 0.2e1;
t241 = t142 / 0.2e1;
t240 = pkin(2) * t145;
t236 = t263 * pkin(2) * t146;
t235 = pkin(2) * qJD(2);
t137 = pkin(9) + qJ(5);
t131 = sin(t137);
t133 = cos(t137);
t175 = Icges(6,5) * t133 - Icges(6,6) * t131;
t208 = qJD(5) * t132;
t233 = t134 * ((-Icges(6,5) * t131 - Icges(6,6) * t133) * t208 + (Icges(6,3) * t132 + t175 * t134) * qJD(2));
t81 = -Icges(6,3) * t134 + t175 * t132;
t232 = t134 * t81;
t219 = t134 * t142;
t106 = -t131 * t219 + t133 * t140;
t107 = t131 * t140 + t133 * t219;
t221 = t132 * t142;
t220 = t134 * t140;
t104 = -t131 * t220 - t133 * t142;
t105 = -t131 * t142 + t133 * t220;
t222 = t132 * t140;
t49 = Icges(6,5) * t105 + Icges(6,6) * t104 + Icges(6,3) * t222;
t51 = Icges(6,4) * t105 + Icges(6,2) * t104 + Icges(6,6) * t222;
t53 = Icges(6,1) * t105 + Icges(6,4) * t104 + Icges(6,5) * t222;
t18 = t106 * t51 + t107 * t53 + t49 * t221;
t231 = t140 * t18;
t50 = Icges(6,5) * t107 + Icges(6,6) * t106 + Icges(6,3) * t221;
t52 = Icges(6,4) * t107 + Icges(6,2) * t106 + Icges(6,6) * t221;
t54 = Icges(6,1) * t107 + Icges(6,4) * t106 + Icges(6,5) * t221;
t17 = t104 * t52 + t105 * t54 + t50 * t222;
t230 = t142 * t17;
t225 = Icges(6,4) * t131;
t224 = Icges(6,4) * t133;
t218 = t139 * t140;
t217 = t139 * t142;
t216 = t140 * t141;
t215 = t141 * t142;
t213 = t263 * t145 * t235;
t212 = qJD(2) * t132;
t211 = qJD(2) * t134;
t210 = qJD(2) * t140;
t209 = qJD(2) * t142;
t206 = t146 * t235;
t205 = t132 * t210;
t204 = t132 * t209;
t203 = t134 * t210;
t202 = t134 * t209;
t125 = t132 * pkin(3) - t134 * qJ(4);
t201 = -t125 - t240;
t126 = rSges(4,1) * t132 + rSges(4,2) * t134;
t200 = -t126 - t240;
t190 = pkin(3) * t134 + qJ(4) * t132;
t199 = t190 * t263 + t236;
t198 = -t213 + t263 * (-t125 * qJD(2) + qJD(4) * t132);
t197 = t201 + t252;
t196 = -t190 * qJD(2) + qJD(4) * t134 - t206;
t193 = rSges(4,1) * t134 - rSges(4,2) * t132;
t195 = -t193 * qJD(2) - t206;
t191 = rSges(6,1) * t133 - rSges(6,2) * t131;
t189 = -t131 * t51 + t133 * t53;
t188 = -t131 * t52 + t133 * t54;
t176 = -Icges(6,2) * t131 + t224;
t82 = -Icges(6,6) * t134 + t176 * t132;
t179 = Icges(6,1) * t133 - t225;
t83 = -Icges(6,5) * t134 + t179 * t132;
t187 = t131 * t82 - t133 * t83;
t21 = t189 * t132 - t134 * t49;
t22 = t188 * t132 - t134 * t50;
t184 = t21 * t140 + t22 * t142;
t55 = rSges(6,1) * t105 + rSges(6,2) * t104 + rSges(6,3) * t222;
t56 = rSges(6,1) * t107 + rSges(6,2) * t106 + rSges(6,3) * t221;
t183 = -t140 * t56 + t142 * t55;
t85 = -rSges(6,3) * t134 + t191 * t132;
t182 = t201 + t249 - t85;
t169 = t196 - (rSges(5,3) * t132 + t192 * t134) * qJD(2);
t67 = -t105 * qJD(5) + t131 * t205;
t68 = t104 * qJD(5) - t133 * t205;
t34 = Icges(6,5) * t68 + Icges(6,6) * t67 + Icges(6,3) * t203;
t168 = t132 * t34 + t49 * t211;
t69 = -t107 * qJD(5) + t131 * t204;
t70 = t106 * qJD(5) - t133 * t204;
t35 = Icges(6,5) * t70 + Icges(6,6) * t69 + Icges(6,3) * t202;
t167 = t132 * t35 + t50 * t211;
t48 = (-rSges(6,1) * t131 - rSges(6,2) * t133) * t208 + (rSges(6,3) * t132 + t191 * t134) * qJD(2);
t166 = t196 - t48 - (-t214 * t132 - t238 * t134) * qJD(2);
t42 = -t126 * t259 - t213;
t165 = t42 * t193;
t150 = t129 * t134 - t132 * t144 - t190;
t149 = qJD(2) * (Icges(5,5) * t134 + (-Icges(5,1) * t141 + Icges(5,4) * t139) * t132);
t148 = qJD(2) * (Icges(5,6) * t134 + (-Icges(5,4) * t141 + Icges(5,2) * t139) * t132);
t115 = t134 * t215 + t218;
t114 = -t134 * t217 + t216;
t113 = t134 * t216 - t217;
t112 = -t134 * t218 - t215;
t87 = t195 * t142;
t86 = t195 * t140;
t76 = t142 * t149;
t75 = t140 * t149;
t74 = t142 * t148;
t73 = t140 * t148;
t58 = t197 * t142;
t57 = t197 * t140;
t47 = (-Icges(6,1) * t131 - t224) * t208 + (Icges(6,5) * t132 + t179 * t134) * qJD(2);
t46 = (-Icges(6,2) * t133 - t225) * t208 + (Icges(6,6) * t132 + t176 * t134) * qJD(2);
t44 = t169 * t142;
t43 = t169 * t140;
t41 = rSges(6,1) * t70 + rSges(6,2) * t69 + rSges(6,3) * t202;
t40 = rSges(6,1) * t68 + rSges(6,2) * t67 + rSges(6,3) * t203;
t39 = Icges(6,1) * t70 + Icges(6,4) * t69 + Icges(6,5) * t202;
t38 = Icges(6,1) * t68 + Icges(6,4) * t67 + Icges(6,5) * t203;
t37 = Icges(6,4) * t70 + Icges(6,2) * t69 + Icges(6,6) * t202;
t36 = Icges(6,4) * t68 + Icges(6,2) * t67 + Icges(6,6) * t203;
t33 = t182 * t142;
t32 = t182 * t140;
t31 = -t134 * t56 - t85 * t221;
t30 = t134 * t55 + t85 * t222;
t29 = -t187 * t132 - t232;
t28 = t183 * t132;
t27 = t166 * t142;
t26 = t166 * t140;
t25 = t106 * t82 + t107 * t83 + t81 * t221;
t24 = t104 * t82 + t105 * t83 + t81 * t222;
t23 = t252 * t259 + t198;
t20 = t140 * (rSges(5,1) * t113 + rSges(5,2) * t112 + rSges(5,3) * t222) + t142 * (rSges(5,1) * t115 + rSges(5,2) * t114 + rSges(5,3) * t221) + t199;
t19 = t106 * t52 + t107 * t54 + t50 * t221;
t16 = t104 * t51 + t105 * t53 + t49 * t222;
t15 = -t48 * t221 - t134 * t41 + (t132 * t56 - t85 * t219) * qJD(2);
t14 = t48 * t222 + t134 * t40 + (-t132 * t55 + t85 * t220) * qJD(2);
t13 = (t150 * t142 + t56) * t142 + (t150 * t140 + t55) * t140 + t199;
t12 = (-t140 * t41 + t142 * t40) * t132 + t183 * t211;
t11 = t140 * t40 + t142 * t41 + t249 * t259 + t198;
t10 = t106 * t37 + t107 * t39 + t167 * t142 + t52 * t69 + t54 * t70;
t9 = t106 * t36 + t107 * t38 + t168 * t142 + t51 * t69 + t53 * t70;
t8 = t104 * t37 + t105 * t39 + t167 * t140 + t52 * t67 + t54 * t68;
t7 = t104 * t36 + t105 * t38 + t168 * t140 + t51 * t67 + t53 * t68;
t6 = (t188 * qJD(2) - t35) * t134 + (qJD(2) * t50 - t131 * t37 + t133 * t39 + (-t131 * t54 - t133 * t52) * qJD(5)) * t132;
t5 = (t189 * qJD(2) - t34) * t134 + (qJD(2) * t49 - t131 * t36 + t133 * t38 + (-t131 * t53 - t133 * t51) * qJD(5)) * t132;
t4 = t10 * t140 - t142 * t9;
t3 = t140 * t8 - t142 * t7;
t2 = -(t106 * t46 + t107 * t47 + t69 * t82 + t70 * t83) * t134 + (t9 * t140 + (t10 - t233) * t142) * t132 + (t25 * t132 + (t231 + (t19 - t232) * t142) * t134) * qJD(2);
t1 = -(t104 * t46 + t105 * t47 + t67 * t82 + t68 * t83) * t134 + (t8 * t142 + (t7 - t233) * t140) * t132 + (t24 * t132 + (t230 + (t16 - t232) * t140) * t134) * qJD(2);
t45 = [0; m(4) * t42 + m(5) * t23 + m(6) * t11 - t266; -t142 * t3 + (t11 * t13 + t26 * t32 + t27 * t33) * t245 + 0.2e1 * m(5) * (t20 * t23 + t43 * t57 + t44 * t58) + (t236 * t42 + (t142 * t165 + t200 * t87) * t142) * t262 + (t200 * t86 * t262 + t4 + (t114 * t74 + t115 * t76 + t165 * t262) * t140 + t260 * t135 + (-t112 * t74 - t113 * t76 - t114 * t73 - t115 * t75 - t261 * t140 + t260 * t142) * t142) * t140 + (t112 * t73 + t113 * t75 - t261 * t142) * t136 + 0.2e1 * t265 * (rSges(3,1) * t146 - rSges(3,2) * t145) * t266; 0; m(6) * (t140 * t27 - t142 * t26) + m(5) * (t140 * t44 - t142 * t43) + m(4) * (t140 * t87 - t142 * t86); 0; (m(5) + m(6)) * t212; m(6) * (-t11 * t134 + t27 * t221 + t26 * t222) + m(5) * (-t134 * t23 + t44 * t221 + t43 * t222) + 0.2e1 * ((t13 * t132 + t33 * t219 + t32 * t220) * t243 + (t132 * t20 + t58 * t219 + t57 * t220) * t244) * qJD(2); 0; -0.4e1 * (t244 + t243) * t265 * t132 * t211; m(6) * t12; t2 * t242 - t142 * t1 / 0.2e1 - t134 * (t6 * t140 - t5 * t142) / 0.2e1 + m(6) * (t11 * t28 + t12 * t13 + t14 * t33 + t15 * t32 + t26 * t31 + t27 * t30) + (t4 * t241 + t3 * t242) * t132 + (t132 * (t140 * t22 - t142 * t21) / 0.2e1 + ((t19 * t140 - t18 * t142) * t241 + (t17 * t140 - t16 * t142) * t242) * t134) * qJD(2); m(6) * (t14 * t140 - t142 * t15); m(6) * (-t12 * t134 + (t14 * t142 + t140 * t15) * t132 + (t132 * t28 + (t140 * t31 + t142 * t30) * t134) * qJD(2)); (t12 * t28 + t14 * t30 + t15 * t31) * t245 + (-t134 * t25 + (t142 * t19 + t231) * t132) * t202 + t2 * t221 + (-t134 * t24 + (t140 * t16 + t230) * t132) * t203 + t1 * t222 + (t184 * t132 - t29 * t134) * t212 - t134 * ((t233 + (t187 * t134 + t184) * qJD(2)) * t134 + (t6 * t142 + t5 * t140 - (qJD(2) * t81 - t131 * t46 + t133 * t47 + (-t131 * t83 - t133 * t82) * qJD(5)) * t134 + t29 * qJD(2)) * t132);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t45(1), t45(2), t45(4), t45(7), t45(11); t45(2), t45(3), t45(5), t45(8), t45(12); t45(4), t45(5), t45(6), t45(9), t45(13); t45(7), t45(8), t45(9), t45(10), t45(14); t45(11), t45(12), t45(13), t45(14), t45(15);];
Mq = res;
