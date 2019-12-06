% Calculate time derivative of joint inertia matrix for
% S5RPPRR3
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:38
% EndTime: 2019-12-05 17:41:44
% DurationCPUTime: 2.33s
% Computational Cost: add. (6403->304), mult. (4524->433), div. (0->0), fcn. (3474->10), ass. (0->178)
t121 = pkin(9) + qJ(4);
t114 = sin(t121);
t123 = qJ(1) + pkin(8);
t115 = sin(t123);
t116 = cos(t121);
t179 = qJD(4) * t116;
t117 = cos(t123);
t182 = qJD(1) * t117;
t239 = t114 * t182 + t115 * t179;
t118 = qJ(5) + t121;
t108 = sin(t118);
t109 = cos(t118);
t86 = rSges(6,1) * t108 + rSges(6,2) * t109;
t178 = qJD(4) * t117;
t93 = rSges(5,1) * t114 + rSges(5,2) * t116;
t238 = t93 * t178;
t139 = Icges(6,5) * t109 - Icges(6,6) * t108;
t60 = Icges(6,3) * t115 + t117 * t139;
t237 = qJD(1) * t60;
t203 = Icges(5,4) * t116;
t144 = -Icges(5,2) * t114 + t203;
t69 = Icges(5,6) * t117 - t115 * t144;
t204 = Icges(5,4) * t114;
t147 = Icges(5,1) * t116 - t204;
t71 = Icges(5,5) * t117 - t115 * t147;
t151 = t114 * t69 - t116 * t71;
t236 = t117 * t151;
t201 = Icges(6,4) * t109;
t142 = -Icges(6,2) * t108 + t201;
t61 = Icges(6,6) * t117 - t115 * t142;
t202 = Icges(6,4) * t108;
t145 = Icges(6,1) * t109 - t202;
t63 = Icges(6,5) * t117 - t115 * t145;
t156 = t108 * t61 - t109 * t63;
t235 = t117 * t156;
t107 = t117 * rSges(5,3);
t220 = rSges(5,2) * t114;
t185 = t115 * t220 + t107;
t126 = -pkin(6) - qJ(3);
t125 = cos(pkin(9));
t110 = pkin(3) * t125 + pkin(2);
t223 = rSges(5,1) * t116;
t162 = -t110 - t223;
t226 = sin(qJ(1)) * pkin(1);
t50 = t115 * t162 - t117 * t126 + t185 - t226;
t104 = t115 * t126;
t158 = -t220 + t223;
t215 = t115 * rSges(5,3);
t225 = cos(qJ(1)) * pkin(1);
t51 = -t215 - t225 + t104 + (-t110 - t158) * t117;
t234 = -t115 * t51 + t117 * t50;
t141 = Icges(5,5) * t116 - Icges(5,6) * t114;
t68 = Icges(5,3) * t115 + t117 * t141;
t62 = Icges(6,6) * t115 + t117 * t142;
t70 = Icges(5,6) * t115 + t117 * t144;
t64 = Icges(6,5) * t115 + t117 * t145;
t72 = Icges(5,5) * t115 + t117 * t147;
t138 = rSges(4,1) * t125 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t209 = rSges(4,3) + qJ(3);
t233 = t115 * t138 - t117 * t209;
t122 = qJD(4) + qJD(5);
t76 = t142 * t122;
t77 = t145 * t122;
t83 = Icges(6,5) * t108 + Icges(6,6) * t109;
t84 = Icges(6,2) * t109 + t202;
t85 = Icges(6,1) * t108 + t201;
t232 = (t122 * t85 + t76) * t108 + (t122 * t84 - t77) * t109 - qJD(1) * t83;
t231 = 2 * m(5);
t230 = 2 * m(6);
t112 = t115 ^ 2;
t113 = t117 ^ 2;
t229 = m(5) * t93;
t228 = t115 / 0.2e1;
t227 = t117 / 0.2e1;
t120 = -pkin(7) + t126;
t100 = pkin(4) * t116 + t110;
t186 = t100 - t110;
t106 = t117 * rSges(6,3);
t219 = rSges(6,2) * t108;
t205 = t115 * t219 + t106;
t221 = rSges(6,1) * t109;
t65 = -t115 * t221 + t205;
t224 = -(-t120 + t126) * t117 + t186 * t115 - t65;
t217 = t114 * t71;
t216 = t114 * t72;
t214 = t115 * rSges(6,3);
t212 = t116 * t69;
t211 = t116 * t70;
t208 = -rSges(6,3) + t120;
t180 = qJD(4) * t115;
t173 = t114 * t180;
t183 = qJD(1) * t115;
t207 = -pkin(4) * t173 - t120 * t183;
t206 = t110 * t183 + t126 * t182;
t59 = Icges(6,3) * t117 - t115 * t139;
t194 = qJD(1) * t59;
t67 = Icges(5,3) * t117 - t115 * t141;
t193 = qJD(1) * t67;
t191 = t108 * t122;
t190 = t109 * t122;
t189 = t115 * t122;
t188 = t117 * t120;
t187 = t117 * t122;
t184 = t112 + t113;
t181 = qJD(4) * t114;
t177 = t182 * t219 + t189 * t86;
t176 = -rSges(5,1) * t173 - rSges(5,2) * t239;
t175 = pkin(4) * t181;
t171 = pkin(4) * t114 + t86;
t35 = qJD(1) * t61 - t187 * t84;
t168 = t122 * t64 + t35;
t36 = -qJD(1) * t62 + t84 * t189;
t167 = t122 * t63 + t36;
t37 = qJD(1) * t63 - t187 * t85;
t166 = -t122 * t62 + t37;
t38 = -qJD(1) * t64 + t85 * t189;
t165 = t122 * t61 - t38;
t161 = -t100 - t221;
t160 = qJD(1) * t226 - qJD(3) * t115;
t157 = -t219 + t221;
t155 = t108 * t62 - t109 * t64;
t154 = t108 * t84 - t109 * t85;
t150 = t114 * t70 - t116 * t72;
t40 = t115 * t161 - t188 + t205 - t226;
t136 = t100 + t157;
t41 = t115 * t208 - t117 * t136 - t225;
t149 = t115 * t41 - t117 * t40;
t57 = t171 * t115;
t58 = t171 * t117;
t148 = t115 * t57 + t117 * t58;
t146 = Icges(5,1) * t114 + t203;
t143 = Icges(5,2) * t116 + t204;
t140 = Icges(5,5) * t114 + Icges(5,6) * t116;
t13 = t115 * t156 + t117 * t59;
t134 = t155 * t115;
t14 = t117 * t60 + t134;
t15 = t115 * t59 - t235;
t16 = t60 * t115 - t117 * t155;
t33 = -t187 * t83 + t194;
t34 = t83 * t189 - t237;
t137 = -t183 * (t115 * t14 + t117 * t13) + t115 * ((t115 * t33 + (-t15 + t134) * qJD(1)) * t115 + (t16 * qJD(1) + (-t108 * t36 + t109 * t38 - t190 * t61 - t191 * t63 + t194) * t117 + (t34 + t166 * t109 - t168 * t108 + (t156 + t60) * qJD(1)) * t115) * t117) + t117 * ((t117 * t34 + (t14 + t235) * qJD(1)) * t117 + (-t13 * qJD(1) + (t108 * t35 - t109 * t37 + t190 * t62 + t191 * t64 - t237) * t115 + (t33 + t165 * t109 + t167 * t108 + (t155 - t59) * qJD(1)) * t117) * t115) + (t115 * t16 + t117 * t15) * t182;
t129 = qJD(1) * t154 + t122 * t139;
t135 = (t108 * t166 + t109 * t168 + t129 * t115 - t117 * t232) * t228 + (-t108 * t165 + t109 * t167 + t115 * t232 + t129 * t117) * t227 - (t108 * t63 + t109 * t61 + t115 * t154 + t117 * t83) * t183 / 0.2e1 + (t108 * t64 + t109 * t62 + t115 * t83 - t117 * t154) * t182 / 0.2e1;
t133 = t150 * t115;
t132 = qJD(4) * t146;
t131 = qJD(4) * t143;
t55 = -t115 * t209 - t117 * t138 - t225;
t105 = qJD(3) * t117;
t82 = t158 * qJD(4);
t78 = t157 * t122;
t74 = t117 * t158 + t215;
t73 = -t115 * t223 + t185;
t66 = t117 * t157 + t214;
t56 = t117 * t66;
t54 = -t226 - t233;
t53 = -t115 * t120 + t117 * t186 + t104;
t49 = qJD(1) * t55 + t105;
t48 = qJD(1) * t233 + t160;
t43 = -qJD(1) * t68 + t140 * t180;
t42 = -t140 * t178 + t193;
t39 = (-t117 * t221 - t214) * qJD(1) + t177;
t32 = t117 * (-t86 * t187 + (-t115 * t157 + t106) * qJD(1));
t31 = -t115 * t65 + t56;
t30 = pkin(4) * t239 + t115 * t78 + t86 * t182;
t29 = t86 * t183 - t117 * t78 + (t114 * t183 - t116 * t178) * pkin(4);
t24 = t105 + (-t225 + t162 * t117 + (-rSges(5,3) + t126) * t115) * qJD(1) - t176;
t23 = t238 + (t115 * t158 - t107) * qJD(1) + t160 + t206;
t22 = t115 * t68 - t117 * t150;
t21 = t115 * t67 - t236;
t20 = t117 * t68 + t133;
t19 = t115 * t151 + t117 * t67;
t18 = t105 + (t117 * t161 - t214 - t225) * qJD(1) + t177 - t207;
t17 = (t122 * t86 + t175) * t117 + (t115 * t136 + t117 * t208) * qJD(1) + t160;
t12 = t115 * t224 + t117 * t53 + t56;
t11 = ((-t74 + t215) * qJD(1) + t176) * t115 + (-t238 + (-t73 + t185) * qJD(1)) * t117;
t10 = -t115 * t39 + t32 + (-t115 * t66 - t117 * t65) * qJD(1);
t3 = t32 + (-t117 * t175 + (-t188 + t224) * qJD(1) + t206) * t117 + (-t39 + (-t110 * t117 + t104 - t53 - t66) * qJD(1) + t207) * t115;
t1 = [0.2e1 * m(4) * (t48 * t55 + t49 * t54) + t116 * t132 + t147 * t181 - t114 * t131 + t144 * t179 + (t23 * t51 + t24 * t50) * t231 + t85 * t190 + t108 * t77 - t84 * t191 + t109 * t76 + (t17 * t41 + t18 * t40) * t230; 0; 0; m(4) * (t115 * t49 + t117 * t48 + (-t115 * t55 + t117 * t54) * qJD(1)) + m(5) * (qJD(1) * t234 + t115 * t24 + t117 * t23) + m(6) * (-qJD(1) * t149 + t115 * t18 + t117 * t17); 0; 0; (-qJD(4) * t151 + t114 * (-qJD(1) * t72 + t115 * t132) + t116 * (-qJD(1) * t70 + t115 * t131)) * t227 + (-qJD(4) * t150 + t114 * (qJD(1) * t71 - t146 * t178) + t116 * (qJD(1) * t69 - t143 * t178)) * t228 + m(5) * ((t115 * t23 - t117 * t24) * t93 - t234 * t82) + m(6) * (t17 * t57 - t18 * t58 + t29 * t40 + t30 * t41) + (t113 / 0.2e1 + t112 / 0.2e1) * t141 * qJD(4) + ((-t217 / 0.2e1 - t212 / 0.2e1 + t50 * t229) * t115 + (t216 / 0.2e1 + t211 / 0.2e1 + t51 * t229) * t117) * qJD(1) + t135; m(5) * t11 + m(6) * t3; m(6) * (-qJD(1) * t148 + t115 * t29 + t117 * t30); ((-t115 * t73 + t117 * t74) * t11 + t184 * t93 * t82) * t231 - (t115 * t20 + t117 * t19) * t183 + t117 * ((t117 * t43 + (t20 + t236) * qJD(1)) * t117 + (-t19 * qJD(1) + (t179 * t70 + t181 * t72) * t115 + (t42 + (t212 + t217) * qJD(4) + (t150 - t67) * qJD(1)) * t117) * t115) + (t115 * t22 + t117 * t21) * t182 + t115 * ((t115 * t42 + (-t21 + t133) * qJD(1)) * t115 + (t22 * qJD(1) + (-t179 * t69 - t181 * t71 + t193) * t117 + (t43 + (-t211 - t216) * qJD(4) + t151 * qJD(1)) * t115) * t117) + (t12 * t3 - t29 * t58 + t30 * t57) * t230 + t137; m(6) * (t149 * t78 + (t115 * t17 - t117 * t18 + (t115 * t40 + t117 * t41) * qJD(1)) * t86) + t135; m(6) * t10; 0; m(6) * (t10 * t12 + t31 * t3 + t148 * t78 + (t115 * t30 - t117 * t29 + (-t115 * t58 + t117 * t57) * qJD(1)) * t86) + t137; (t184 * t78 * t86 + t31 * t10) * t230 + t137;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
