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
% m [6x1]
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:14:08
% EndTime: 2022-01-23 09:14:14
% DurationCPUTime: 2.77s
% Computational Cost: add. (6403->300), mult. (4524->428), div. (0->0), fcn. (3474->10), ass. (0->175)
t123 = qJ(1) + pkin(8);
t114 = sin(t123);
t116 = cos(t123);
t121 = pkin(9) + qJ(4);
t117 = qJ(5) + t121;
t108 = sin(t117);
t212 = rSges(6,2) * t108;
t109 = cos(t117);
t214 = rSges(6,1) * t109;
t161 = -t212 + t214;
t68 = t114 * rSges(6,3) + t161 * t116;
t113 = sin(t121);
t115 = cos(t121);
t94 = rSges(5,1) * t113 + rSges(5,2) * t115;
t138 = qJD(4) * t94;
t231 = t114 * t138;
t122 = qJD(4) + qJD(5);
t211 = rSges(6,2) * t109;
t89 = rSges(6,1) * t108 + t211;
t139 = t89 * t122;
t177 = qJD(4) * t113;
t173 = pkin(4) * t177;
t229 = t139 + t173;
t144 = Icges(6,5) * t109 - Icges(6,6) * t108;
t61 = -Icges(6,3) * t116 + t144 * t114;
t228 = qJD(1) * t61;
t196 = Icges(5,4) * t115;
t148 = -Icges(5,2) * t113 + t196;
t72 = Icges(5,6) * t114 + t148 * t116;
t197 = Icges(5,4) * t113;
t151 = Icges(5,1) * t115 - t197;
t74 = Icges(5,5) * t114 + t151 * t116;
t154 = t113 * t72 - t115 * t74;
t227 = t154 * t114;
t71 = -Icges(5,6) * t116 + t148 * t114;
t73 = -Icges(5,5) * t116 + t151 * t114;
t155 = t113 * t71 - t115 * t73;
t226 = t155 * t116;
t194 = Icges(6,4) * t109;
t146 = -Icges(6,2) * t108 + t194;
t64 = Icges(6,6) * t114 + t146 * t116;
t195 = Icges(6,4) * t108;
t149 = Icges(6,1) * t109 - t195;
t66 = Icges(6,5) * t114 + t149 * t116;
t159 = t108 * t64 - t109 * t66;
t225 = t159 * t114;
t63 = -Icges(6,6) * t116 + t146 * t114;
t65 = -Icges(6,5) * t116 + t149 * t114;
t160 = t108 * t63 - t109 * t65;
t224 = t160 * t116;
t126 = -pkin(6) - qJ(3);
t125 = cos(pkin(9));
t110 = t125 * pkin(3) + pkin(2);
t213 = rSges(5,2) * t113;
t215 = rSges(5,1) * t115;
t162 = -t213 + t215;
t141 = -t110 - t162;
t216 = sin(qJ(1)) * pkin(1);
t50 = -t216 + (rSges(5,3) - t126) * t116 + t141 * t114;
t119 = cos(qJ(1)) * pkin(1);
t183 = t114 * t126;
t107 = t114 * rSges(5,3);
t198 = t116 * t215 + t107;
t51 = -t183 + t119 + (t110 - t213) * t116 + t198;
t223 = t114 * t51 + t116 * t50;
t87 = Icges(6,2) * t109 + t195;
t88 = Icges(6,1) * t108 + t194;
t158 = t108 * t87 - t109 * t88;
t222 = t158 * qJD(1) + t144 * t122;
t145 = Icges(5,5) * t115 - Icges(5,6) * t113;
t69 = -Icges(5,3) * t116 + t145 * t114;
t143 = rSges(4,1) * t125 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t202 = rSges(4,3) + qJ(3);
t55 = t202 * t114 + t143 * t116 + t119;
t221 = 2 * m(5);
t220 = 2 * m(6);
t111 = t114 ^ 2;
t112 = t116 ^ 2;
t219 = m(5) * t94;
t218 = t114 / 0.2e1;
t217 = -t116 / 0.2e1;
t67 = -rSges(6,3) * t116 + t161 * t114;
t31 = t114 * t67 + t116 * t68;
t210 = t113 * t73;
t209 = t113 * t74;
t207 = t115 * t71;
t206 = t115 * t72;
t204 = t122 * t87;
t203 = t122 * t88;
t120 = pkin(7) - t126;
t201 = rSges(6,3) + t120;
t101 = t120 * t114;
t97 = pkin(4) * t115 + t110;
t200 = t116 * t97 + t101;
t178 = qJD(1) * t116;
t179 = qJD(1) * t114;
t199 = rSges(6,3) * t178 + t179 * t212;
t62 = Icges(6,3) * t114 + t144 * t116;
t187 = qJD(1) * t62;
t70 = Icges(5,3) * t114 + t145 * t116;
t186 = qJD(1) * t70;
t185 = t108 * t122;
t184 = t109 * t122;
t182 = t116 * t122;
t181 = t116 * t126;
t180 = t111 + t112;
t176 = qJD(4) * t115;
t175 = t114 * (t68 * qJD(1) - t114 * t139) + t116 * (-t182 * t211 + (-t108 * t182 - t109 * t179) * rSges(6,1) + t199) + t67 * t178;
t174 = t116 * t213;
t172 = t113 * t179;
t171 = -pkin(4) * t113 - t89;
t36 = -qJD(1) * t63 - t116 * t204;
t168 = t122 * t66 + t36;
t37 = t64 * qJD(1) - t114 * t204;
t167 = t122 * t65 + t37;
t38 = -qJD(1) * t65 - t116 * t203;
t166 = -t122 * t64 + t38;
t39 = t66 * qJD(1) - t114 * t203;
t165 = t122 * t63 - t39;
t13 = -t160 * t114 - t61 * t116;
t14 = -t116 * t62 - t225;
t15 = t114 * t61 - t224;
t16 = t114 * t62 - t159 * t116;
t86 = Icges(6,5) * t108 + Icges(6,6) * t109;
t135 = t86 * t122;
t34 = -t116 * t135 - t228;
t35 = -t114 * t135 + t187;
t164 = -t116 * ((t116 * t35 + (t14 + t224) * qJD(1)) * t116 + (t13 * qJD(1) + (-t108 * t36 + t109 * t38 - t64 * t184 - t66 * t185 + t187) * t114 + (-t34 + t165 * t109 + t167 * t108 + (-t159 - t61) * qJD(1)) * t116) * t114) + t114 * ((t114 * t34 + (t15 + t225) * qJD(1)) * t114 + (t16 * qJD(1) + (t108 * t37 - t109 * t39 + t63 * t184 + t65 * t185 - t228) * t116 + (-t35 + t166 * t109 - t168 * t108 + (-t160 + t62) * qJD(1)) * t114) * t116) + (t114 * t14 - t116 * t13) * t179 + (t114 * t16 - t116 * t15) * t178;
t142 = -t161 - t97;
t40 = t142 * t114 + t201 * t116 - t216;
t41 = t119 + t68 + t200;
t153 = t114 * t41 + t116 * t40;
t59 = t171 * t114;
t60 = t171 * t116;
t152 = t114 * t59 + t116 * t60;
t150 = Icges(5,1) * t113 + t196;
t147 = Icges(5,2) * t115 + t197;
t78 = t146 * t122;
t79 = t149 * t122;
t129 = qJD(1) * t86 + (t79 - t204) * t109 + (-t78 - t203) * t108;
t140 = (t166 * t108 + t168 * t109 + t114 * t222 + t129 * t116) * t218 + (-t165 * t108 + t167 * t109 + t129 * t114 - t116 * t222) * t217 + (t108 * t65 + t109 * t63 - t158 * t114 - t116 * t86) * t179 / 0.2e1 + (t108 * t66 + t109 * t64 + t114 * t86 - t158 * t116) * t178 / 0.2e1;
t134 = qJD(4) * t150;
t133 = qJD(4) * t147;
t132 = qJD(4) * (-Icges(5,5) * t113 - Icges(5,6) * t115);
t130 = rSges(5,2) * t172 + rSges(5,3) * t178 - t116 * t138;
t54 = -t143 * t114 + t202 * t116 - t216;
t105 = qJD(3) * t116;
t104 = qJD(3) * t114;
t100 = t126 * t179;
t99 = t120 * t178;
t85 = t162 * qJD(4);
t80 = t161 * t122;
t76 = -t174 + t198;
t75 = -rSges(5,3) * t116 + t162 * t114;
t53 = -t110 * t116 + t183 + t200;
t52 = (-t120 - t126) * t116 + t114 * (-t110 + t97);
t49 = -qJD(1) * t55 + t105;
t48 = t54 * qJD(1) + t104;
t43 = t114 * t132 + t186;
t42 = -qJD(1) * t69 + t116 * t132;
t30 = -t89 * t178 - t114 * t80 + (-t113 * t178 - t114 * t176) * pkin(4);
t29 = t89 * t179 - t116 * t80 + (-t116 * t176 + t172) * pkin(4);
t24 = t100 + t105 + t231 + (t141 * t116 - t107 - t119) * qJD(1);
t23 = t104 + (-t216 - t181 + (-t110 - t215) * t114) * qJD(1) + t130;
t22 = t70 * t114 - t154 * t116;
t21 = t114 * t69 - t226;
t20 = -t116 * t70 - t227;
t19 = -t155 * t114 - t116 * t69;
t18 = t105 + t229 * t114 + (-t201 * t114 + t142 * t116 - t119) * qJD(1);
t17 = t104 + t99 - t229 * t116 + (-t216 + (-t97 - t214) * t114) * qJD(1) + t199;
t12 = t114 * t52 + t116 * t53 + t31;
t11 = (qJD(1) * t75 + t130) * t116 + (-t231 + (-t174 - t76 + t107) * qJD(1)) * t114;
t10 = -t68 * t179 + t175;
t3 = t114 * (-t114 * t173 + t100) + t116 * (-t116 * t173 + t99) + ((t52 + t181) * t116 + (-t53 - t68 + t101) * t114) * qJD(1) + t175;
t1 = [(t17 * t41 + t18 * t40) * t220 + t88 * t184 + t108 * t79 - t87 * t185 + t109 * t78 + (t23 * t51 + t24 * t50) * t221 + 0.2e1 * m(4) * (t48 * t55 + t49 * t54) + (t151 - t147) * t177 + (t150 + t148) * t176; 0; 0; m(6) * (t153 * qJD(1) + t114 * t18 - t116 * t17) + m(5) * (t223 * qJD(1) + t114 * t24 - t116 * t23) + m(4) * (t114 * t49 - t116 * t48 + (t114 * t55 + t116 * t54) * qJD(1)); 0; 0; (-t154 * qJD(4) + t113 * (-qJD(1) * t73 - t116 * t134) + t115 * (-qJD(1) * t71 - t116 * t133)) * t218 + (-t155 * qJD(4) + t113 * (t74 * qJD(1) - t114 * t134) + t115 * (t72 * qJD(1) - t114 * t133)) * t217 + m(6) * (t59 * t17 + t18 * t60 + t29 * t40 + t30 * t41) + m(5) * ((-t114 * t23 - t116 * t24) * t94 - t223 * t85) + (t111 / 0.2e1 + t112 / 0.2e1) * t145 * qJD(4) + ((-t51 * t219 + t209 / 0.2e1 + t206 / 0.2e1) * t116 + (t50 * t219 + t210 / 0.2e1 + t207 / 0.2e1) * t114) * qJD(1) + t140; m(5) * t11 + m(6) * t3; m(6) * (t152 * qJD(1) + t114 * t29 - t116 * t30); ((t114 * t75 + t116 * t76) * t11 + t180 * t94 * t85) * t221 + (t114 * t22 - t116 * t21) * t178 + t114 * ((t114 * t42 + (t21 + t227) * qJD(1)) * t114 + (t22 * qJD(1) + (t71 * t176 + t73 * t177) * t116 + (-t43 + (-t206 - t209) * qJD(4) + (-t155 + t70) * qJD(1)) * t114) * t116) + (t114 * t20 - t116 * t19) * t179 - t116 * ((t116 * t43 + (t20 + t226) * qJD(1)) * t116 + (t19 * qJD(1) + (-t72 * t176 - t74 * t177 + t186) * t114 + (-t42 + (t207 + t210) * qJD(4) - t154 * qJD(1)) * t116) * t114) + (t12 * t3 + t29 * t60 + t59 * t30) * t220 + t164; m(6) * (-t153 * t80 + (-t114 * t17 - t116 * t18 + (t114 * t40 - t116 * t41) * qJD(1)) * t89) + t140; m(6) * t10; 0; m(6) * (t10 * t12 + t31 * t3 - t152 * t80 + (-t114 * t30 - t116 * t29 + (t114 * t60 - t116 * t59) * qJD(1)) * t89) + t164; (t180 * t89 * t80 + t31 * t10) * t220 + t164;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
