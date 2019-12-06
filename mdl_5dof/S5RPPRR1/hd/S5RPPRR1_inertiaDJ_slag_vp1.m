% Calculate time derivative of joint inertia matrix for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:08
% DurationCPUTime: 3.34s
% Computational Cost: add. (3040->326), mult. (4660->469), div. (0->0), fcn. (3536->6), ass. (0->187)
t120 = sin(qJ(1));
t122 = cos(qJ(1));
t118 = qJ(4) + qJ(5);
t109 = sin(t118);
t110 = cos(t118);
t159 = rSges(6,1) * t109 + rSges(6,2) * t110;
t63 = -rSges(6,3) * t120 + t159 * t122;
t115 = qJD(4) + qJD(5);
t85 = rSges(6,1) * t110 - rSges(6,2) * t109;
t242 = t85 * t115;
t144 = Icges(6,5) * t109 + Icges(6,6) * t110;
t56 = Icges(6,3) * t122 + t144 * t120;
t241 = qJD(1) * t56;
t210 = t122 * rSges(6,3);
t62 = t159 * t120 + t210;
t240 = qJD(1) * t62;
t119 = sin(qJ(4));
t121 = cos(qJ(4));
t145 = Icges(5,5) * t119 + Icges(5,6) * t121;
t70 = Icges(5,3) * t122 + t145 * t120;
t239 = qJD(1) * t70;
t116 = t120 ^ 2;
t117 = t122 ^ 2;
t184 = t116 + t117;
t113 = t122 * qJ(2);
t220 = rSges(5,2) * t121;
t160 = rSges(5,1) * t119 + t220;
t226 = -pkin(1) - qJ(3);
t136 = -t160 + t226;
t229 = -rSges(5,3) - pkin(6);
t125 = t136 * t120 + t229 * t122;
t49 = t113 + t125;
t185 = t122 * pkin(1) + t120 * qJ(2);
t173 = t122 * qJ(3) + t185;
t189 = t119 * t122;
t188 = rSges(5,1) * t189 + t122 * t220;
t50 = t229 * t120 + t173 + t188;
t238 = t120 * t50 + t122 * t49;
t178 = qJD(4) * t122;
t172 = t119 * t178;
t182 = qJD(1) * t122;
t186 = qJ(2) * t182 + qJD(2) * t120;
t174 = qJD(3) * t122 + t186;
t171 = t121 * t178;
t97 = rSges(5,1) * t171;
t22 = -rSges(5,2) * t172 + t125 * qJD(1) + t174 + t97;
t183 = qJD(1) * t120;
t105 = pkin(6) * t183;
t108 = qJD(2) * t122;
t221 = rSges(5,2) * t119;
t94 = rSges(5,1) * t121 - t221;
t23 = t105 + t108 + (-t94 * qJD(4) - qJD(3)) * t120 + ((rSges(5,3) - qJ(2)) * t120 + t136 * t122) * qJD(1);
t237 = t120 * t22 + t122 * t23;
t201 = Icges(6,4) * t110;
t83 = -Icges(6,2) * t109 + t201;
t202 = Icges(6,4) * t109;
t84 = Icges(6,1) * t110 - t202;
t156 = t109 * t84 + t110 * t83;
t236 = t156 * qJD(1) - t144 * t115;
t146 = Icges(6,2) * t110 + t202;
t58 = Icges(6,6) * t122 + t146 * t120;
t204 = Icges(5,4) * t119;
t147 = Icges(5,2) * t121 + t204;
t72 = Icges(5,6) * t122 + t147 * t120;
t148 = Icges(6,1) * t109 + t201;
t60 = Icges(6,5) * t122 + t148 * t120;
t203 = Icges(5,4) * t121;
t149 = Icges(5,1) * t119 + t203;
t74 = Icges(5,5) * t122 + t149 * t120;
t235 = 2 * m(5);
t234 = 2 * m(6);
t233 = m(5) * t94;
t232 = -t120 / 0.2e1;
t231 = t122 / 0.2e1;
t230 = rSges(3,2) - pkin(1);
t228 = pkin(4) * t119;
t59 = -Icges(6,6) * t120 + t146 * t122;
t61 = -Icges(6,5) * t120 + t148 * t122;
t157 = t109 * t61 + t110 * t59;
t139 = t157 * t120;
t158 = t109 * t60 + t110 * t58;
t140 = t158 * t122;
t16 = -t120 * t56 + t140;
t135 = t115 * t84;
t35 = -qJD(1) * t60 + t122 * t135;
t165 = t115 * t59 - t35;
t134 = t115 * t83;
t33 = -qJD(1) * t58 + t122 * t134;
t167 = -t115 * t61 - t33;
t57 = -Icges(6,3) * t120 + t144 * t122;
t17 = -t120 * t57 + t157 * t122;
t191 = t110 * t115;
t192 = t109 * t115;
t82 = Icges(6,5) * t110 - Icges(6,6) * t109;
t133 = t115 * t82;
t31 = t122 * t133 - t241;
t194 = qJD(1) * t57;
t32 = t120 * t133 + t194;
t34 = t59 * qJD(1) + t120 * t134;
t36 = t61 * qJD(1) + t120 * t135;
t2 = (t120 * t31 + (-t16 + t139) * qJD(1)) * t120 + (-t17 * qJD(1) + (t109 * t36 + t110 * t34 + t60 * t191 - t58 * t192 - t241) * t122 + (-t32 + t167 * t110 + t165 * t109 + (-t158 + t57) * qJD(1)) * t120) * t122;
t227 = t120 * t2;
t123 = -pkin(7) - pkin(6);
t225 = -t62 - t120 * t228 - (-pkin(6) - t123) * t122;
t181 = qJD(1) * t123;
t224 = -pkin(4) * t171 - t122 * t181;
t216 = t119 * t72;
t73 = -Icges(5,6) * t120 + t147 * t122;
t215 = t119 * t73;
t212 = t121 * t74;
t75 = -Icges(5,5) * t120 + t149 * t122;
t211 = t121 * t75;
t69 = t159 * t115;
t207 = t122 * t69;
t28 = -t85 * t183 - t207 + (-t121 * t183 - t172) * pkin(4);
t206 = t28 * t122;
t205 = -rSges(6,3) + t123;
t71 = -Icges(5,3) * t120 + t145 * t122;
t193 = qJD(1) * t71;
t187 = pkin(4) * t189 + t120 * t123;
t180 = qJD(4) * t119;
t179 = qJD(4) * t121;
t177 = -rSges(4,3) + t226;
t176 = pkin(4) * t179;
t175 = t85 * t182;
t170 = pkin(4) * t121 + t85;
t166 = t115 * t60 + t34;
t164 = -t115 * t58 + t36;
t163 = t184 * t69;
t90 = t160 * qJD(4);
t162 = t184 * t90;
t14 = t158 * t120 + t122 * t56;
t15 = t122 * t57 + t139;
t1 = t122 * ((t122 * t32 + (-t15 + t140) * qJD(1)) * t122 + (-t14 * qJD(1) + (-t109 * t35 - t110 * t33 - t61 * t191 + t59 * t192 + t194) * t120 + (-t31 + t166 * t110 + t164 * t109 + (-t157 - t56) * qJD(1)) * t122) * t120);
t7 = -t120 * t17 + t122 * t16;
t161 = -t7 * t182 + t1;
t153 = t119 * t74 + t121 * t72;
t152 = t119 * t75 + t121 * t73;
t128 = -t159 + t226 - t228;
t127 = t128 * t120;
t39 = t205 * t122 + t113 + t127;
t40 = t63 + t173 + t187;
t151 = t120 * t40 + t122 * t39;
t143 = t122 * t242;
t142 = rSges(3,3) * t122 + t230 * t120;
t67 = t146 * t115;
t68 = t148 * t115;
t124 = -qJD(1) * t82 + (-t67 + t135) * t110 + (-t68 - t134) * t109;
t141 = (t167 * t109 - t165 * t110 - t120 * t236 + t124 * t122) * t232 + (-t166 * t109 + t164 * t110 + t124 * t120 + t122 * t236) * t231 - (-t109 * t58 + t110 * t60 + t156 * t120 + t122 * t82) * t183 / 0.2e1 - (-t109 * t59 + t110 * t61 - t120 * t82 + t156 * t122) * t182 / 0.2e1;
t138 = t153 * t122;
t137 = t152 * t120;
t132 = qJD(4) * (Icges(5,1) * t121 - t204);
t131 = qJD(4) * (-Icges(5,2) * t119 + t203);
t130 = qJD(4) * (Icges(5,5) * t121 - Icges(5,6) * t119);
t129 = rSges(4,2) * t122 + t177 * t120;
t11 = (t127 - t210) * qJD(1) + t143 + t174 - t224;
t12 = t108 + (-qJD(3) - t176 - t242) * t120 + ((-qJ(2) - t205) * t120 + t128 * t122) * qJD(1);
t126 = t11 * t120 + t12 * t122 + (-t120 * t39 + t122 * t40) * qJD(1);
t80 = pkin(6) * t120 + t187;
t79 = -rSges(3,2) * t122 + rSges(3,3) * t120 + t185;
t78 = t113 + t142;
t77 = -rSges(5,3) * t120 + t188;
t76 = t122 * rSges(5,3) + t120 * t160;
t65 = rSges(4,2) * t120 + rSges(4,3) * t122 + t173;
t64 = t113 + t129;
t55 = t170 * t122;
t54 = t170 * t120;
t53 = t63 * t183;
t52 = t108 + (t230 * t122 + (-rSges(3,3) - qJ(2)) * t120) * qJD(1);
t51 = t142 * qJD(1) + t186;
t48 = -qJD(3) * t120 + t108 + ((-rSges(4,2) - qJ(2)) * t120 + t177 * t122) * qJD(1);
t47 = t129 * qJD(1) + t174;
t42 = t120 * t130 + t193;
t41 = t122 * t130 - t239;
t38 = t63 * qJD(1) + t120 * t242;
t37 = t143 - t240;
t30 = -t120 * t62 - t122 * t63;
t29 = t175 - t120 * t69 + (-t120 * t180 + t121 * t182) * pkin(4);
t21 = -t120 * t71 + t152 * t122;
t20 = -t120 * t70 + t138;
t19 = t122 * t71 + t137;
t18 = t153 * t120 + t122 * t70;
t13 = (-t63 - t80) * t122 + t225 * t120;
t10 = -t120 * t38 + t53 + (-t37 - t240) * t122;
t6 = -t120 * t15 + t122 * t14;
t3 = t53 + (-t37 + (-t122 * pkin(6) + t225) * qJD(1) + t224) * t122 + (qJD(1) * t80 - t105 - t38 + (-t176 - t181) * t120) * t120;
t4 = [(t11 * t40 + t12 * t39) * t234 + 0.2e1 * m(3) * (t51 * t79 + t52 * t78) + 0.2e1 * m(4) * (t47 * t65 + t48 * t64) - t119 * t132 - t149 * t179 - t121 * t131 + t147 * t180 + (t22 * t50 + t23 * t49) * t235 - t84 * t192 - t110 * t68 - t83 * t191 + t109 * t67; m(6) * (t151 * qJD(1) - t11 * t122 + t12 * t120) + m(3) * (t120 * t52 - t122 * t51 + (t120 * t79 + t122 * t78) * qJD(1)) + m(4) * (t120 * t48 - t122 * t47 + (t120 * t65 + t122 * t64) * qJD(1)) + m(5) * (qJD(1) * t238 + t120 * t23 - t122 * t22); 0; m(6) * t126 + m(4) * (t120 * t47 + t122 * t48 + (-t120 * t64 + t122 * t65) * qJD(1)) + m(5) * ((-t120 * t49 + t122 * t50) * qJD(1) + t237); 0; 0; (-t152 * qJD(4) - t119 * (-qJD(1) * t72 + t122 * t131) + t121 * (-qJD(1) * t74 + t122 * t132)) * t232 + (-t153 * qJD(4) - t119 * (t73 * qJD(1) + t120 * t131) + t121 * (t75 * qJD(1) + t120 * t132)) * t231 + m(6) * (t11 * t54 + t12 * t55 + t28 * t39 + t29 * t40) + m(5) * (t237 * t94 - t238 * t90) - (t116 / 0.2e1 + t117 / 0.2e1) * t145 * qJD(4) + ((t215 / 0.2e1 - t211 / 0.2e1 + t50 * t233) * t122 + (t216 / 0.2e1 - t212 / 0.2e1 - t49 * t233) * t120) * qJD(1) + t141; m(6) * (t120 * t28 - t122 * t29 + (t120 * t54 + t122 * t55) * qJD(1)); m(6) * (t29 * t120 + t206 + (-t120 * t55 + t122 * t54) * qJD(1)) - m(5) * t162; (-t94 * t162 + (-t122 * t97 + (-t94 * t116 + t117 * t221) * qJD(4) + (rSges(5,3) * t184 + t120 * t77 - t122 * t76) * qJD(1)) * (-t120 * t76 - t122 * t77)) * t235 - (-t120 * t21 + t122 * t20) * t182 - t120 * ((t120 * t41 + (-t20 + t137) * qJD(1)) * t120 + (-t21 * qJD(1) + (t74 * t179 - t72 * t180 - t239) * t122 + (-t42 + (-t211 + t215) * qJD(4) + (-t153 + t71) * qJD(1)) * t120) * t122) + t122 * ((t122 * t42 + (-t19 + t138) * qJD(1)) * t122 + (-t18 * qJD(1) + (-t75 * t179 + t73 * t180 + t193) * t120 + (-t41 + (t212 - t216) * qJD(4) + (-t152 - t70) * qJD(1)) * t122) * t120) + (t13 * t3 + t28 * t55 + t29 * t54) * t234 - t227 + t161 + (t120 * t19 - t122 * t18 - t6) * t183; m(6) * (t126 * t85 - t151 * t69) + t141; 0; -m(6) * t163; m(6) * (t10 * t13 + t54 * t175 + t85 * t206 - t55 * t207 + t3 * t30) + (m(6) * (-t54 * t69 + (-t55 * qJD(1) + t29) * t85) - t2 - qJD(1) * t6) * t120 + t161; (t10 * t30 - t85 * t163) * t234 - t227 + t1 + (-t120 * t6 - t122 * t7) * qJD(1);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
