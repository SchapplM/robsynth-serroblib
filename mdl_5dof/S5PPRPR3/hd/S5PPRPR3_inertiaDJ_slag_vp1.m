% Calculate time derivative of joint inertia matrix for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:48
% EndTime: 2019-12-05 15:04:57
% DurationCPUTime: 3.42s
% Computational Cost: add. (8836->426), mult. (14452->676), div. (0->0), fcn. (15146->10), ass. (0->191)
t160 = qJ(3) + pkin(9);
t158 = sin(t160);
t159 = cos(t160);
t164 = cos(pkin(7));
t162 = sin(pkin(7));
t163 = cos(pkin(8));
t186 = t162 * t163;
t132 = t158 * t186 + t159 * t164;
t122 = t132 * qJD(3);
t133 = -t158 * t164 + t159 * t186;
t123 = t133 * qJD(3);
t73 = -Icges(5,5) * t122 - Icges(5,6) * t123;
t169 = cos(qJ(3));
t181 = t164 * t169;
t167 = sin(qJ(3));
t185 = t162 * t167;
t150 = -t163 * t185 - t181;
t140 = t150 * qJD(3);
t182 = t164 * t167;
t184 = t162 * t169;
t151 = t163 * t184 - t182;
t141 = t151 * qJD(3);
t99 = Icges(4,5) * t140 - Icges(4,6) * t141;
t210 = t99 + t73;
t152 = -t163 * t182 + t184;
t142 = t152 * qJD(3);
t153 = t163 * t181 + t185;
t143 = t153 * qJD(3);
t100 = Icges(4,5) * t142 - Icges(4,6) * t143;
t183 = t163 * t164;
t134 = t158 * t183 - t162 * t159;
t124 = t134 * qJD(3);
t135 = t158 * t162 + t159 * t183;
t125 = t135 * qJD(3);
t74 = -Icges(5,5) * t124 - Icges(5,6) * t125;
t209 = t100 + t74;
t161 = sin(pkin(8));
t178 = qJD(3) * t161;
t126 = (-Icges(5,5) * t158 - Icges(5,6) * t159) * t178;
t144 = (-Icges(4,5) * t167 - Icges(4,6) * t169) * t178;
t208 = -t126 - t144;
t207 = 2 * m(6);
t206 = t163 ^ 2;
t205 = pkin(3) * t169;
t168 = cos(qJ(5));
t166 = sin(qJ(5));
t188 = t161 * t166;
t109 = t133 * t168 + t162 * t188;
t63 = -qJD(5) * t109 + t122 * t166;
t187 = t161 * t168;
t108 = -t133 * t166 + t162 * t187;
t64 = qJD(5) * t108 - t122 * t168;
t44 = rSges(6,1) * t64 + rSges(6,2) * t63 + rSges(6,3) * t123;
t203 = -pkin(4) * t122 + pkin(6) * t123 + t44;
t55 = rSges(6,1) * t109 + rSges(6,2) * t108 + rSges(6,3) * t132;
t202 = pkin(4) * t133 + pkin(6) * t132 + t55;
t189 = t161 * t164;
t170 = qJ(4) * t161 + t163 * t205;
t96 = pkin(3) * t185 + t164 * t170;
t201 = -rSges(5,1) * t135 + rSges(5,2) * t134 - rSges(5,3) * t189 - t96;
t118 = -qJ(4) * t163 + t161 * t205;
t190 = t161 * t162;
t95 = -pkin(3) * t182 + t162 * t170;
t200 = t118 * t190 + t163 * t95;
t177 = qJD(4) * t161;
t117 = pkin(3) * t142 + t164 * t177;
t199 = rSges(5,1) * t124 + rSges(5,2) * t125 - t117;
t149 = t159 * t187 - t163 * t166;
t173 = t158 * t178;
t112 = -qJD(5) * t149 + t166 * t173;
t148 = -t159 * t188 - t163 * t168;
t113 = qJD(5) * t148 - t168 * t173;
t172 = t159 * t178;
t60 = rSges(6,1) * t113 + rSges(6,2) * t112 + rSges(6,3) * t172;
t198 = (-pkin(4) * t158 + pkin(6) * t159) * t178 + t60;
t192 = t158 * t161;
t87 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t192;
t197 = (pkin(4) * t159 + pkin(6) * t158) * t161 + t87;
t196 = Icges(4,4) * t167;
t195 = Icges(4,4) * t169;
t194 = Icges(5,4) * t158;
t193 = Icges(5,4) * t159;
t116 = pkin(3) * t140 + t162 * t177;
t154 = -pkin(3) * t167 * t178 - qJD(4) * t163;
t180 = t163 * t116 + t154 * t190;
t179 = qJD(3) * t159;
t110 = -t135 * t166 + t164 * t187;
t111 = t135 * t168 + t164 * t188;
t56 = rSges(6,1) * t111 + rSges(6,2) * t110 + rSges(6,3) * t134;
t176 = -pkin(4) * t135 - pkin(6) * t134 - t56 - t96;
t65 = -qJD(5) * t111 + t124 * t166;
t66 = qJD(5) * t110 - t124 * t168;
t45 = rSges(6,1) * t66 + rSges(6,2) * t65 + rSges(6,3) * t125;
t175 = pkin(4) * t124 - pkin(6) * t125 - t117 - t45;
t105 = rSges(4,1) * t140 - rSges(4,2) * t141;
t147 = (-rSges(4,1) * t167 - rSges(4,2) * t169) * t178;
t61 = t105 * t163 + t147 * t190;
t106 = rSges(4,1) * t142 - rSges(4,2) * t143;
t62 = -t106 * t163 - t147 * t189;
t171 = t162 * t61 - t164 * t62;
t146 = (-Icges(4,1) * t167 - t195) * t178;
t145 = (-Icges(4,2) * t169 - t196) * t178;
t137 = -Icges(4,5) * t163 + (Icges(4,1) * t169 - t196) * t161;
t136 = -Icges(4,6) * t163 + (-Icges(4,2) * t167 + t195) * t161;
t129 = (-rSges(5,1) * t158 - rSges(5,2) * t159) * t178;
t128 = (-Icges(5,1) * t158 - t193) * t178;
t127 = (-Icges(5,2) * t159 - t194) * t178;
t121 = -rSges(5,3) * t163 + (rSges(5,1) * t159 - rSges(5,2) * t158) * t161;
t120 = -Icges(5,5) * t163 + (Icges(5,1) * t159 - t194) * t161;
t119 = -Icges(5,6) * t163 + (-Icges(5,2) * t158 + t193) * t161;
t107 = t116 * t189;
t104 = Icges(4,1) * t142 - Icges(4,4) * t143;
t103 = Icges(4,1) * t140 - Icges(4,4) * t141;
t102 = Icges(4,4) * t142 - Icges(4,2) * t143;
t101 = Icges(4,4) * t140 - Icges(4,2) * t141;
t98 = rSges(4,1) * t153 + rSges(4,2) * t152 + rSges(4,3) * t189;
t97 = rSges(4,1) * t151 + rSges(4,2) * t150 + rSges(4,3) * t190;
t94 = Icges(4,1) * t153 + Icges(4,4) * t152 + Icges(4,5) * t189;
t93 = Icges(4,1) * t151 + Icges(4,4) * t150 + Icges(4,5) * t190;
t92 = Icges(4,4) * t153 + Icges(4,2) * t152 + Icges(4,6) * t189;
t91 = Icges(4,4) * t151 + Icges(4,2) * t150 + Icges(4,6) * t190;
t86 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t192;
t85 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t192;
t84 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t192;
t81 = t95 * t189;
t79 = -rSges(5,1) * t122 - rSges(5,2) * t123;
t78 = -Icges(5,1) * t124 - Icges(5,4) * t125;
t77 = -Icges(5,1) * t122 - Icges(5,4) * t123;
t76 = -Icges(5,4) * t124 - Icges(5,2) * t125;
t75 = -Icges(5,4) * t122 - Icges(5,2) * t123;
t71 = rSges(5,1) * t133 - rSges(5,2) * t132 + rSges(5,3) * t190;
t70 = Icges(5,1) * t135 - Icges(5,4) * t134 + Icges(5,5) * t189;
t69 = Icges(5,1) * t133 - Icges(5,4) * t132 + Icges(5,5) * t190;
t68 = Icges(5,4) * t135 - Icges(5,2) * t134 + Icges(5,6) * t189;
t67 = Icges(5,4) * t133 - Icges(5,2) * t132 + Icges(5,6) * t190;
t59 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t172;
t58 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t172;
t57 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t172;
t54 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t134;
t53 = Icges(6,1) * t109 + Icges(6,4) * t108 + Icges(6,5) * t132;
t52 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t134;
t51 = Icges(6,4) * t109 + Icges(6,2) * t108 + Icges(6,6) * t132;
t50 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t134;
t49 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t132;
t48 = (t105 * t164 - t106 * t162) * t161;
t47 = t199 * t163 + (-t129 - t154) * t189;
t46 = t129 * t190 + t163 * t79 + t180;
t43 = Icges(6,1) * t66 + Icges(6,4) * t65 + Icges(6,5) * t125;
t42 = Icges(6,1) * t64 + Icges(6,4) * t63 + Icges(6,5) * t123;
t41 = Icges(6,4) * t66 + Icges(6,2) * t65 + Icges(6,6) * t125;
t40 = Icges(6,4) * t64 + Icges(6,2) * t63 + Icges(6,6) * t123;
t39 = Icges(6,5) * t66 + Icges(6,6) * t65 + Icges(6,3) * t125;
t38 = Icges(6,5) * t64 + Icges(6,6) * t63 + Icges(6,3) * t123;
t37 = -t134 * t87 + t192 * t56;
t36 = t132 * t87 - t192 * t55;
t35 = t107 + (t162 * t199 + t164 * t79) * t161;
t34 = t148 * t85 + t149 * t86 + t192 * t84;
t33 = -t132 * t56 + t134 * t55;
t32 = t110 * t85 + t111 * t86 + t134 * t84;
t31 = t108 * t85 + t109 * t86 + t132 * t84;
t30 = t176 * t163 + (-t118 - t197) * t189;
t29 = t163 * t202 + t190 * t197 + t200;
t28 = t148 * t52 + t149 * t54 + t192 * t50;
t27 = t148 * t51 + t149 * t53 + t192 * t49;
t26 = t110 * t52 + t111 * t54 + t134 * t50;
t25 = t110 * t51 + t111 * t53 + t134 * t49;
t24 = t108 * t52 + t109 * t54 + t132 * t50;
t23 = t108 * t51 + t109 * t53 + t132 * t49;
t22 = t175 * t163 + (-t154 - t198) * t189;
t21 = t163 * t203 + t190 * t198 + t180;
t20 = t81 + (t162 * t176 + t164 * t202) * t161;
t19 = -t125 * t87 - t134 * t60 + (t158 * t45 + t179 * t56) * t161;
t18 = t123 * t87 + t132 * t60 + (-t158 * t44 - t55 * t179) * t161;
t17 = t107 + (t162 * t175 + t164 * t203) * t161;
t16 = t112 * t85 + t113 * t86 + t148 * t58 + t149 * t59 + (t158 * t57 + t179 * t84) * t161;
t15 = -t123 * t56 + t125 * t55 - t132 * t45 + t134 * t44;
t14 = t110 * t58 + t111 * t59 + t125 * t84 + t134 * t57 + t65 * t85 + t66 * t86;
t13 = t108 * t58 + t109 * t59 + t123 * t84 + t132 * t57 + t63 * t85 + t64 * t86;
t12 = t112 * t52 + t113 * t54 + t148 * t41 + t149 * t43 + (t158 * t39 + t179 * t50) * t161;
t11 = t112 * t51 + t113 * t53 + t148 * t40 + t149 * t42 + (t158 * t38 + t179 * t49) * t161;
t10 = t110 * t41 + t111 * t43 + t125 * t50 + t134 * t39 + t52 * t65 + t54 * t66;
t9 = t110 * t40 + t111 * t42 + t125 * t49 + t134 * t38 + t51 * t65 + t53 * t66;
t8 = t108 * t41 + t109 * t43 + t123 * t50 + t132 * t39 + t52 * t63 + t54 * t64;
t7 = t108 * t40 + t109 * t42 + t123 * t49 + t132 * t38 + t51 * t63 + t53 * t64;
t6 = -t16 * t163 + (t11 * t162 + t12 * t164) * t161;
t5 = -t14 * t163 + (t10 * t164 + t162 * t9) * t161;
t4 = -t13 * t163 + (t162 * t7 + t164 * t8) * t161;
t3 = t11 * t132 + t12 * t134 + t123 * t27 + t125 * t28 + (t158 * t16 + t179 * t34) * t161;
t2 = t10 * t134 + t123 * t25 + t26 * t125 + t132 * t9 + (t14 * t158 + t179 * t32) * t161;
t1 = t123 * t23 + t125 * t24 + t132 * t7 + t134 * t8 + (t13 * t158 + t179 * t31) * t161;
t72 = [0; 0; 0; m(4) * t48 + m(5) * t35 + m(6) * t17; m(4) * t171 + m(5) * (t162 * t46 - t164 * t47) + m(6) * (t162 * t21 - t164 * t22); -t163 * t6 + (t17 * t20 + t21 * t29 + t22 * t30) * t207 + 0.2e1 * m(5) * (t200 * t46 + t81 * t35 + (t201 * t47 + t71 * t46) * t163 + (((-t118 - t121) * t47 + t71 * t35) * t164 + (t121 * t46 + t201 * t35) * t162) * t161) + 0.2e1 * m(4) * ((t61 * t97 - t62 * t98) * t163 + ((-t162 * t98 + t164 * t97) * t48 + t171 * (-t163 * rSges(4,3) + (rSges(4,1) * t169 - rSges(4,2) * t167) * t161)) * t161) - t163 * (t206 * t126 + (((-t158 * t76 + t159 * t78) * t164 + (-t158 * t75 + t159 * t77) * t162 + ((-t158 * t70 - t159 * t68) * t164 + (-t158 * t69 - t159 * t67) * t162) * qJD(3)) * t161 + (t127 * t158 - t128 * t159 - t73 * t162 - t74 * t164 + (t119 * t159 + t120 * t158) * qJD(3)) * t163) * t161) - t163 * (t206 * t144 + (((-t102 * t167 + t104 * t169) * t164 + (-t101 * t167 + t103 * t169) * t162 + ((-t167 * t94 - t169 * t92) * t164 + (-t167 * t93 - t169 * t91) * t162) * qJD(3)) * t161 + (-t100 * t164 + t145 * t167 - t146 * t169 - t99 * t162 + (t136 * t169 + t137 * t167) * qJD(3)) * t163) * t161) + (t4 + (t101 * t150 + t103 * t151 - t122 * t69 - t123 * t67 - t132 * t75 + t133 * t77 + t140 * t93 - t141 * t91 + t210 * t190) * t190 + (t119 * t123 + t120 * t122 + t127 * t132 - t128 * t133 + t136 * t141 - t140 * t137 - t145 * t150 - t151 * t146 + t208 * t190) * t163) * t190 + (t5 + (t102 * t152 + t104 * t153 - t124 * t70 - t125 * t68 - t134 * t76 + t135 * t78 + t142 * t94 - t143 * t92 + t209 * t189) * t189 + (t119 * t125 + t120 * t124 + t127 * t134 - t128 * t135 + t136 * t143 - t137 * t142 - t145 * t152 - t146 * t153 + t208 * t189) * t163 + (t101 * t152 + t103 * t153 + t142 * t93 - t143 * t91 + t102 * t150 + t104 * t151 + t140 * t94 - t141 * t92 - t124 * t69 - t125 * t67 - t134 * t75 + t135 * t77 - t122 * t70 - t123 * t68 - t132 * t76 + t133 * t78 + t210 * t189 + t209 * t190) * t190) * t189; 0; 0; m(6) * (-t163 * t17 + (t162 * t22 + t164 * t21) * t161) + m(5) * (-t163 * t35 + (t162 * t47 + t164 * t46) * t161); 0; m(6) * t15; m(6) * (t162 * t18 - t164 * t19); t132 * t4 / 0.2e1 + t134 * t5 / 0.2e1 + m(6) * (t15 * t20 + t17 * t33 + t18 * t29 + t19 * t30 + t21 * t36 + t22 * t37) + (-t3 / 0.2e1 - t123 * t31 / 0.2e1 - t125 * t32 / 0.2e1) * t163 + (t164 * t2 / 0.2e1 + t162 * t1 / 0.2e1 + t123 * (t162 * t23 + t164 * t24) / 0.2e1 + t125 * (t162 * t25 + t164 * t26) / 0.2e1 + t158 * t6 / 0.2e1 + (-t34 * t163 / 0.2e1 + (t162 * t27 + t164 * t28) * t161 / 0.2e1) * t179) * t161; m(6) * (-t15 * t163 + (t162 * t19 + t164 * t18) * t161); (t15 * t33 + t18 * t36 + t19 * t37) * t207 + t125 * (t132 * t25 + t134 * t26 + t192 * t32) + t134 * t2 + t123 * (t132 * t23 + t134 * t24 + t192 * t31) + t132 * t1 + (t132 * t27 + t134 * t28 + t192 * t34) * t172 + t3 * t192;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t72(1), t72(2), t72(4), t72(7), t72(11); t72(2), t72(3), t72(5), t72(8), t72(12); t72(4), t72(5), t72(6), t72(9), t72(13); t72(7), t72(8), t72(9), t72(10), t72(14); t72(11), t72(12), t72(13), t72(14), t72(15);];
Mq = res;
