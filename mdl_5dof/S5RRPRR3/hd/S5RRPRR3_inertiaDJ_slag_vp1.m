% Calculate time derivative of joint inertia matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:59
% EndTime: 2022-01-20 10:34:04
% DurationCPUTime: 1.53s
% Computational Cost: add. (4960->203), mult. (2858->285), div. (0->0), fcn. (1986->10), ass. (0->125)
t112 = sin(qJ(5));
t144 = qJD(5) * t112;
t114 = cos(qJ(5));
t143 = qJD(5) * t114;
t110 = qJD(1) + qJD(2);
t105 = qJD(4) + t110;
t177 = t105 * t112;
t176 = t105 * t114;
t111 = qJ(1) + qJ(2);
t107 = sin(t111);
t108 = cos(t111);
t75 = t108 * rSges(3,1) - rSges(3,2) * t107;
t156 = rSges(6,2) * t112;
t159 = rSges(6,1) * t114;
t175 = -t156 + t159;
t106 = pkin(9) + t111;
t101 = sin(t106);
t102 = cos(t106);
t103 = pkin(2) * t108;
t174 = -t102 * rSges(4,1) + rSges(4,2) * t101 - t103;
t145 = Icges(6,4) * t114;
t125 = -Icges(6,2) * t112 + t145;
t104 = qJ(4) + t106;
t98 = sin(t104);
t99 = cos(t104);
t43 = Icges(6,6) * t98 + t125 * t99;
t146 = Icges(6,4) * t112;
t126 = Icges(6,1) * t114 - t146;
t45 = Icges(6,5) * t98 + t126 * t99;
t128 = t112 * t43 - t114 * t45;
t173 = t128 * t98;
t42 = -Icges(6,6) * t99 + t125 * t98;
t44 = -Icges(6,5) * t99 + t126 * t98;
t130 = t112 * t42 - t114 * t44;
t172 = t130 * t99;
t147 = pkin(3) * t102 + t103;
t124 = Icges(6,5) * t114 - Icges(6,6) * t112;
t86 = Icges(6,2) * t114 + t146;
t87 = Icges(6,1) * t112 + t145;
t127 = t112 * t86 - t114 * t87;
t170 = t124 * qJD(5) + t127 * t105;
t40 = -Icges(6,3) * t99 + t124 * t98;
t169 = 2 * m(3);
t168 = 2 * m(4);
t167 = 2 * m(5);
t166 = 2 * m(6);
t113 = sin(qJ(1));
t163 = pkin(1) * t113;
t162 = pkin(2) * t107;
t90 = t98 * rSges(6,3);
t141 = t105 * t156;
t149 = t105 * t99;
t161 = rSges(6,3) * t149 + t98 * t141;
t160 = t99 * rSges(6,3) + t98 * t156;
t152 = pkin(1) * qJD(1);
t41 = Icges(6,3) * t98 + t124 * t99;
t151 = t105 * t41;
t150 = t105 * t98;
t148 = t114 * t98;
t138 = rSges(6,2) * t143;
t142 = -t99 * t141 + (-t144 * rSges(6,1) - t138) * t98;
t140 = t113 * t152;
t115 = cos(qJ(1));
t139 = t115 * t152;
t135 = -pkin(4) - t159;
t63 = t99 * rSges(5,1) - rSges(5,2) * t98;
t53 = -rSges(5,1) * t149 + rSges(5,2) * t150;
t61 = t75 * t110;
t62 = -rSges(5,1) * t98 - rSges(5,2) * t99;
t134 = -pkin(3) * t101 - t162;
t74 = -rSges(3,1) * t107 - rSges(3,2) * t108;
t89 = rSges(6,1) * t112 + rSges(6,2) * t114;
t131 = t112 * t44 + t114 * t42;
t129 = t112 * t45 + t114 * t43;
t47 = t175 * t99 + t90;
t51 = t63 + t147;
t85 = Icges(6,5) * t112 + Icges(6,6) * t114;
t52 = t62 * t105;
t123 = t134 * t110;
t122 = t147 * t110;
t121 = (t126 - t86) * t144 + (t125 + t87) * t143;
t120 = (-t128 * qJD(5) + t170 * t98 - t42 * t176 - t44 * t177) * t98 / 0.2e1 - (-t130 * qJD(5) - t170 * t99 + t43 * t176 + t45 * t177) * t99 / 0.2e1 + (-t127 * t98 - t85 * t99 + t131) * t150 / 0.2e1 + (-t127 * t99 + t85 * t98 + t129) * t149 / 0.2e1;
t60 = t74 * t110;
t33 = t99 * pkin(4) + t98 * pkin(8) + t47;
t117 = qJD(5) * t85;
t58 = -rSges(4,1) * t101 - rSges(4,2) * t102 - t162;
t49 = t174 * t110;
t32 = t99 * pkin(8) + t135 * t98 + t160;
t29 = t33 + t147;
t48 = t58 * t110;
t50 = t134 + t62;
t35 = -t122 + t53;
t28 = t134 + t32;
t34 = t123 + t52;
t15 = (t135 * t99 + (-rSges(6,3) - pkin(8)) * t98) * t105 - t142;
t14 = -t99 * t138 - pkin(4) * t150 + pkin(8) * t149 + (-t105 * t148 - t99 * t144) * rSges(6,1) + t161;
t13 = -t122 + t15;
t12 = t123 + t14;
t109 = t115 * pkin(1);
t81 = t175 * qJD(5);
t65 = t109 + t75;
t64 = t74 - t163;
t57 = -t61 - t139;
t56 = t60 - t140;
t55 = t109 - t174;
t54 = t58 - t163;
t46 = rSges(6,1) * t148 - t160;
t39 = t109 + t51;
t38 = t50 - t163;
t37 = t49 - t139;
t36 = t48 - t140;
t31 = t35 - t139;
t30 = t34 - t140;
t27 = t109 + t29;
t26 = t28 - t163;
t19 = -t98 * t117 + t151;
t18 = -t105 * t40 - t99 * t117;
t11 = t13 - t139;
t10 = t12 - t140;
t9 = -t128 * t99 + t41 * t98;
t8 = t40 * t98 - t172;
t7 = -t41 * t99 - t173;
t6 = -t130 * t98 - t40 * t99;
t1 = ((-t47 + t90) * t105 + t142) * t98 + (-t89 * t99 * qJD(5) + t105 * t46 + t161) * t99;
t2 = [(t10 * t27 + t11 * t26) * t166 + (t30 * t39 + t31 * t38) * t167 + (t56 * t65 + t57 * t64) * t169 + (t36 * t55 + t37 * t54) * t168 + t121; m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + m(5) * (t30 * t51 + t31 * t50 + t34 * t39 + t35 * t38) + m(3) * (t56 * t75 + t57 * t74 + t60 * t65 - t61 * t64) + m(4) * (-t174 * t36 + t37 * t58 + t48 * t55 + t49 * t54) + t121; (t12 * t29 + t13 * t28) * t166 + (t34 * t51 + t35 * t50) * t167 + (-t174 * t48 + t49 * t58) * t168 + (t60 * t75 - t61 * t74) * t169 + t121; 0; 0; 0; m(6) * (t10 * t33 + t11 * t32 + t14 * t27 + t15 * t26) + m(5) * (t30 * t63 + t31 * t62 + t38 * t53 + t39 * t52) + t121; m(6) * (t12 * t33 + t13 * t32 + t14 * t29 + t15 * t28) + m(5) * (t34 * t63 + t35 * t62 + t50 * t53 + t51 * t52) + t121; 0; (t52 * t63 + t53 * t62) * t167 + (t14 * t33 + t15 * t32) * t166 + t121; m(6) * ((-t26 * t99 - t27 * t98) * t81 + (-t10 * t98 - t11 * t99 + (t26 * t98 - t27 * t99) * t105) * t89) + t120; m(6) * ((-t28 * t99 - t29 * t98) * t81 + (-t12 * t98 - t13 * t99 + (t28 * t98 - t29 * t99) * t105) * t89) + t120; m(6) * t1; m(6) * ((-t32 * t99 - t33 * t98) * t81 + (-t14 * t98 - t15 * t99 + (t32 * t98 - t33 * t99) * t105) * t89) + t120; ((t46 * t98 + t47 * t99) * t1 + (t98 ^ 2 + t99 ^ 2) * t89 * t81) * t166 + (-t8 * t99 + t9 * t98) * t149 + t98 * ((t98 * t18 + (t8 + t173) * t105) * t98 + (t9 * t105 + (t42 * t143 + t44 * t144) * t99 + (-t19 - t129 * qJD(5) + (-t130 + t41) * t105) * t98) * t99) + (-t6 * t99 + t7 * t98) * t150 - t99 * ((t99 * t19 + (t7 + t172) * t105) * t99 + (t6 * t105 + (-t43 * t143 - t45 * t144 + t151) * t98 + (t131 * qJD(5) - t128 * t105 - t18) * t99) * t98);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
