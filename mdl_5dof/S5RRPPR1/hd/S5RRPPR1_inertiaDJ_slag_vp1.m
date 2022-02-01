% Calculate time derivative of joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:19
% EndTime: 2022-01-20 09:51:25
% DurationCPUTime: 1.58s
% Computational Cost: add. (3820->195), mult. (2408->274), div. (0->0), fcn. (1712->10), ass. (0->123)
t109 = pkin(9) + qJ(5);
t103 = sin(t109);
t146 = qJD(5) * t103;
t104 = cos(t109);
t145 = qJD(5) * t104;
t111 = qJ(1) + qJ(2);
t105 = pkin(8) + t111;
t100 = cos(t105);
t110 = qJD(1) + qJD(2);
t147 = t100 * t110;
t107 = cos(t111);
t102 = pkin(2) * t107;
t177 = rSges(5,3) + qJ(4);
t99 = sin(t105);
t178 = t177 * t99 + t102;
t113 = cos(pkin(9));
t137 = -rSges(5,1) * t113 - pkin(3);
t106 = sin(t111);
t76 = t107 * rSges(3,1) - rSges(3,2) * t106;
t114 = -pkin(7) - qJ(4);
t176 = -t114 * t99 + t102;
t157 = rSges(6,2) * t103;
t160 = rSges(6,1) * t104;
t175 = -t157 + t160;
t148 = Icges(6,4) * t104;
t127 = -Icges(6,2) * t103 + t148;
t41 = Icges(6,6) * t99 + t127 * t100;
t149 = Icges(6,4) * t103;
t128 = Icges(6,1) * t104 - t149;
t43 = Icges(6,5) * t99 + t128 * t100;
t130 = t103 * t41 - t104 * t43;
t174 = t130 * t99;
t173 = -t100 * rSges(4,1) - t102;
t40 = -Icges(6,6) * t100 + t127 * t99;
t42 = -Icges(6,5) * t100 + t128 * t99;
t131 = t103 * t40 - t104 * t42;
t172 = t131 * t100;
t126 = Icges(6,5) * t104 - Icges(6,6) * t103;
t68 = Icges(6,2) * t104 + t149;
t69 = Icges(6,1) * t103 + t148;
t129 = t103 * t68 - t104 * t69;
t171 = t126 * qJD(5) + t129 * t110;
t170 = 2 * m(3);
t169 = 2 * m(4);
t168 = 2 * m(5);
t167 = 2 * m(6);
t115 = sin(qJ(1));
t164 = pkin(1) * t115;
t163 = pkin(2) * t106;
t92 = t99 * rSges(6,3);
t162 = t100 * rSges(6,3) + t99 * t157;
t158 = rSges(5,2) * sin(pkin(9));
t156 = pkin(1) * qJD(1);
t153 = t110 * t41;
t152 = t110 * t43;
t151 = t110 * t99;
t142 = t110 * t157;
t144 = -t100 * t142 + (-rSges(6,1) * t146 - rSges(6,2) * t145) * t99;
t143 = t110 * t158;
t141 = t115 * t156;
t116 = cos(qJ(1));
t140 = t116 * t156;
t101 = pkin(4) * t113 + pkin(3);
t132 = -t101 - t160;
t118 = -t100 * t114 + t132 * t99 - t163;
t74 = rSges(6,1) * t103 + rSges(6,2) * t104;
t119 = -t74 * t100 * qJD(5) + rSges(6,3) * t147 + t99 * t142;
t88 = qJD(4) * t99;
t12 = t118 * t110 + t119 + t88;
t10 = t12 - t141;
t30 = t118 + t162;
t28 = t30 - t164;
t136 = t110 * t28 - t10;
t89 = qJD(4) * t100;
t13 = t89 - t144 + (t132 * t100 - t176 - t92) * t110;
t11 = t13 - t140;
t108 = t116 * pkin(1);
t45 = t100 * t175 + t92;
t31 = t100 * t101 + t176 + t45;
t29 = t108 + t31;
t135 = t110 * t29 + t11;
t134 = t110 * t30 - t12;
t133 = t110 * t31 + t13;
t55 = -rSges(4,2) * t99 - t173;
t58 = t76 * t110;
t75 = -rSges(3,1) * t106 - rSges(3,2) * t107;
t67 = Icges(6,5) * t103 + Icges(6,6) * t104;
t125 = (t128 - t68) * t146 + (t127 + t69) * t145;
t124 = (-t130 * qJD(5) + t103 * (Icges(6,5) * t147 - t128 * t151) + t104 * (Icges(6,6) * t147 - t127 * t151) + t171 * t99) * t99 / 0.2e1 - (-t131 * qJD(5) - t171 * t100 + t103 * t152 + t104 * t153) * t100 / 0.2e1 + (-t100 * t67 + t103 * t42 + t104 * t40 - t129 * t99) * t151 / 0.2e1 + (-t129 * t100 + t103 * t43 + t104 * t41 + t67 * t99) * t147 / 0.2e1;
t57 = t75 * t110;
t121 = t67 * qJD(5);
t54 = -rSges(4,1) * t99 - rSges(4,2) * t100 - t163;
t120 = t137 * t99 - t163;
t47 = rSges(4,2) * t151 + t173 * t110;
t35 = (-t158 - t137) * t100 + t178;
t46 = t54 * t110;
t39 = Icges(6,3) * t99 + t126 * t100;
t34 = t100 * t177 + t99 * t158 + t120;
t20 = t120 * t110 + t99 * t143 + t147 * t177 + t88;
t21 = t100 * t143 + t89 + (t137 * t100 - t178) * t110;
t64 = t175 * qJD(5);
t60 = t108 + t76;
t59 = t75 - t164;
t51 = -t58 - t140;
t50 = t57 - t141;
t49 = t108 + t55;
t48 = t54 - t164;
t44 = t99 * t160 - t162;
t38 = -Icges(6,3) * t100 + t126 * t99;
t37 = t47 - t140;
t36 = t46 - t141;
t33 = t108 + t35;
t32 = t34 - t164;
t23 = t110 * t39 - t99 * t121;
t22 = -t126 * t151 + (Icges(6,3) * t110 - t121) * t100;
t17 = t21 - t140;
t16 = t20 - t141;
t9 = -t130 * t100 + t39 * t99;
t8 = t38 * t99 - t172;
t7 = -t100 * t39 - t174;
t6 = -t100 * t38 - t131 * t99;
t3 = ((-t45 + t92) * t110 + t144) * t99 + (t110 * t44 + t119) * t100;
t1 = [(t10 * t29 + t11 * t28) * t167 + (t36 * t49 + t37 * t48) * t169 + (t16 * t33 + t17 * t32) * t168 + (t50 * t60 + t51 * t59) * t170 + t125; m(6) * (t10 * t31 + t11 * t30 + t12 * t29 + t13 * t28) + m(4) * (t36 * t55 + t37 * t54 + t46 * t49 + t47 * t48) + m(5) * (t16 * t35 + t17 * t34 + t20 * t33 + t21 * t32) + m(3) * (t50 * t76 + t51 * t75 + t57 * t60 - t58 * t59) + t125; (t12 * t31 + t13 * t30) * t167 + (t46 * t55 + t47 * t54) * t169 + (t20 * t35 + t21 * t34) * t168 + (t57 * t76 - t58 * t75) * t170 + t125; 0; 0; 0; m(6) * (t136 * t100 + t135 * t99) + m(5) * ((t110 * t33 + t17) * t99 + (t110 * t32 - t16) * t100); m(6) * (t134 * t100 + t133 * t99) + m(5) * ((t110 * t35 + t21) * t99 + (t110 * t34 - t20) * t100); 0; 0; m(6) * ((-t100 * t28 - t29 * t99) * t64 + (-t135 * t100 + t136 * t99) * t74) + t124; m(6) * ((-t100 * t30 - t31 * t99) * t64 + (-t133 * t100 + t134 * t99) * t74) + t124; m(6) * t3; 0; ((t100 * t45 + t44 * t99) * t3 + (t100 ^ 2 + t99 ^ 2) * t74 * t64) * t167 + (-t100 * t8 + t9 * t99) * t147 + t99 * ((t99 * t22 + (t8 + t174) * t110) * t99 + (t9 * t110 + (t40 * t145 + t42 * t146) * t100 + (-t23 + (-qJD(5) * t41 + t110 * t42) * t104 + (-qJD(5) * t43 - t110 * t40) * t103) * t99) * t100) + (-t100 * t6 + t7 * t99) * t151 - t100 * ((t100 * t23 + (t7 + t172) * t110) * t100 + (t6 * t110 + (-t41 * t145 - t43 * t146) * t99 + (-t22 + (qJD(5) * t40 + t152) * t104 + (qJD(5) * t42 - t153) * t103) * t100) * t99);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
