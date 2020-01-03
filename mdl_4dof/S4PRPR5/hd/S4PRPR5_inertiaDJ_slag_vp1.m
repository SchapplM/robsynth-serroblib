% Calculate time derivative of joint inertia matrix for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:23:03
% DurationCPUTime: 2.30s
% Computational Cost: add. (3704->227), mult. (5864->387), div. (0->0), fcn. (5543->8), ass. (0->132)
t101 = qJ(2) + pkin(7);
t100 = cos(t101);
t106 = sin(qJ(2));
t108 = cos(qJ(2));
t99 = sin(t101);
t197 = (-Icges(3,5) * t106 - Icges(4,5) * t99 - Icges(3,6) * t108 - Icges(4,6) * t100) * qJD(2);
t103 = cos(pkin(6));
t181 = t103 ^ 2;
t102 = sin(pkin(6));
t182 = t102 ^ 2;
t196 = t181 + t182;
t195 = t197 * t102;
t194 = t197 * t103;
t186 = t196 * qJD(2);
t183 = 2 * m(5);
t180 = t102 / 0.2e1;
t179 = t103 / 0.2e1;
t177 = pkin(2) * t106;
t174 = t196 * pkin(2) * t108;
t171 = pkin(2) * qJD(2);
t173 = t196 * t106 * t171;
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t121 = Icges(5,5) * t107 - Icges(5,6) * t105;
t157 = qJD(4) * t99;
t170 = t100 * ((-Icges(5,5) * t105 - Icges(5,6) * t107) * t157 + (Icges(5,3) * t99 + t121 * t100) * qJD(2));
t62 = -Icges(5,3) * t100 + t121 * t99;
t169 = t100 * t62;
t136 = rSges(5,1) * t107 - rSges(5,2) * t105;
t65 = -rSges(5,3) * t100 + t136 * t99;
t168 = t100 * t65;
t164 = t103 * t99;
t166 = t102 * t99;
t153 = t103 * t107;
t156 = t102 * t105;
t86 = -t100 * t156 - t153;
t154 = t103 * t105;
t155 = t102 * t107;
t87 = t100 * t155 - t154;
t45 = Icges(5,5) * t87 + Icges(5,6) * t86 + Icges(5,3) * t166;
t47 = Icges(5,4) * t87 + Icges(5,2) * t86 + Icges(5,6) * t166;
t49 = Icges(5,1) * t87 + Icges(5,4) * t86 + Icges(5,5) * t166;
t88 = -t100 * t154 + t155;
t89 = t100 * t153 + t156;
t18 = t45 * t164 + t47 * t88 + t49 * t89;
t167 = t102 * t18;
t46 = Icges(5,5) * t89 + Icges(5,6) * t88 + Icges(5,3) * t164;
t48 = Icges(5,4) * t89 + Icges(5,2) * t88 + Icges(5,6) * t164;
t50 = Icges(5,1) * t89 + Icges(5,4) * t88 + Icges(5,5) * t164;
t17 = t46 * t166 + t48 * t86 + t50 * t87;
t165 = t103 * t17;
t160 = Icges(5,4) * t105;
t159 = Icges(5,4) * t107;
t158 = qJD(2) * t99;
t152 = qJD(2) * t100;
t150 = t108 * t171;
t149 = t105 * t158;
t148 = t107 * t158;
t147 = t102 * t152;
t146 = t103 * t152;
t94 = rSges(4,1) * t99 + rSges(4,2) * t100;
t145 = -t94 - t177;
t95 = pkin(3) * t99 - pkin(5) * t100;
t144 = -t65 - t95 - t177;
t141 = rSges(4,1) * t100 - rSges(4,2) * t99;
t143 = -t141 * qJD(2) - t150;
t142 = pkin(3) * t100 + pkin(5) * t99;
t97 = t106 * rSges(3,1) + rSges(3,2) * t108;
t133 = -t105 * t47 + t107 * t49;
t20 = -t100 * t45 + t133 * t99;
t132 = -t105 * t48 + t107 * t50;
t21 = -t100 * t46 + t132 * t99;
t135 = t20 * t102 + t21 * t103;
t51 = rSges(5,1) * t87 + rSges(5,2) * t86 + rSges(5,3) * t166;
t52 = rSges(5,1) * t89 + rSges(5,2) * t88 + rSges(5,3) * t164;
t134 = -t102 * t52 + t103 * t51;
t122 = -Icges(5,2) * t105 + t159;
t63 = -Icges(5,6) * t100 + t122 * t99;
t124 = Icges(5,1) * t107 - t160;
t64 = -Icges(5,5) * t100 + t124 * t99;
t131 = t105 * t63 - t107 * t64;
t44 = (-rSges(5,1) * t105 - rSges(5,2) * t107) * t157 + (rSges(5,3) * t99 + t136 * t100) * qJD(2);
t120 = -t142 * qJD(2) - t150 - t44;
t54 = -t87 * qJD(4) + t102 * t149;
t55 = t86 * qJD(4) - t102 * t148;
t31 = Icges(5,5) * t55 + Icges(5,6) * t54 + Icges(5,3) * t147;
t119 = t45 * t152 + t31 * t99;
t56 = -t89 * qJD(4) + t103 * t149;
t57 = t88 * qJD(4) - t103 * t148;
t32 = Icges(5,5) * t57 + Icges(5,6) * t56 + Icges(5,3) * t146;
t118 = t46 * t152 + t32 * t99;
t28 = -t94 * t186 - t173;
t117 = t141 * t28;
t59 = t143 * t103;
t58 = t143 * t102;
t53 = t97 * t186;
t43 = (-Icges(5,1) * t105 - t159) * t157 + (Icges(5,5) * t99 + t124 * t100) * qJD(2);
t42 = (-Icges(5,2) * t107 - t160) * t157 + (Icges(5,6) * t99 + t122 * t100) * qJD(2);
t40 = t144 * t103;
t39 = t144 * t102;
t38 = rSges(5,1) * t57 + rSges(5,2) * t56 + rSges(5,3) * t146;
t37 = rSges(5,1) * t55 + rSges(5,2) * t54 + rSges(5,3) * t147;
t36 = Icges(5,1) * t57 + Icges(5,4) * t56 + Icges(5,5) * t146;
t35 = Icges(5,1) * t55 + Icges(5,4) * t54 + Icges(5,5) * t147;
t34 = Icges(5,4) * t57 + Icges(5,2) * t56 + Icges(5,6) * t146;
t33 = Icges(5,4) * t55 + Icges(5,2) * t54 + Icges(5,6) * t147;
t30 = t120 * t103;
t29 = t120 * t102;
t27 = -t100 * t52 - t65 * t164;
t26 = t100 * t51 + t65 * t166;
t25 = -t131 * t99 - t169;
t24 = t134 * t99;
t23 = t62 * t164 + t63 * t88 + t64 * t89;
t22 = t62 * t166 + t63 * t86 + t64 * t87;
t19 = t46 * t164 + t48 * t88 + t50 * t89;
t16 = t45 * t166 + t47 * t86 + t49 * t87;
t15 = (t142 * t103 + t52) * t103 + (t142 * t102 + t51) * t102 + t174;
t14 = -t44 * t164 - t100 * t38 + (-t103 * t168 + t52 * t99) * qJD(2);
t13 = t44 * t166 + t100 * t37 + (t102 * t168 - t51 * t99) * qJD(2);
t12 = t102 * t37 + t103 * t38 - t95 * t186 - t173;
t11 = (-t102 * t38 + t103 * t37) * t99 + t134 * t152;
t10 = t118 * t103 + t34 * t88 + t36 * t89 + t48 * t56 + t50 * t57;
t9 = t119 * t103 + t33 * t88 + t35 * t89 + t47 * t56 + t49 * t57;
t8 = t118 * t102 + t34 * t86 + t36 * t87 + t48 * t54 + t50 * t55;
t7 = t119 * t102 + t33 * t86 + t35 * t87 + t47 * t54 + t49 * t55;
t6 = (t132 * qJD(2) - t32) * t100 + (qJD(2) * t46 - t105 * t34 + t107 * t36 + (-t105 * t50 - t107 * t48) * qJD(4)) * t99;
t5 = (t133 * qJD(2) - t31) * t100 + (qJD(2) * t45 - t105 * t33 + t107 * t35 + (-t105 * t49 - t107 * t47) * qJD(4)) * t99;
t4 = t10 * t102 - t103 * t9;
t3 = t102 * t8 - t103 * t7;
t2 = -(t42 * t88 + t43 * t89 + t56 * t63 + t57 * t64) * t100 + (t9 * t102 + (t10 - t170) * t103) * t99 + (t23 * t99 + (t167 + (t19 - t169) * t103) * t100) * qJD(2);
t1 = -(t42 * t86 + t43 * t87 + t54 * t63 + t55 * t64) * t100 + (t8 * t103 + (t7 - t170) * t102) * t99 + (t22 * t99 + (t165 + (t16 - t169) * t102) * t100) * qJD(2);
t41 = [0; -m(3) * t53 + m(4) * t28 + m(5) * t12; 0.2e1 * m(4) * (t174 * t28 + (t103 * t117 + t145 * t59) * t103) - t103 * t3 + (t12 * t15 + t29 * t39 + t30 * t40) * t183 + 0.2e1 * m(3) * (qJD(2) * t97 - t53) * t196 * (rSges(3,1) * t108 - rSges(3,2) * t106) - t195 * t103 * t181 + (t4 + t194 * t182 + 0.2e1 * (t102 * t117 + t145 * t58) * m(4) + (-t195 * t102 + t194 * t103) * t103) * t102; 0; m(4) * (t102 * t59 - t103 * t58) + m(5) * (t102 * t30 - t103 * t29); 0; m(5) * t11; -t100 * (t102 * t6 - t103 * t5) / 0.2e1 + t2 * t180 - t103 * t1 / 0.2e1 + m(5) * (t11 * t15 + t12 * t24 + t13 * t40 + t14 * t39 + t26 * t30 + t27 * t29) + (t4 * t179 + t3 * t180) * t99 + (t99 * (t102 * t21 - t103 * t20) / 0.2e1 + ((t102 * t19 - t103 * t18) * t179 + (t102 * t17 - t103 * t16) * t180) * t100) * qJD(2); m(5) * (t102 * t13 - t103 * t14); (t11 * t24 + t13 * t26 + t14 * t27) * t183 + (-t100 * t23 + (t103 * t19 + t167) * t99) * t146 + t2 * t164 + (-t100 * t22 + (t102 * t16 + t165) * t99) * t147 + t1 * t166 + (-t100 * t25 + t135 * t99) * t158 - t100 * ((t170 + (t131 * t100 + t135) * qJD(2)) * t100 + (t6 * t103 + t5 * t102 - (qJD(2) * t62 - t105 * t42 + t107 * t43 + (-t105 * t64 - t107 * t63) * qJD(4)) * t100 + t25 * qJD(2)) * t99);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t41(1), t41(2), t41(4), t41(7); t41(2), t41(3), t41(5), t41(8); t41(4), t41(5), t41(6), t41(9); t41(7), t41(8), t41(9), t41(10);];
Mq = res;
