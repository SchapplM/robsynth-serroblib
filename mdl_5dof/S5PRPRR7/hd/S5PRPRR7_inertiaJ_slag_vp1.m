% Calculate joint inertia matrix for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:45
% EndTime: 2019-12-05 15:59:49
% DurationCPUTime: 1.50s
% Computational Cost: add. (3170->232), mult. (5191->381), div. (0->0), fcn. (5625->8), ass. (0->122)
t165 = Icges(4,1) + Icges(3,3);
t117 = sin(qJ(2));
t119 = cos(qJ(2));
t164 = (-Icges(4,4) + Icges(3,5)) * t119 + (Icges(4,5) - Icges(3,6)) * t117;
t114 = sin(pkin(8));
t110 = t114 ^ 2;
t115 = cos(pkin(8));
t111 = t115 ^ 2;
t142 = t110 + t111;
t163 = t165 * t114 + t164 * t115;
t162 = -t164 * t114 + t165 * t115;
t161 = t117 ^ 2;
t160 = t114 / 0.2e1;
t159 = -t115 / 0.2e1;
t158 = t117 / 0.2e1;
t118 = cos(qJ(4));
t157 = pkin(4) * t118;
t146 = t115 * t119;
t113 = qJ(4) + qJ(5);
t108 = sin(t113);
t109 = cos(t113);
t147 = t115 * t117;
t85 = -t108 * t114 + t109 * t147;
t86 = t108 * t147 + t109 * t114;
t51 = rSges(6,1) * t86 + rSges(6,2) * t85 + rSges(6,3) * t146;
t116 = sin(qJ(4));
t145 = t116 * t117;
t121 = pkin(4) * t145 + pkin(7) * t119;
t64 = t114 * t157 + t115 * t121;
t155 = t51 + t64;
t148 = t114 * t119;
t149 = t114 * t117;
t87 = t108 * t115 + t109 * t149;
t88 = t108 * t149 - t109 * t115;
t52 = rSges(6,1) * t88 + rSges(6,2) * t87 + rSges(6,3) * t148;
t65 = t114 * t121 - t115 * t157;
t154 = t52 + t65;
t72 = rSges(6,3) * t117 + (-rSges(6,1) * t108 - rSges(6,2) * t109) * t119;
t95 = -t119 * t116 * pkin(4) + pkin(7) * t117;
t153 = -t72 - t95;
t152 = t142 * (pkin(2) * t119 + qJ(3) * t117);
t69 = Icges(6,3) * t117 + (-Icges(6,5) * t108 - Icges(6,6) * t109) * t119;
t151 = t117 * t69;
t89 = Icges(5,3) * t117 + (-Icges(5,5) * t116 - Icges(5,6) * t118) * t119;
t150 = t117 * t89;
t144 = t117 * t118;
t104 = pkin(2) * t117 - qJ(3) * t119;
t143 = rSges(4,2) * t117 + rSges(4,3) * t119 - t104;
t141 = -m(4) - m(5) - m(6);
t45 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t146;
t47 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t146;
t49 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t146;
t23 = t117 * t45 + (-t108 * t49 - t109 * t47) * t119;
t46 = Icges(6,5) * t88 + Icges(6,6) * t87 + Icges(6,3) * t148;
t48 = Icges(6,4) * t88 + Icges(6,2) * t87 + Icges(6,6) * t148;
t50 = Icges(6,1) * t88 + Icges(6,4) * t87 + Icges(6,5) * t148;
t24 = t117 * t46 + (-t108 * t50 - t109 * t48) * t119;
t19 = t146 * t45 + t47 * t85 + t49 * t86;
t20 = t146 * t46 + t48 * t85 + t50 * t86;
t70 = Icges(6,6) * t117 + (-Icges(6,4) * t108 - Icges(6,2) * t109) * t119;
t71 = Icges(6,5) * t117 + (-Icges(6,1) * t108 - Icges(6,4) * t109) * t119;
t5 = (t70 * t85 + t71 * t86) * t117 + (t20 * t114 + (t19 + t151) * t115) * t119;
t21 = t148 * t45 + t47 * t87 + t49 * t88;
t22 = t148 * t46 + t48 * t87 + t50 * t88;
t6 = (t70 * t87 + t71 * t88) * t117 + (t21 * t115 + (t22 + t151) * t114) * t119;
t140 = t5 * t146 + t6 * t148 + t117 * (t161 * t69 + (t23 * t115 + t24 * t114 + (-t108 * t71 - t109 * t70) * t117) * t119);
t139 = -t117 * pkin(6) - t104;
t138 = t152 + (t114 * t148 + t115 * t146) * pkin(6);
t12 = t114 * t19 - t115 * t20;
t13 = t114 * t21 - t115 * t22;
t137 = t13 * t148 / 0.2e1 + t6 * t159 + t5 * t160 + t12 * t146 / 0.2e1 + (t114 * t23 - t115 * t24) * t158;
t92 = rSges(5,3) * t117 + (-rSges(5,1) * t116 - rSges(5,2) * t118) * t119;
t136 = t139 - t92;
t128 = t139 + t153;
t106 = rSges(3,1) * t117 + rSges(3,2) * t119;
t102 = t114 * t145 - t115 * t118;
t101 = t114 * t144 + t115 * t116;
t100 = t114 * t118 + t115 * t145;
t99 = -t114 * t116 + t115 * t144;
t91 = Icges(5,5) * t117 + (-Icges(5,1) * t116 - Icges(5,4) * t118) * t119;
t90 = Icges(5,6) * t117 + (-Icges(5,4) * t116 - Icges(5,2) * t118) * t119;
t68 = t143 * t115;
t67 = t143 * t114;
t66 = t72 * t148;
t63 = rSges(5,1) * t102 + rSges(5,2) * t101 + rSges(5,3) * t148;
t62 = rSges(5,1) * t100 + rSges(5,2) * t99 + rSges(5,3) * t146;
t61 = t136 * t115;
t60 = t136 * t114;
t59 = Icges(5,1) * t102 + Icges(5,4) * t101 + Icges(5,5) * t148;
t58 = Icges(5,1) * t100 + Icges(5,4) * t99 + Icges(5,5) * t146;
t57 = Icges(5,4) * t102 + Icges(5,2) * t101 + Icges(5,6) * t148;
t56 = Icges(5,4) * t100 + Icges(5,2) * t99 + Icges(5,6) * t146;
t55 = Icges(5,5) * t102 + Icges(5,6) * t101 + Icges(5,3) * t148;
t54 = Icges(5,5) * t100 + Icges(5,6) * t99 + Icges(5,3) * t146;
t53 = t142 * (rSges(3,1) * t119 - rSges(3,2) * t117);
t44 = t117 * t51;
t43 = t52 * t146;
t42 = t128 * t115;
t41 = t128 * t114;
t40 = t117 * t62 - t146 * t92;
t39 = -t117 * t63 + t148 * t92;
t38 = -t146 * t72 + t44;
t37 = -t117 * t52 + t66;
t36 = t152 + t142 * (-rSges(4,2) * t119 + rSges(4,3) * t117);
t35 = (-t114 * t62 + t115 * t63) * t119;
t34 = -t148 * t51 + t43;
t33 = t117 * t55 + (-t116 * t59 - t118 * t57) * t119;
t32 = t117 * t54 + (-t116 * t58 - t118 * t56) * t119;
t31 = t117 * t64 + t146 * t153 + t44;
t30 = -t117 * t154 + t148 * t95 + t66;
t29 = t114 * t63 + t115 * t62 + t138;
t28 = t101 * t57 + t102 * t59 + t148 * t55;
t27 = t101 * t56 + t102 * t58 + t148 * t54;
t26 = t100 * t59 + t146 * t55 + t57 * t99;
t25 = t100 * t58 + t146 * t54 + t56 * t99;
t18 = t43 + (-t114 * t155 + t115 * t65) * t119;
t17 = t114 * t154 + t115 * t155 + t138;
t16 = t114 * t27 - t115 * t28;
t15 = t114 * t25 - t115 * t26;
t8 = (t101 * t90 + t102 * t91) * t117 + (t27 * t115 + (t28 + t150) * t114) * t119;
t7 = (t100 * t91 + t90 * t99) * t117 + (t26 * t114 + (t25 + t150) * t115) * t119;
t1 = [m(2) + m(3) - t141; m(3) * t53 + m(4) * t36 + m(5) * t29 + m(6) * t17; m(6) * (t17 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t29 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t36 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(3) * (t106 ^ 2 * t142 + t53 ^ 2) + (t162 * t111 - t13 - t16) * t115 + (t12 + t15 + t163 * t110 + (t162 * t114 + t163 * t115) * t115) * t114; t141 * t119; m(6) * (-t119 * t17 + (t114 * t41 + t115 * t42) * t117) + m(5) * (-t119 * t29 + (t114 * t60 + t115 * t61) * t117) + m(4) * (-t119 * t36 + (t114 * t67 + t115 * t68) * t117); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t119 ^ 2 + t142 * t161); m(5) * t35 + m(6) * t18; t7 * t160 + t8 * t159 + (t114 * t32 - t115 * t33) * t158 + (t16 * t160 + t115 * t15 / 0.2e1) * t119 + m(6) * (t17 * t18 + t30 * t42 + t31 * t41) + m(5) * (t29 * t35 + t39 * t61 + t40 * t60) + t137; m(5) * (-t119 * t35 + (t114 * t40 + t115 * t39) * t117) + m(6) * (-t119 * t18 + (t114 * t31 + t115 * t30) * t117); t8 * t148 + t7 * t146 + m(6) * (t18 ^ 2 + t30 ^ 2 + t31 ^ 2) + t117 * (t161 * t89 + (t32 * t115 + t33 * t114 + (-t116 * t91 - t118 * t90) * t117) * t119) + m(5) * (t35 ^ 2 + t39 ^ 2 + t40 ^ 2) + t140; m(6) * t34; m(6) * (t17 * t34 + t37 * t42 + t38 * t41) + t137; m(6) * (-t119 * t34 + (t114 * t38 + t115 * t37) * t117); m(6) * (t18 * t34 + t30 * t37 + t31 * t38) + t140; m(6) * (t34 ^ 2 + t37 ^ 2 + t38 ^ 2) + t140;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
