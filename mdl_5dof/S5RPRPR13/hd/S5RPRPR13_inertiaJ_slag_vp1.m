% Calculate joint inertia matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR13_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR13_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:58
% EndTime: 2019-12-31 18:32:02
% DurationCPUTime: 1.50s
% Computational Cost: add. (2148->231), mult. (2759->354), div. (0->0), fcn. (2817->8), ass. (0->112)
t156 = -Icges(5,4) + Icges(4,5);
t155 = Icges(5,5) - Icges(4,6);
t154 = Icges(5,1) + Icges(4,3);
t93 = pkin(8) + qJ(3);
t88 = sin(t93);
t89 = cos(t93);
t153 = t155 * t88 + t156 * t89;
t100 = sin(qJ(1));
t152 = -t100 / 0.2e1;
t144 = t100 / 0.2e1;
t102 = cos(qJ(1));
t151 = -t102 / 0.2e1;
t150 = t102 / 0.2e1;
t149 = t154 * t100 + t153 * t102;
t148 = -t153 * t100 + t154 * t102;
t130 = t102 * t89;
t101 = cos(qJ(5));
t126 = t102 * t101;
t99 = sin(qJ(5));
t132 = t100 * t99;
t62 = t88 * t126 - t132;
t127 = t100 * t101;
t129 = t102 * t99;
t63 = t88 * t129 + t127;
t32 = t63 * rSges(6,1) + t62 * rSges(6,2) + rSges(6,3) * t130;
t147 = t100 * pkin(4) + pkin(7) * t130 + t32;
t94 = t100 ^ 2;
t95 = t102 ^ 2;
t146 = m(5) / 0.2e1;
t145 = m(6) / 0.2e1;
t134 = qJ(4) * t88;
t131 = t102 * t88;
t141 = pkin(3) * t130 + qJ(4) * t131;
t143 = t94 * (pkin(3) * t89 + t134) + t102 * t141;
t73 = t88 * pkin(3) - t89 * qJ(4);
t142 = t88 * rSges(5,2) + t89 * rSges(5,3) - t73;
t139 = t94 + t95;
t138 = Icges(4,4) * t88;
t137 = Icges(4,4) * t89;
t136 = Icges(5,6) * t88;
t135 = Icges(5,6) * t89;
t133 = t100 * t89;
t128 = rSges(3,3) + qJ(2);
t125 = t146 + t145;
t26 = Icges(6,5) * t63 + Icges(6,6) * t62 + Icges(6,3) * t130;
t28 = Icges(6,4) * t63 + Icges(6,2) * t62 + Icges(6,6) * t130;
t30 = Icges(6,1) * t63 + Icges(6,4) * t62 + Icges(6,5) * t130;
t10 = t88 * t26 + (-t101 * t28 - t30 * t99) * t89;
t41 = Icges(6,3) * t88 + (-Icges(6,5) * t99 - Icges(6,6) * t101) * t89;
t42 = Icges(6,6) * t88 + (-Icges(6,4) * t99 - Icges(6,2) * t101) * t89;
t43 = Icges(6,5) * t88 + (-Icges(6,1) * t99 - Icges(6,4) * t101) * t89;
t12 = t41 * t130 + t62 * t42 + t63 * t43;
t124 = t10 / 0.2e1 + t12 / 0.2e1;
t64 = t88 * t127 + t129;
t65 = t88 * t132 - t126;
t27 = Icges(6,5) * t65 + Icges(6,6) * t64 + Icges(6,3) * t133;
t29 = Icges(6,4) * t65 + Icges(6,2) * t64 + Icges(6,6) * t133;
t31 = Icges(6,1) * t65 + Icges(6,4) * t64 + Icges(6,5) * t133;
t11 = t88 * t27 + (-t101 * t29 - t31 * t99) * t89;
t13 = t41 * t133 + t64 * t42 + t65 * t43;
t123 = t11 / 0.2e1 + t13 / 0.2e1;
t97 = cos(pkin(8));
t86 = t97 * pkin(2) + pkin(1);
t98 = -pkin(6) - qJ(2);
t122 = -t100 * t98 + t102 * t86;
t44 = t88 * rSges(6,3) + (-rSges(6,1) * t99 - rSges(6,2) * t101) * t89;
t121 = -pkin(7) * t88 - t44 - t73;
t119 = rSges(4,1) * t89 - rSges(4,2) * t88;
t118 = -t65 * rSges(6,1) - t64 * rSges(6,2);
t113 = -t101 * t42 - t99 * t43;
t112 = Icges(4,1) * t89 - t138;
t111 = -Icges(4,2) * t88 + t137;
t108 = -Icges(5,2) * t89 + t136;
t107 = Icges(5,3) * t88 - t135;
t106 = rSges(4,1) * t130 - rSges(4,2) * t131 + t100 * rSges(4,3);
t105 = t100 * rSges(5,1) - rSges(5,2) * t130 + rSges(5,3) * t131;
t104 = t122 + t141;
t96 = sin(pkin(8));
t103 = rSges(3,1) * t97 - rSges(3,2) * t96 + pkin(1);
t77 = t102 * rSges(2,1) - t100 * rSges(2,2);
t76 = -t100 * rSges(2,1) - t102 * rSges(2,2);
t75 = t88 * rSges(4,1) + t89 * rSges(4,2);
t40 = t128 * t100 + t103 * t102;
t39 = -t103 * t100 + t128 * t102;
t38 = t88 * t41;
t37 = t142 * t102;
t36 = t142 * t100;
t35 = t106 + t122;
t34 = (rSges(4,3) - t98) * t102 + (-t119 - t86) * t100;
t33 = rSges(6,3) * t133 - t118;
t25 = t102 * t106 + (-t102 * rSges(4,3) + t119 * t100) * t100;
t24 = t121 * t102;
t23 = t121 * t100;
t22 = t104 + t105;
t21 = (rSges(5,1) - t98) * t102 + (-t86 + (rSges(5,2) - pkin(3)) * t89 + (-rSges(5,3) - qJ(4)) * t88) * t100;
t20 = -t44 * t130 + t88 * t32;
t19 = t44 * t133 - t88 * t33;
t18 = t104 + t147;
t17 = (pkin(4) - t98) * t102 + (-t134 - t86 + (-rSges(6,3) - pkin(3) - pkin(7)) * t89) * t100 + t118;
t16 = t102 * t105 + (-t102 * rSges(5,1) + (-rSges(5,2) * t89 + rSges(5,3) * t88) * t100) * t100 + t143;
t15 = (t113 * t89 + t38) * t88;
t14 = (-t100 * t32 + t102 * t33) * t89;
t9 = t147 * t102 + (-t102 * pkin(4) + pkin(7) * t133 + t33) * t100 + t143;
t8 = t27 * t133 + t64 * t29 + t65 * t31;
t7 = t26 * t133 + t64 * t28 + t65 * t30;
t6 = t27 * t130 + t62 * t29 + t63 * t31;
t5 = t26 * t130 + t62 * t28 + t63 * t30;
t4 = t7 * t100 - t8 * t102;
t3 = t5 * t100 - t6 * t102;
t2 = t13 * t88 + (t100 * t8 + t102 * t7) * t89;
t1 = t12 * t88 + (t100 * t6 + t102 * t5) * t89;
t45 = [Icges(3,2) * t97 ^ 2 + Icges(2,3) + t38 + (Icges(3,1) * t96 + 0.2e1 * Icges(3,4) * t97) * t96 + m(6) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + m(2) * (t76 ^ 2 + t77 ^ 2) + (t113 + t136 + t138 + (Icges(4,2) + Icges(5,3)) * t89) * t89 + (t135 + t137 + (Icges(4,1) + Icges(5,2)) * t88) * t88; m(6) * (t100 * t17 - t102 * t18) + m(5) * (t100 * t21 - t102 * t22) + m(4) * (t100 * t34 - t102 * t35) + m(3) * (t100 * t39 - t102 * t40); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t125) * t139; ((Icges(5,5) * t151 + Icges(4,6) * t150 + t107 * t144 + t111 * t152) * t89 + (Icges(5,4) * t151 + Icges(4,5) * t150 + t108 * t144 + t112 * t152) * t88 - t123) * t102 + ((Icges(5,5) * t152 + Icges(4,6) * t144 + t107 * t151 + t111 * t150) * t89 + (Icges(5,4) * t152 + Icges(4,5) * t144 + t108 * t151 + t112 * t150) * t88 + t124) * t100 + m(6) * (t24 * t17 + t23 * t18) + m(5) * (t37 * t21 + t36 * t22) + m(4) * (-t100 * t35 - t102 * t34) * t75 + (-t155 * t89 + t156 * t88) * (t95 / 0.2e1 + t94 / 0.2e1); m(5) * (t37 * t100 - t36 * t102) + m(6) * (t24 * t100 - t23 * t102); m(6) * (t23 ^ 2 + t24 ^ 2 + t9 ^ 2) + m(5) * (t16 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(4) * (t139 * t75 ^ 2 + t25 ^ 2) + (t148 * t95 - t4) * t102 + (t3 + t149 * t94 + (t148 * t100 + t149 * t102) * t102) * t100; 0.2e1 * ((t100 * t18 + t102 * t17) * t145 + (t100 * t22 + t102 * t21) * t146) * t88; 0; m(6) * (-t89 * t9 + (t100 * t23 + t102 * t24) * t88) + m(5) * (-t89 * t16 + (t100 * t36 + t102 * t37) * t88); 0.2e1 * t125 * (t139 * t88 ^ 2 + t89 ^ 2); m(6) * (t19 * t17 + t20 * t18) + t15 + (t123 * t100 + t124 * t102) * t89; m(6) * (t19 * t100 - t20 * t102); t88 * (t10 * t100 - t11 * t102) / 0.2e1 + t1 * t144 + m(6) * (t14 * t9 + t19 * t24 + t20 * t23) + t2 * t151 + (t4 * t144 + t3 * t150) * t89; m(6) * (-t14 * t89 + (t100 * t20 + t102 * t19) * t88); m(6) * (t14 ^ 2 + t19 ^ 2 + t20 ^ 2) + t88 * t15 + (t102 * t1 + t100 * t2 + t88 * (t10 * t102 + t100 * t11)) * t89;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t45(1), t45(2), t45(4), t45(7), t45(11); t45(2), t45(3), t45(5), t45(8), t45(12); t45(4), t45(5), t45(6), t45(9), t45(13); t45(7), t45(8), t45(9), t45(10), t45(14); t45(11), t45(12), t45(13), t45(14), t45(15);];
Mq = res;
