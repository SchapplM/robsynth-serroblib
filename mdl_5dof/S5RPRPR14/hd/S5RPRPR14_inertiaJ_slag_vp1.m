% Calculate joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR14_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR14_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:34:35
% DurationCPUTime: 1.28s
% Computational Cost: add. (1888->215), mult. (2609->324), div. (0->0), fcn. (2662->8), ass. (0->111)
t97 = qJ(3) + pkin(8);
t91 = cos(t97);
t165 = Icges(5,5) * t91;
t90 = sin(t97);
t164 = Icges(5,6) * t90;
t105 = cos(qJ(3));
t163 = Icges(4,5) * t105;
t102 = sin(qJ(3));
t162 = Icges(4,6) * t102;
t161 = Icges(4,3) + Icges(5,3);
t160 = Icges(4,5) * t102 + Icges(5,5) * t90 + Icges(4,6) * t105 + Icges(5,6) * t91;
t159 = t163 / 0.2e1 + t165 / 0.2e1 - t162 / 0.2e1 - t164 / 0.2e1;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t158 = t161 * t103 - t160 * t106;
t157 = t160 * t103 + t161 * t106;
t156 = (rSges(4,1) * t102 + rSges(4,2) * t105) * t106;
t98 = t103 ^ 2;
t99 = t106 ^ 2;
t153 = pkin(4) * t90;
t150 = t106 / 0.2e1;
t149 = pkin(3) * t105;
t104 = cos(qJ(5));
t101 = sin(qJ(5));
t38 = Icges(6,3) * t90 + (Icges(6,5) * t104 - Icges(6,6) * t101) * t91;
t40 = Icges(6,5) * t90 + (Icges(6,1) * t104 - Icges(6,4) * t101) * t91;
t148 = t91 * t104 * t40 + t90 * t38;
t41 = t90 * rSges(6,3) + (rSges(6,1) * t104 - rSges(6,2) * t101) * t91;
t147 = t91 * pkin(4) + t90 * pkin(7) + t41;
t129 = t106 * t104;
t133 = t103 * t101;
t62 = -t90 * t133 + t129;
t130 = t106 * t101;
t132 = t103 * t104;
t63 = t90 * t132 + t130;
t146 = t63 * rSges(6,1) + t62 * rSges(6,2);
t83 = t98 + t99;
t145 = (m(5) + m(6)) * t83;
t100 = -qJ(4) - pkin(6);
t144 = t106 * t102 * pkin(3) + t103 * t100;
t143 = t106 * pkin(1) + t103 * qJ(2);
t39 = Icges(6,6) * t90 + (Icges(6,4) * t104 - Icges(6,2) * t101) * t91;
t140 = t101 * t39;
t139 = t103 * t90;
t138 = t103 * t91;
t94 = t106 * rSges(5,3);
t137 = t106 * t91;
t134 = t102 * t103;
t131 = t103 * t105;
t128 = -rSges(5,1) * t139 - rSges(5,2) * t138 - t94;
t127 = rSges(4,1) * t134 + rSges(4,2) * t131 + t106 * rSges(4,3);
t93 = t106 * qJ(2);
t126 = t93 + t144;
t23 = Icges(6,5) * t63 + Icges(6,6) * t62 - Icges(6,3) * t138;
t25 = Icges(6,4) * t63 + Icges(6,2) * t62 - Icges(6,6) * t138;
t27 = Icges(6,1) * t63 + Icges(6,4) * t62 - Icges(6,5) * t138;
t10 = t90 * t23 + (-t101 * t25 + t104 * t27) * t91;
t12 = -t38 * t138 + t62 * t39 + t63 * t40;
t125 = -t12 / 0.2e1 - t10 / 0.2e1;
t64 = t90 * t130 + t132;
t65 = -t90 * t129 + t133;
t24 = Icges(6,5) * t65 + Icges(6,6) * t64 + Icges(6,3) * t137;
t26 = Icges(6,4) * t65 + Icges(6,2) * t64 + Icges(6,6) * t137;
t28 = Icges(6,1) * t65 + Icges(6,4) * t64 + Icges(6,5) * t137;
t11 = t90 * t24 + (-t101 * t26 + t104 * t28) * t91;
t13 = t38 * t137 + t64 * t39 + t65 * t40;
t124 = t13 / 0.2e1 + t11 / 0.2e1;
t123 = (-rSges(6,3) - pkin(7)) * t91;
t121 = rSges(5,1) * t90 + rSges(5,2) * t91;
t120 = -t65 * rSges(6,1) - t64 * rSges(6,2);
t86 = pkin(3) * t134;
t108 = -t106 * t100 + t143 + t86;
t35 = t93 + t156 + (-rSges(4,3) - pkin(1) - pkin(6)) * t103;
t36 = t106 * pkin(6) + t127 + t143;
t107 = m(4) * (t103 * t35 - t106 * t36);
t87 = pkin(3) * t131;
t82 = pkin(4) * t139;
t79 = t106 * rSges(2,1) - t103 * rSges(2,2);
t78 = t105 * rSges(4,1) - t102 * rSges(4,2);
t77 = -t103 * rSges(2,1) - t106 * rSges(2,2);
t69 = t91 * rSges(5,1) - t90 * rSges(5,2);
t61 = t86 + (-pkin(6) - t100) * t106;
t60 = -t106 * rSges(3,2) + t103 * rSges(3,3) + t143;
t59 = t106 * rSges(3,3) + t93 + (rSges(3,2) - pkin(1)) * t103;
t50 = t106 * (-t103 * pkin(6) - t144);
t43 = (-t69 - t149) * t106;
t42 = t103 * t69 + t87;
t33 = t108 - t128;
t32 = t121 * t106 + (-rSges(5,3) - pkin(1)) * t103 + t126;
t31 = -t103 * t127 + (t103 * rSges(4,3) - t156) * t106;
t30 = rSges(6,3) * t137 - t120;
t29 = -rSges(6,3) * t138 + t146;
t22 = (-t147 - t149) * t106;
t21 = t147 * t103 + t87;
t20 = t103 * t123 + t108 + t146 + t82;
t19 = -t103 * pkin(1) + (t123 + t153) * t106 + t120 + t126;
t18 = t41 * t137 - t90 * t30;
t17 = t41 * t138 + t90 * t29;
t16 = t50 - t121 * t99 + (t128 - t61 + t94) * t103;
t15 = (-t91 * t140 + t148) * t90;
t14 = (-t103 * t30 - t106 * t29) * t91;
t9 = t50 + (t30 + (pkin(7) * t91 - t153) * t106) * t106 + (pkin(7) * t138 - t29 - t61 - t82) * t103;
t8 = t24 * t137 + t64 * t26 + t65 * t28;
t7 = t23 * t137 + t64 * t25 + t65 * t27;
t6 = -t24 * t138 + t62 * t26 + t63 * t28;
t5 = -t23 * t138 + t62 * t25 + t63 * t27;
t4 = t8 * t103 + t7 * t106;
t3 = t6 * t103 + t5 * t106;
t2 = t13 * t90 + (-t103 * t7 + t106 * t8) * t91;
t1 = t12 * t90 + (-t103 * t5 + t106 * t6) * t91;
t34 = [Icges(4,1) * t105 ^ 2 + Icges(3,1) + Icges(2,3) + (Icges(5,1) * t91 - t140) * t91 + m(6) * (t19 ^ 2 + t20 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t59 ^ 2 + t60 ^ 2) + m(2) * (t77 ^ 2 + t79 ^ 2) + t148 + (-0.2e1 * Icges(5,4) * t91 + Icges(5,2) * t90) * t90 + (-0.2e1 * Icges(4,4) * t105 + Icges(4,2) * t102) * t102; m(6) * (t103 * t19 - t106 * t20) + m(5) * (t103 * t32 - t106 * t33) + t107 + m(3) * (t103 * t59 - t106 * t60); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t83 + t145; m(6) * (t21 * t19 + t22 * t20) + m(5) * (t42 * t32 + t43 * t33) + t78 * t107 + (-t162 + t163 - t164 + t165) * (t99 / 0.2e1 + t98 / 0.2e1) + (t159 * t106 - t125) * t106 + (t159 * t103 + t124) * t103; m(5) * (t42 * t103 - t43 * t106) + m(6) * (t21 * t103 - t22 * t106) + m(4) * t83 * t78; m(6) * (t21 ^ 2 + t22 ^ 2 + t9 ^ 2) + m(5) * (t16 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(4) * (t83 * t78 ^ 2 + t31 ^ 2) + (t158 * t98 + t4) * t103 + (t3 + t157 * t99 + (t157 * t103 + t158 * t106) * t103) * t106; m(6) * (t103 * t20 + t106 * t19) + m(5) * (t103 * t33 + t106 * t32); 0; m(6) * (t103 * t22 + t106 * t21) + m(5) * (t103 * t43 + t106 * t42); t145; m(6) * (t17 * t20 + t18 * t19) + t15 + (t125 * t103 + t124 * t106) * t91; m(6) * (t18 * t103 - t17 * t106); m(6) * (t14 * t9 + t17 * t22 + t18 * t21) + t90 * (t10 * t106 + t11 * t103) / 0.2e1 + t1 * t150 + t103 * t2 / 0.2e1 + (t4 * t150 - t103 * t3 / 0.2e1) * t91; m(6) * (t17 * t103 + t18 * t106); m(6) * (t14 ^ 2 + t17 ^ 2 + t18 ^ 2) + t90 * t15 + (-t103 * t1 + t106 * t2 + t90 * (-t10 * t103 + t106 * t11)) * t91;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t34(1), t34(2), t34(4), t34(7), t34(11); t34(2), t34(3), t34(5), t34(8), t34(12); t34(4), t34(5), t34(6), t34(9), t34(13); t34(7), t34(8), t34(9), t34(10), t34(14); t34(11), t34(12), t34(13), t34(14), t34(15);];
Mq = res;
