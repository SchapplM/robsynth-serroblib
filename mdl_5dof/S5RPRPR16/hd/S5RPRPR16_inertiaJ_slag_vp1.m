% Calculate joint inertia matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR16_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR16_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:48
% DurationCPUTime: 1.45s
% Computational Cost: add. (1163->218), mult. (2769->333), div. (0->0), fcn. (2813->6), ass. (0->113)
t167 = -Icges(4,4) - Icges(5,6);
t166 = Icges(4,1) + Icges(5,2);
t165 = Icges(4,2) + Icges(5,3);
t103 = cos(qJ(3));
t164 = t167 * t103;
t100 = sin(qJ(3));
t163 = t167 * t100;
t162 = Icges(5,1) + Icges(4,3);
t161 = t165 * t103 - t163;
t160 = t166 * t100 - t164;
t159 = (-Icges(5,5) + Icges(4,6)) * t103 + (-Icges(5,4) + Icges(4,5)) * t100;
t101 = sin(qJ(1));
t158 = -t101 / 0.2e1;
t146 = t101 / 0.2e1;
t104 = cos(qJ(1));
t157 = -t104 / 0.2e1;
t156 = t104 / 0.2e1;
t155 = t103 / 0.2e1;
t141 = rSges(5,2) * t100;
t105 = -t141 + (-rSges(5,3) - qJ(4)) * t103;
t135 = t100 * t104;
t88 = pkin(3) * t135;
t90 = t104 * qJ(2);
t144 = t88 + t90;
t148 = -pkin(1) - pkin(6);
t21 = t105 * t104 + (-rSges(5,1) + t148) * t101 + t144;
t143 = t104 * pkin(1) + t101 * qJ(2);
t128 = t104 * pkin(6) + t143;
t136 = t100 * t101;
t87 = pkin(3) * t136;
t22 = t104 * rSges(5,1) + t105 * t101 + t128 + t87;
t154 = m(5) * (t101 * t21 - t104 * t22);
t102 = cos(qJ(5));
t132 = t104 * t102;
t99 = sin(qJ(5));
t63 = -t101 * t99 + t103 * t132;
t133 = t103 * t104;
t64 = t101 * t102 + t99 * t133;
t123 = -t64 * rSges(6,1) - t63 * rSges(6,2);
t16 = (-qJ(4) * t103 + (rSges(6,3) + pkin(7)) * t100) * t104 + (-pkin(4) + t148) * t101 + t123 + t144;
t134 = t101 * t103;
t65 = -t102 * t134 - t104 * t99;
t66 = -t99 * t134 + t132;
t33 = t66 * rSges(6,1) + t65 * rSges(6,2) + rSges(6,3) * t136;
t67 = -qJ(4) * t134 + t87;
t94 = t104 * pkin(4);
t149 = -pkin(7) * t136 - t33 - t67 - t94;
t17 = t128 - t149;
t153 = m(6) * (t101 * t16 - t104 * t17);
t152 = t162 * t101 - t159 * t104;
t151 = t159 * t101 + t162 * t104;
t150 = (rSges(4,1) * t100 + rSges(4,2) * t103) * t104;
t96 = t101 ^ 2;
t98 = t104 ^ 2;
t80 = t103 * rSges(4,1) - t100 * rSges(4,2);
t147 = m(4) * t80;
t142 = t96 + t98;
t131 = m(5) / 0.2e1 + m(6) / 0.2e1;
t41 = Icges(6,3) * t103 + (Icges(6,5) * t99 + Icges(6,6) * t102) * t100;
t44 = Icges(6,6) * t103 + (Icges(6,4) * t99 + Icges(6,2) * t102) * t100;
t47 = Icges(6,5) * t103 + (Icges(6,1) * t99 + Icges(6,4) * t102) * t100;
t130 = t103 * t41 + (t102 * t44 + t47 * t99) * t100;
t129 = rSges(4,1) * t136 + rSges(4,2) * t134 + t104 * rSges(4,3);
t26 = Icges(6,5) * t64 + Icges(6,6) * t63 - Icges(6,3) * t135;
t28 = Icges(6,4) * t64 + Icges(6,2) * t63 - Icges(6,6) * t135;
t30 = Icges(6,1) * t64 + Icges(6,4) * t63 - Icges(6,5) * t135;
t10 = t103 * t26 + (t102 * t28 + t30 * t99) * t100;
t12 = -t41 * t135 + t63 * t44 + t64 * t47;
t127 = -t12 / 0.2e1 - t10 / 0.2e1;
t27 = Icges(6,5) * t66 + Icges(6,6) * t65 + Icges(6,3) * t136;
t29 = Icges(6,4) * t66 + Icges(6,2) * t65 + Icges(6,6) * t136;
t31 = Icges(6,1) * t66 + Icges(6,4) * t65 + Icges(6,5) * t136;
t11 = t103 * t27 + (t102 * t29 + t31 * t99) * t100;
t13 = t41 * t136 + t65 * t44 + t66 * t47;
t126 = t13 / 0.2e1 + t11 / 0.2e1;
t125 = Icges(4,5) * t155 - Icges(5,4) * t103 / 0.2e1 + (-Icges(4,6) / 0.2e1 + Icges(5,5) / 0.2e1) * t100;
t56 = t103 * rSges(6,3) + (rSges(6,1) * t99 + rSges(6,2) * t102) * t100;
t124 = pkin(7) * t103 + t56;
t121 = rSges(5,3) * t103 + t141;
t19 = t103 * t33 - t56 * t136;
t32 = -rSges(6,3) * t135 - t123;
t20 = -t103 * t32 - t56 * t135;
t115 = t20 * t101 - t19 * t104;
t78 = t103 * pkin(3) + t100 * qJ(4);
t68 = t101 * t78;
t24 = t124 * t101 + t68;
t25 = (-t124 - t78) * t104;
t113 = t24 * t101 - t25 * t104;
t79 = -t103 * rSges(5,2) + t100 * rSges(5,3);
t36 = t101 * t79 + t68;
t37 = (-t78 - t79) * t104;
t112 = t36 * t101 - t37 * t104;
t81 = t104 * rSges(2,1) - t101 * rSges(2,2);
t77 = -t101 * rSges(2,1) - t104 * rSges(2,2);
t59 = t104 * (qJ(4) * t133 - t88);
t58 = -t104 * rSges(3,2) + t101 * rSges(3,3) + t143;
t57 = t104 * rSges(3,3) + t90 + (rSges(3,2) - pkin(1)) * t101;
t35 = t128 + t129;
t34 = t90 + t150 + (-rSges(4,3) + t148) * t101;
t23 = -t101 * t129 + (t101 * rSges(4,3) - t150) * t104;
t18 = t130 * t103;
t15 = t59 + t121 * t98 + (t121 * t101 - t67) * t101;
t14 = (t101 * t32 + t104 * t33) * t100;
t9 = t59 + (-pkin(7) * t135 + t32) * t104 + (t94 + t149) * t101;
t8 = t27 * t136 + t65 * t29 + t66 * t31;
t7 = t26 * t136 + t65 * t28 + t66 * t30;
t6 = -t27 * t135 + t63 * t29 + t64 * t31;
t5 = -t26 * t135 + t63 * t28 + t64 * t30;
t4 = t7 * t101 + t8 * t104;
t3 = t5 * t101 + t6 * t104;
t2 = t13 * t103 + (t101 * t8 - t104 * t7) * t100;
t1 = t12 * t103 + (t101 * t6 - t104 * t5) * t100;
t38 = [Icges(3,1) + Icges(2,3) + m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + m(2) * (t77 ^ 2 + t81 ^ 2) + t130 + (t166 * t103 + t163) * t103 + (t165 * t100 + t164) * t100; t153 + t154 + m(4) * (t101 * t34 - t104 * t35) + m(3) * (t101 * t57 - t104 * t58); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t131) * t142; m(6) * (t24 * t16 + t25 * t17) + m(5) * (t36 * t21 + t37 * t22) + (-t35 * t147 + t125 * t104 + (Icges(5,4) * t157 + Icges(4,5) * t156 + t160 * t146) * t103 + t126) * t104 + (t34 * t147 + t125 * t101 + (Icges(5,4) * t158 + Icges(4,5) * t146 + t160 * t157) * t103 - t127) * t101 + ((Icges(5,5) * t156 + Icges(4,6) * t157 + t161 * t158) * t104 + (Icges(5,5) * t146 + Icges(4,6) * t158 + t161 * t156) * t101) * t100; m(5) * t112 + m(6) * t113 + t142 * t147; m(6) * (t24 ^ 2 + t25 ^ 2 + t9 ^ 2) + m(5) * (t15 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(4) * (t142 * t80 ^ 2 + t23 ^ 2) + (t151 * t98 + t4) * t104 + (t3 + t152 * t96 + (t151 * t101 + t152 * t104) * t104) * t101; 0.2e1 * (-t153 / 0.2e1 - t154 / 0.2e1) * t103; -0.2e1 * t131 * t142 * t103; m(6) * (t100 * t9 - t113 * t103) + m(5) * (t100 * t15 - t112 * t103); 0.2e1 * t131 * (t142 * t103 ^ 2 + t100 ^ 2); m(6) * (t20 * t16 + t19 * t17) + t18 + (t126 * t101 + t127 * t104) * t100; m(6) * t115; m(6) * (t14 * t9 + t19 * t25 + t20 * t24) + t2 * t156 + t1 * t146 + (t10 * t101 + t11 * t104) * t155 + (t4 * t146 + t3 * t157) * t100; m(6) * (t14 * t100 - t115 * t103); m(6) * (t14 ^ 2 + t19 ^ 2 + t20 ^ 2) + t103 * t18 + (t101 * t2 - t104 * t1 + t103 * (-t10 * t104 + t101 * t11)) * t100;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t38(1), t38(2), t38(4), t38(7), t38(11); t38(2), t38(3), t38(5), t38(8), t38(12); t38(4), t38(5), t38(6), t38(9), t38(13); t38(7), t38(8), t38(9), t38(10), t38(14); t38(11), t38(12), t38(13), t38(14), t38(15);];
Mq = res;
