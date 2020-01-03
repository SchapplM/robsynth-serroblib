% Calculate joint inertia matrix for
% S5RPRPR9
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:23:47
% DurationCPUTime: 1.31s
% Computational Cost: add. (2226->216), mult. (2637->326), div. (0->0), fcn. (2717->8), ass. (0->118)
t98 = sin(qJ(3));
t156 = t98 / 0.2e1;
t167 = Icges(5,1) + Icges(4,3);
t101 = cos(qJ(3));
t166 = (Icges(5,5) - Icges(4,6)) * t98 + (-Icges(5,4) + Icges(4,5)) * t101;
t165 = Icges(4,6) / 0.2e1 - Icges(5,5) / 0.2e1;
t164 = -Icges(5,4) * t98 + 0.2e1 * Icges(4,5) * t156 + t165 * t101;
t95 = qJ(1) + pkin(8);
t92 = sin(t95);
t163 = -t92 / 0.2e1;
t150 = t92 / 0.2e1;
t93 = cos(t95);
t162 = -t93 / 0.2e1;
t161 = t93 / 0.2e1;
t155 = -t166 * t92 + t167 * t93;
t154 = t166 * t93 + t167 * t92;
t132 = t101 * t93;
t100 = cos(qJ(5));
t130 = t93 * t100;
t97 = sin(qJ(5));
t56 = t98 * t130 - t92 * t97;
t131 = t92 * t100;
t142 = t97 * t98;
t57 = t93 * t142 + t131;
t30 = t57 * rSges(6,1) + t56 * rSges(6,2) + rSges(6,3) * t132;
t153 = t92 * pkin(4) + pkin(7) * t132 + t30;
t90 = t92 ^ 2;
t91 = t93 ^ 2;
t152 = m(5) / 0.2e1;
t151 = m(6) / 0.2e1;
t149 = -m(5) - m(6);
t75 = t98 * rSges(4,1) + t101 * rSges(4,2);
t148 = m(4) * t75;
t147 = t93 * pkin(4);
t99 = sin(qJ(1));
t146 = t99 * pkin(1);
t145 = t93 * rSges(5,1);
t144 = t93 * rSges(4,3);
t143 = t93 * t98;
t134 = qJ(4) * t98;
t139 = pkin(3) * t132 + t93 * t134;
t141 = t90 * (pkin(3) * t101 + t134) + t93 * t139;
t73 = t98 * pkin(3) - t101 * qJ(4);
t140 = t98 * rSges(5,2) + t101 * rSges(5,3) - t73;
t137 = t90 + t91;
t136 = Icges(4,4) * t98;
t135 = Icges(5,6) * t98;
t133 = t101 * t92;
t129 = Icges(4,4) * t101;
t128 = Icges(6,5) * t101;
t127 = Icges(5,6) * t101;
t126 = Icges(6,6) * t101;
t125 = Icges(6,3) * t101;
t102 = cos(qJ(1));
t94 = t102 * pkin(1);
t124 = t93 * pkin(2) + t92 * pkin(6) + t94;
t58 = t98 * t131 + t93 * t97;
t59 = t92 * t142 - t130;
t25 = Icges(6,5) * t59 + Icges(6,6) * t58 + t92 * t125;
t27 = Icges(6,4) * t59 + Icges(6,2) * t58 + t92 * t126;
t29 = Icges(6,1) * t59 + Icges(6,4) * t58 + t92 * t128;
t11 = t98 * t25 + (-t100 * t27 - t29 * t97) * t101;
t60 = Icges(6,3) * t98 + (-Icges(6,5) * t97 - Icges(6,6) * t100) * t101;
t61 = Icges(6,6) * t98 + (-Icges(6,4) * t97 - Icges(6,2) * t100) * t101;
t62 = Icges(6,5) * t98 + (-Icges(6,1) * t97 - Icges(6,4) * t100) * t101;
t15 = t60 * t133 + t58 * t61 + t59 * t62;
t123 = t11 / 0.2e1 + t15 / 0.2e1;
t24 = Icges(6,5) * t57 + Icges(6,6) * t56 + t93 * t125;
t26 = Icges(6,4) * t57 + Icges(6,2) * t56 + t93 * t126;
t28 = Icges(6,1) * t57 + Icges(6,4) * t56 + t93 * t128;
t10 = t98 * t24 + (-t100 * t26 - t28 * t97) * t101;
t14 = t60 * t132 + t56 * t61 + t57 * t62;
t122 = t14 / 0.2e1 + t10 / 0.2e1;
t120 = t93 * pkin(6) - t146;
t63 = t98 * rSges(6,3) + (-rSges(6,1) * t97 - rSges(6,2) * t100) * t101;
t119 = -pkin(7) * t98 - t63 - t73;
t118 = t124 + t139;
t117 = -t59 * rSges(6,1) - t58 * rSges(6,2);
t116 = rSges(4,1) * t101 - rSges(4,2) * t98;
t115 = -t100 * t61 - t97 * t62;
t110 = rSges(4,1) * t132 - rSges(4,2) * t143 + t92 * rSges(4,3);
t109 = Icges(4,1) * t101 - t136;
t108 = -Icges(4,2) * t98 + t129;
t105 = -Icges(5,2) * t101 + t135;
t104 = Icges(5,3) * t98 - t127;
t103 = t92 * rSges(5,1) - rSges(5,2) * t132 + rSges(5,3) * t143;
t77 = t102 * rSges(2,1) - t99 * rSges(2,2);
t76 = -t99 * rSges(2,1) - t102 * rSges(2,2);
t65 = t93 * rSges(3,1) - t92 * rSges(3,2) + t94;
t64 = -t92 * rSges(3,1) - t93 * rSges(3,2) - t146;
t52 = t98 * t60;
t37 = t140 * t93;
t36 = t140 * t92;
t35 = t110 + t124;
t34 = t144 + (-pkin(2) - t116) * t92 + t120;
t33 = t119 * t93;
t32 = t119 * t92;
t31 = rSges(6,3) * t133 - t117;
t23 = t103 + t118;
t22 = t145 + (-pkin(2) + (-rSges(5,3) - qJ(4)) * t98 + (rSges(5,2) - pkin(3)) * t101) * t92 + t120;
t21 = t93 * t110 + (t116 * t92 - t144) * t92;
t20 = (t115 * t101 + t52) * t98;
t19 = -t63 * t132 + t98 * t30;
t18 = t63 * t133 - t98 * t31;
t17 = t118 + t153;
t16 = t147 + (-t134 - pkin(2) + (-rSges(6,3) - pkin(3) - pkin(7)) * t101) * t92 + t117 + t120;
t13 = t93 * t103 + (-t145 + (-rSges(5,2) * t101 + rSges(5,3) * t98) * t92) * t92 + t141;
t12 = (-t30 * t92 + t31 * t93) * t101;
t9 = t25 * t133 + t58 * t27 + t59 * t29;
t8 = t24 * t133 + t58 * t26 + t59 * t28;
t7 = t25 * t132 + t56 * t27 + t57 * t29;
t6 = t24 * t132 + t56 * t26 + t57 * t28;
t5 = t153 * t93 + (pkin(7) * t133 - t147 + t31) * t92 + t141;
t4 = t8 * t92 - t9 * t93;
t3 = t6 * t92 - t7 * t93;
t2 = t15 * t98 + (t8 * t93 + t9 * t92) * t101;
t1 = t14 * t98 + (t6 * t93 + t7 * t92) * t101;
t38 = [Icges(2,3) + Icges(3,3) + t52 + m(6) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t64 ^ 2 + t65 ^ 2) + m(2) * (t76 ^ 2 + t77 ^ 2) + (t127 + t129 + (Icges(4,1) + Icges(5,2)) * t98) * t98 + (t115 + t135 + t136 + (Icges(4,2) + Icges(5,3)) * t101) * t101; 0; m(3) + m(4) - t149; m(6) * (t33 * t16 + t32 * t17) + m(5) * (t37 * t22 + t36 * t23) + (-t34 * t148 + (t105 * t150 + t109 * t163) * t98 - t123 + t164 * t93) * t93 + (-t35 * t148 + (t105 * t162 + t109 * t161) * t98 + t122 + t164 * t92) * t92 + ((t104 * t150 + t108 * t163 + t165 * t93) * t93 + (t104 * t162 + t108 * t161 + t165 * t92) * t92) * t101; m(4) * t21 + m(5) * t13 + m(6) * t5; m(6) * (t32 ^ 2 + t33 ^ 2 + t5 ^ 2) + m(5) * (t13 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(4) * (t137 * t75 ^ 2 + t21 ^ 2) + (t154 * t90 + t3) * t92 + (-t4 + t155 * t91 + (t154 * t93 + t155 * t92) * t92) * t93; 0.2e1 * ((t16 * t93 + t17 * t92) * t151 + (t22 * t93 + t23 * t92) * t152) * t98; t149 * t101; m(6) * (-t101 * t5 + (t32 * t92 + t33 * t93) * t98) + m(5) * (-t101 * t13 + (t36 * t92 + t37 * t93) * t98); 0.2e1 * (t152 + t151) * (t137 * t98 ^ 2 + t101 ^ 2); t20 + m(6) * (t18 * t16 + t19 * t17) + (t122 * t93 + t123 * t92) * t101; m(6) * t12; t2 * t162 + t1 * t150 + (t10 * t92 - t11 * t93) * t156 + m(6) * (t12 * t5 + t18 * t33 + t19 * t32) + (t4 * t150 + t3 * t161) * t101; m(6) * (-t12 * t101 + (t18 * t93 + t19 * t92) * t98); m(6) * (t12 ^ 2 + t18 ^ 2 + t19 ^ 2) + t98 * t20 + (t93 * t1 + t92 * t2 + t98 * (t10 * t93 + t11 * t92)) * t101;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t38(1), t38(2), t38(4), t38(7), t38(11); t38(2), t38(3), t38(5), t38(8), t38(12); t38(4), t38(5), t38(6), t38(9), t38(13); t38(7), t38(8), t38(9), t38(10), t38(14); t38(11), t38(12), t38(13), t38(14), t38(15);];
Mq = res;
