% Calculate joint inertia matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:55
% EndTime: 2019-12-31 18:18:58
% DurationCPUTime: 1.26s
% Computational Cost: add. (2903->213), mult. (2497->315), div. (0->0), fcn. (2578->10), ass. (0->115)
t95 = qJ(3) + pkin(9);
t92 = cos(t95);
t96 = qJ(1) + pkin(8);
t93 = cos(t96);
t142 = t92 * t93;
t164 = Icges(5,6) * t92;
t102 = cos(qJ(3));
t160 = Icges(4,6) * t102;
t90 = sin(t95);
t161 = Icges(5,5) * t90;
t99 = sin(qJ(3));
t162 = Icges(4,5) * t99;
t163 = t160 / 0.2e1 + t161 / 0.2e1 + t162 / 0.2e1;
t159 = Icges(4,3) + Icges(5,3);
t158 = Icges(4,5) * t102 + Icges(5,5) * t92 - Icges(4,6) * t99 - Icges(5,6) * t90;
t91 = sin(t96);
t154 = -t158 * t91 + t159 * t93;
t153 = t158 * t93 + t159 * t91;
t144 = t90 * t93;
t101 = cos(qJ(5));
t126 = t91 * t101;
t98 = sin(qJ(5));
t139 = t93 * t98;
t61 = -t92 * t139 + t126;
t125 = t93 * t101;
t143 = t91 * t98;
t62 = t92 * t125 + t143;
t31 = t62 * rSges(6,1) + t61 * rSges(6,2) + rSges(6,3) * t144;
t152 = pkin(4) * t142 + pkin(7) * t144 + t31;
t88 = t91 ^ 2;
t89 = t93 ^ 2;
t151 = t91 / 0.2e1;
t150 = -t92 / 0.2e1;
t149 = pkin(3) * t99;
t148 = pkin(4) * t92;
t147 = rSges(4,2) * t99;
t100 = sin(qJ(1));
t146 = t100 * pkin(1);
t145 = t90 * t91;
t141 = t93 * rSges(4,3);
t97 = -qJ(4) - pkin(6);
t140 = t93 * t97;
t50 = -Icges(6,6) * t92 + (Icges(6,4) * t101 - Icges(6,2) * t98) * t90;
t138 = t98 * t50;
t87 = t102 * pkin(3) + pkin(2);
t75 = t93 * t87;
t86 = t93 * pkin(6);
t137 = t91 * (t140 + t86 + (-pkin(2) + t87) * t91) + t93 * (-t93 * pkin(2) + t75 + (-pkin(6) - t97) * t91);
t133 = rSges(4,1) * t102;
t135 = t91 * rSges(4,3) + t93 * t133;
t134 = t88 + t89;
t130 = Icges(5,4) * t92;
t129 = Icges(6,5) * t90;
t128 = Icges(6,6) * t90;
t127 = Icges(6,3) * t90;
t59 = -t92 * t143 - t125;
t60 = t92 * t126 - t139;
t24 = Icges(6,5) * t60 + Icges(6,6) * t59 + t91 * t127;
t26 = Icges(6,4) * t60 + Icges(6,2) * t59 + t91 * t128;
t28 = Icges(6,1) * t60 + Icges(6,4) * t59 + t91 * t129;
t10 = -t92 * t24 + (t101 * t28 - t26 * t98) * t90;
t47 = -Icges(6,3) * t92 + (Icges(6,5) * t101 - Icges(6,6) * t98) * t90;
t53 = -Icges(6,5) * t92 + (Icges(6,1) * t101 - Icges(6,4) * t98) * t90;
t14 = t47 * t145 + t59 * t50 + t60 * t53;
t123 = t10 / 0.2e1 + t14 / 0.2e1;
t25 = Icges(6,5) * t62 + Icges(6,6) * t61 + t93 * t127;
t27 = Icges(6,4) * t62 + Icges(6,2) * t61 + t93 * t128;
t29 = Icges(6,1) * t62 + Icges(6,4) * t61 + t93 * t129;
t11 = -t92 * t25 + (t101 * t29 - t27 * t98) * t90;
t15 = t47 * t144 + t61 * t50 + t62 * t53;
t122 = t11 / 0.2e1 + t15 / 0.2e1;
t121 = -t90 * rSges(5,1) - t92 * rSges(5,2) - t149;
t56 = -t92 * rSges(6,3) + (rSges(6,1) * t101 - rSges(6,2) * t98) * t90;
t120 = -t90 * pkin(4) + t92 * pkin(7) - t149 - t56;
t103 = cos(qJ(1));
t94 = t103 * pkin(1);
t118 = -t91 * t97 + t75 + t94;
t117 = rSges(5,1) * t92 - rSges(5,2) * t90;
t116 = -t60 * rSges(6,1) - t59 * rSges(6,2);
t113 = t133 - t147;
t109 = -Icges(5,2) * t90 + t130;
t107 = rSges(5,1) * t142 - rSges(5,2) * t144 + t91 * rSges(5,3);
t82 = t103 * rSges(2,1) - t100 * rSges(2,2);
t81 = -t100 * rSges(2,1) - t103 * rSges(2,2);
t80 = t99 * rSges(4,1) + t102 * rSges(4,2);
t64 = t93 * rSges(3,1) - t91 * rSges(3,2) + t94;
t63 = -t91 * rSges(3,1) - t93 * rSges(3,2) - t146;
t46 = t121 * t93;
t45 = t121 * t91;
t38 = t90 * t101 * t53;
t35 = t91 * pkin(6) + t94 + (pkin(2) - t147) * t93 + t135;
t34 = t141 - t146 + t86 + (-pkin(2) - t113) * t91;
t33 = t107 + t118;
t32 = -t146 + (rSges(5,3) - t97) * t93 + (-t117 - t87) * t91;
t30 = rSges(6,3) * t145 - t116;
t23 = t120 * t93;
t22 = t120 * t91;
t21 = t93 * (-t93 * t147 + t135) + (t113 * t91 - t141) * t91;
t20 = -t90 * t138 - t92 * t47 + t38;
t19 = -t56 * t144 - t92 * t31;
t18 = t56 * t145 + t92 * t30;
t17 = t118 + t152;
t16 = -t146 - t140 + (-t148 - t87 + (-rSges(6,3) - pkin(7)) * t90) * t91 + t116;
t13 = (t30 * t93 - t31 * t91) * t90;
t12 = t93 * t107 + (-t93 * rSges(5,3) + t117 * t91) * t91 + t137;
t9 = t25 * t144 + t61 * t27 + t62 * t29;
t8 = t24 * t144 + t61 * t26 + t62 * t28;
t7 = t25 * t145 + t59 * t27 + t60 * t29;
t6 = t24 * t145 + t59 * t26 + t60 * t28;
t5 = t152 * t93 + (t30 + (pkin(7) * t90 + t148) * t91) * t91 + t137;
t4 = -t8 * t93 + t9 * t91;
t3 = -t6 * t93 + t7 * t91;
t2 = -t15 * t92 + (t8 * t91 + t9 * t93) * t90;
t1 = -t14 * t92 + (t6 * t91 + t7 * t93) * t90;
t36 = [Icges(4,1) * t99 ^ 2 + Icges(2,3) + Icges(3,3) + t38 + (Icges(5,4) * t90 + Icges(5,2) * t92 - t47) * t92 + (Icges(5,1) * t90 + t130 - t138) * t90 + m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(3) * (t63 ^ 2 + t64 ^ 2) + m(2) * (t81 ^ 2 + t82 ^ 2) + (0.2e1 * Icges(4,4) * t99 + Icges(4,2) * t102) * t102; 0; m(3) + m(4) + m(5) + m(6); (t109 * t91 * t150 - t123 + (-Icges(5,6) * t150 + t163) * t93) * t93 + (t109 * t142 / 0.2e1 + t122 + (t164 / 0.2e1 + t163) * t91) * t91 + m(6) * (t23 * t16 + t22 * t17) + m(5) * (t46 * t32 + t45 * t33) + m(4) * (-t34 * t93 - t35 * t91) * t80 + (t160 + t161 + t162 + t164) * (t89 / 0.2e1 + t88 / 0.2e1); m(4) * t21 + m(5) * t12 + m(6) * t5; m(6) * (t22 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(4) * (t134 * t80 ^ 2 + t21 ^ 2) + (t153 * t88 + t4) * t91 + (-t3 + t154 * t89 + (t153 * t93 + t154 * t91) * t91) * t93; m(6) * (t91 * t16 - t93 * t17) + m(5) * (t91 * t32 - t93 * t33); 0; m(6) * (-t93 * t22 + t91 * t23) + m(5) * (-t93 * t45 + t91 * t46); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t134; -t20 * t92 + m(6) * (t18 * t16 + t19 * t17) + (t122 * t93 + t123 * t91) * t90; m(6) * t13; -t93 * t1 / 0.2e1 + m(6) * (t13 * t5 + t18 * t23 + t19 * t22) + t2 * t151 + (-t10 * t93 + t11 * t91) * t150 + (t3 * t151 + t93 * t4 / 0.2e1) * t90; m(6) * (t18 * t91 - t19 * t93); m(6) * (t13 ^ 2 + t18 ^ 2 + t19 ^ 2) + t92 ^ 2 * t20 + (t93 * t2 + t91 * t1 - t92 * (t10 * t91 + t11 * t93)) * t90;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t36(1), t36(2), t36(4), t36(7), t36(11); t36(2), t36(3), t36(5), t36(8), t36(12); t36(4), t36(5), t36(6), t36(9), t36(13); t36(7), t36(8), t36(9), t36(10), t36(14); t36(11), t36(12), t36(13), t36(14), t36(15);];
Mq = res;
