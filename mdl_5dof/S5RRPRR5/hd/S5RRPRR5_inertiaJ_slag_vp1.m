% Calculate joint inertia matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:31
% EndTime: 2020-01-03 12:03:33
% DurationCPUTime: 0.75s
% Computational Cost: add. (2367->167), mult. (1556->249), div. (0->0), fcn. (1342->10), ass. (0->93)
t102 = qJ(1) + qJ(2);
t95 = sin(t102);
t96 = cos(t102);
t146 = t95 * t96;
t150 = rSges(4,3) + qJ(3);
t101 = pkin(9) + qJ(4);
t94 = qJ(5) + t101;
t85 = sin(t94);
t86 = cos(t94);
t149 = -rSges(6,1) * t86 + rSges(6,2) * t85;
t92 = sin(t101);
t93 = cos(t101);
t148 = -rSges(5,1) * t93 + rSges(5,2) * t92;
t103 = sin(pkin(9));
t104 = cos(pkin(9));
t147 = rSges(4,1) * t104 - rSges(4,2) * t103 + pkin(2);
t35 = -t95 * rSges(6,3) + t149 * t96;
t87 = t104 * pkin(3) + pkin(2);
t64 = pkin(4) * t93 + t87;
t145 = t96 * t64 - t35;
t90 = t95 ^ 2;
t91 = t96 ^ 2;
t144 = -t95 / 0.2e1;
t143 = -t96 / 0.2e1;
t61 = t92 * rSges(5,1) + t93 * rSges(5,2);
t142 = m(5) * t61;
t54 = t85 * rSges(6,1) + t86 * rSges(6,2);
t141 = m(6) * t54;
t105 = -pkin(7) - qJ(3);
t100 = -pkin(8) + t105;
t136 = t96 * t100 + t95 * t64;
t135 = -t96 * t105 - t95 * t87;
t62 = t95 * rSges(3,1) + t96 * rSges(3,2);
t134 = t91 + t90;
t131 = Icges(5,4) * t92;
t130 = Icges(5,4) * t93;
t129 = Icges(6,4) * t85;
t128 = Icges(6,4) * t86;
t115 = -Icges(6,2) * t85 + t128;
t117 = Icges(6,1) * t86 - t129;
t52 = Icges(6,2) * t86 + t129;
t53 = Icges(6,1) * t85 + t128;
t120 = t52 * t85 - t53 * t86;
t51 = Icges(6,5) * t85 + Icges(6,6) * t86;
t127 = (t120 * t96 + t86 * (-Icges(6,6) * t95 - t115 * t96) + t85 * (-Icges(6,5) * t95 - t117 * t96) - t95 * t51) * t144 + (-t120 * t95 + t86 * (-Icges(6,6) * t96 + t115 * t95) + t85 * (-Icges(6,5) * t96 + t117 * t95) - t96 * t51) * t143;
t126 = pkin(4) * t92 + t54;
t63 = t96 * rSges(3,1) - t95 * rSges(3,2);
t113 = Icges(6,5) * t86 - Icges(6,6) * t85;
t29 = -Icges(6,3) * t96 + t113 * t95;
t30 = -Icges(6,3) * t95 - t113 * t96;
t125 = -t96 * (t30 * t146 + t91 * t29) - t95 * (t29 * t146 + t90 * t30);
t57 = Icges(5,2) * t93 + t131;
t58 = Icges(5,1) * t92 + t130;
t119 = t57 * t92 - t58 * t93;
t118 = Icges(5,1) * t93 - t131;
t116 = -Icges(5,2) * t92 + t130;
t114 = Icges(5,5) * t93 - Icges(5,6) * t92;
t112 = -t95 * rSges(5,3) + t148 * t96;
t56 = Icges(5,5) * t92 + Icges(5,6) * t93;
t111 = t127 + (t119 * t96 + t93 * (-Icges(5,6) * t95 - t116 * t96) + t92 * (-Icges(5,5) * t95 - t118 * t96) - t95 * t56) * t144 + (-t119 * t95 + t93 * (-Icges(5,6) * t96 + t116 * t95) + t92 * (-Icges(5,5) * t96 + t118 * t95) - t96 * t56) * t143;
t110 = Icges(4,2) * t104 ^ 2 + t86 * t52 + t85 * t53 + t93 * t57 + t92 * t58 + Icges(3,3) + (Icges(4,1) * t103 + 0.2e1 * Icges(4,4) * t104) * t103;
t109 = -t96 * rSges(5,3) - t148 * t95;
t108 = -t96 * rSges(6,3) - t149 * t95;
t25 = t147 * t96 + t150 * t95;
t20 = t109 - t135;
t16 = t108 + t136;
t68 = t96 * t87;
t21 = -t95 * t105 - t112 + t68;
t17 = -t95 * t100 + t145;
t24 = t147 * t95 - t150 * t96;
t107 = cos(qJ(1));
t106 = sin(qJ(1));
t99 = t107 * pkin(1);
t98 = t106 * pkin(1);
t73 = t107 * rSges(2,1) - t106 * rSges(2,2);
t72 = t106 * rSges(2,1) + t107 * rSges(2,2);
t49 = t63 + t99;
t48 = t98 + t62;
t37 = -Icges(5,3) * t95 - t114 * t96;
t36 = -Icges(5,3) * t96 + t114 * t95;
t28 = t126 * t96;
t27 = t126 * t95;
t26 = t95 * t108;
t23 = t25 + t99;
t22 = t98 + t24;
t19 = t21 + t99;
t18 = t20 + t98;
t15 = t17 + t99;
t14 = t16 + t98;
t13 = t95 * t109 - t96 * t112;
t12 = -t96 * t35 + t26;
t3 = t26 + (-t68 + t145) * t96 + (-t96 * (t100 - t105) + t135 + t136) * t95;
t1 = [Icges(2,3) + m(2) * (t72 ^ 2 + t73 ^ 2) + m(3) * (t48 ^ 2 + t49 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t14 ^ 2 + t15 ^ 2) + t110; m(3) * (t62 * t48 + t63 * t49) + m(4) * (t24 * t22 + t25 * t23) + m(5) * (t20 * t18 + t21 * t19) + m(6) * (t16 * t14 + t17 * t15) + t110; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t62 ^ 2 + t63 ^ 2) + t110; m(4) * (-t95 * t22 - t96 * t23) + m(5) * (-t95 * t18 - t96 * t19) + m(6) * (-t95 * t14 - t96 * t15); m(6) * (-t95 * t16 - t96 * t17) + m(5) * (-t95 * t20 - t96 * t21) + m(4) * (-t95 * t24 - t96 * t25); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t134; m(6) * (t28 * t14 - t27 * t15) + (t18 * t96 - t19 * t95) * t142 + t111; m(6) * (t28 * t16 - t27 * t17) + (t20 * t96 - t21 * t95) * t142 + t111; m(6) * (t27 * t96 - t28 * t95); m(5) * (t134 * t61 ^ 2 + t13 ^ 2) - t96 * (t37 * t146 + t91 * t36) - t95 * (t36 * t146 + t90 * t37) + m(6) * (t27 ^ 2 + t28 ^ 2 + t3 ^ 2) + t125; (t14 * t96 - t15 * t95) * t141 + t127; (t16 * t96 - t17 * t95) * t141 + t127; 0; m(6) * (t12 * t3 + (t27 * t95 + t28 * t96) * t54) + t125; m(6) * (t134 * t54 ^ 2 + t12 ^ 2) + t125;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
