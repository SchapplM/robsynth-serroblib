% Calculate joint inertia matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:22
% DurationCPUTime: 0.97s
% Computational Cost: add. (1825->184), mult. (3997->293), div. (0->0), fcn. (4966->8), ass. (0->95)
t82 = sin(qJ(4));
t128 = Icges(5,5) * t82;
t127 = -t128 / 0.2e1;
t101 = sin(pkin(8));
t102 = cos(pkin(8));
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t63 = t101 * t86 - t102 * t83;
t123 = t63 ^ 2;
t62 = -t101 * t83 - t102 * t86;
t124 = t62 ^ 2;
t126 = t123 + t124;
t85 = cos(qJ(4));
t116 = t62 * t85;
t117 = t62 * t82;
t81 = sin(qJ(5));
t114 = t81 * t85;
t84 = cos(qJ(5));
t45 = t114 * t62 + t63 * t84;
t112 = t84 * t85;
t46 = -t112 * t62 + t63 * t81;
t27 = t46 * rSges(6,1) + t45 * rSges(6,2) - rSges(6,3) * t117;
t125 = -pkin(4) * t116 - pkin(7) * t117 + t27;
t122 = -t62 / 0.2e1;
t121 = t85 / 0.2e1;
t69 = -rSges(5,1) * t82 - rSges(5,2) * t85;
t120 = m(5) * t69;
t119 = pkin(4) * t85;
t118 = t62 * rSges(5,3);
t115 = t63 * t82;
t51 = Icges(6,5) * t85 + (-Icges(6,1) * t84 + Icges(6,4) * t81) * t82;
t113 = t84 * t51;
t49 = Icges(6,3) * t85 + (-Icges(6,5) * t84 + Icges(6,6) * t81) * t82;
t50 = Icges(6,6) * t85 + (-Icges(6,4) * t84 + Icges(6,2) * t81) * t82;
t111 = t82 * t81 * t50 + t85 * t49;
t52 = t85 * rSges(6,3) + (-rSges(6,1) * t84 + rSges(6,2) * t81) * t82;
t110 = pkin(4) * t82 - pkin(7) * t85 - t52;
t108 = t86 * pkin(1) + t83 * qJ(2);
t106 = Icges(5,4) * t85;
t105 = Icges(6,5) * t82;
t104 = Icges(6,6) * t82;
t103 = Icges(6,3) * t82;
t100 = t86 * pkin(2) + t108;
t43 = t114 * t63 - t62 * t84;
t44 = -t112 * t63 - t62 * t81;
t20 = Icges(6,5) * t44 + Icges(6,6) * t43 - t103 * t63;
t22 = Icges(6,4) * t44 + Icges(6,2) * t43 - t104 * t63;
t24 = Icges(6,1) * t44 + Icges(6,4) * t43 - t105 * t63;
t10 = t85 * t20 + (t22 * t81 - t24 * t84) * t82;
t15 = -t115 * t49 + t43 * t50 + t44 * t51;
t99 = -t10 / 0.2e1 - t15 / 0.2e1;
t21 = Icges(6,5) * t46 + Icges(6,6) * t45 - t103 * t62;
t23 = Icges(6,4) * t46 + Icges(6,2) * t45 - t104 * t62;
t25 = Icges(6,1) * t46 + Icges(6,4) * t45 - t105 * t62;
t11 = t85 * t21 + (t23 * t81 - t25 * t84) * t82;
t16 = -t117 * t49 + t45 * t50 + t46 * t51;
t98 = -t11 / 0.2e1 - t16 / 0.2e1;
t78 = t86 * qJ(2);
t97 = t78 + (-pkin(1) - pkin(2)) * t83;
t96 = -rSges(5,1) * t85 + rSges(5,2) * t82;
t95 = -t44 * rSges(6,1) - t43 * rSges(6,2);
t91 = Icges(5,2) * t82 - t106;
t90 = -Icges(5,5) * t85 + Icges(5,6) * t82;
t89 = -rSges(5,1) * t116 + rSges(5,2) * t117 + t63 * rSges(5,3);
t88 = t62 * pkin(6) + t97;
t87 = -t62 * pkin(3) + pkin(6) * t63 + t100;
t71 = rSges(2,1) * t86 - rSges(2,2) * t83;
t70 = -rSges(2,1) * t83 - rSges(2,2) * t86;
t54 = rSges(3,1) * t86 + rSges(3,3) * t83 + t108;
t53 = t86 * rSges(3,3) + t78 + (-rSges(3,1) - pkin(1)) * t83;
t40 = -rSges(4,1) * t62 - rSges(4,2) * t63 + t100;
t39 = t63 * rSges(4,1) - t62 * rSges(4,2) + t97;
t34 = Icges(5,3) * t63 + t62 * t90;
t32 = t110 * t62;
t31 = t110 * t63;
t30 = (-t113 * t82 + t111) * t85;
t29 = t87 + t89;
t28 = t118 + (pkin(3) - t96) * t63 + t88;
t26 = -rSges(6,3) * t115 - t95;
t19 = t62 * t89 + (t63 * t96 - t118) * t63;
t18 = t117 * t52 + t27 * t85;
t17 = -t115 * t52 - t26 * t85;
t14 = t87 + t125;
t13 = (t119 + pkin(3) + (rSges(6,3) + pkin(7)) * t82) * t63 + t88 + t95;
t12 = (-t26 * t62 + t27 * t63) * t82;
t9 = t125 * t62 + (t26 + (-pkin(7) * t82 - t119) * t63) * t63;
t8 = -t117 * t21 + t23 * t45 + t25 * t46;
t7 = -t117 * t20 + t22 * t45 + t24 * t46;
t6 = -t115 * t21 + t23 * t43 + t25 * t44;
t5 = -t115 * t20 + t22 * t43 + t24 * t44;
t4 = -t62 * t7 + t63 * t8;
t3 = -t5 * t62 + t6 * t63;
t2 = t16 * t85 + (-t62 * t8 - t63 * t7) * t82;
t1 = t15 * t85 + (-t5 * t63 - t6 * t62) * t82;
t33 = [-t85 * (-Icges(5,4) * t82 - Icges(5,2) * t85) + Icges(3,2) + Icges(2,3) + Icges(4,3) + (Icges(5,1) * t82 + t106 - t113) * t82 + m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t28 ^ 2 + t29 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t53 ^ 2 + t54 ^ 2) + m(2) * (t70 ^ 2 + t71 ^ 2) + t111; m(6) * (t13 * t83 - t14 * t86) + m(5) * (t28 * t83 - t29 * t86) + m(4) * (t39 * t83 - t40 * t86) + m(3) * (t53 * t83 - t54 * t86); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t83 ^ 2 + t86 ^ 2); 0; 0; m(4) + m(5) + m(6); (-t85 * (Icges(5,6) * t63 + t91 * t62) / 0.2e1 + t63 * t127 - t98) * t63 + ((-Icges(5,6) * t62 + t63 * t91) * t121 + t62 * t127 + t99) * t62 + m(6) * (t13 * t32 + t14 * t31) + (-t28 * t62 - t29 * t63) * t120 + (t124 / 0.2e1 + t123 / 0.2e1) * (-Icges(5,6) * t85 - t128); m(6) * (-t31 * t86 + t32 * t83) + (-t62 * t83 + t63 * t86) * t120; -m(5) * t19 - m(6) * t9; m(5) * (t126 * t69 ^ 2 + t19 ^ 2) + m(6) * (t31 ^ 2 + t32 ^ 2 + t9 ^ 2) + (t123 * t34 + t4) * t63 + (t62 * t34 * t63 - t3 - t126 * (-Icges(5,3) * t62 + t63 * t90)) * t62; m(6) * (t13 * t17 + t14 * t18) + t30 + (t62 * t98 + t63 * t99) * t82; m(6) * (t17 * t83 - t18 * t86); -m(6) * t12; m(6) * (t12 * t9 + t17 * t32 + t18 * t31) + t63 * t2 / 0.2e1 + t1 * t122 + (-t10 * t62 + t11 * t63) * t121 + (t4 * t122 - t63 * t3 / 0.2e1) * t82; m(6) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + t85 * t30 + (-t62 * t2 - t63 * t1 + t85 * (-t10 * t63 - t11 * t62)) * t82;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t33(1), t33(2), t33(4), t33(7), t33(11); t33(2), t33(3), t33(5), t33(8), t33(12); t33(4), t33(5), t33(6), t33(9), t33(13); t33(7), t33(8), t33(9), t33(10), t33(14); t33(11), t33(12), t33(13), t33(14), t33(15);];
Mq = res;
