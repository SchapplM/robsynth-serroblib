% Calculate joint inertia matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:14
% EndTime: 2019-03-09 01:37:16
% DurationCPUTime: 1.10s
% Computational Cost: add. (2016->200), mult. (4363->314), div. (0->0), fcn. (5332->8), ass. (0->100)
t89 = sin(qJ(5));
t138 = Icges(6,5) * t89;
t92 = cos(qJ(5));
t129 = -t92 / 0.2e1;
t137 = t138 / 0.2e1 - Icges(6,6) * t129;
t110 = sin(pkin(9));
t87 = cos(pkin(9));
t90 = sin(qJ(1));
t93 = cos(qJ(1));
t64 = t93 * t110 + t90 * t87;
t132 = t64 ^ 2;
t63 = -t90 * t110 + t93 * t87;
t133 = t63 ^ 2;
t136 = -t133 - t132;
t135 = m(5) + m(6) + m(7);
t122 = t64 * t92;
t123 = t64 * t89;
t88 = sin(qJ(6));
t121 = t88 * t92;
t91 = cos(qJ(6));
t43 = -t64 * t121 - t63 * t91;
t120 = t91 * t92;
t44 = t64 * t120 - t63 * t88;
t28 = t44 * rSges(7,1) + t43 * rSges(7,2) + rSges(7,3) * t123;
t134 = pkin(5) * t122 + pkin(8) * t123 + t28;
t131 = -t63 / 0.2e1;
t75 = rSges(6,1) * t89 + rSges(6,2) * t92;
t128 = m(6) * t75;
t127 = pkin(5) * t92;
t51 = -Icges(7,6) * t92 + (Icges(7,4) * t91 - Icges(7,2) * t88) * t89;
t126 = t51 * t88;
t125 = t63 * t89;
t124 = t64 * rSges(6,3);
t119 = -pkin(1) - qJ(3);
t53 = -t92 * rSges(7,3) + (rSges(7,1) * t91 - rSges(7,2) * t88) * t89;
t118 = pkin(5) * t89 - pkin(8) * t92 + t53;
t116 = t93 * pkin(1) + t90 * qJ(2);
t113 = Icges(7,5) * t89;
t112 = Icges(7,6) * t89;
t111 = Icges(7,3) * t89;
t109 = t93 * qJ(3) + t116;
t22 = Icges(7,5) * t44 + Icges(7,6) * t43 + t64 * t111;
t24 = Icges(7,4) * t44 + Icges(7,2) * t43 + t64 * t112;
t26 = Icges(7,1) * t44 + Icges(7,4) * t43 + t64 * t113;
t10 = -t92 * t22 + (-t24 * t88 + t26 * t91) * t89;
t50 = -Icges(7,3) * t92 + (Icges(7,5) * t91 - Icges(7,6) * t88) * t89;
t52 = -Icges(7,5) * t92 + (Icges(7,1) * t91 - Icges(7,4) * t88) * t89;
t15 = t50 * t123 + t43 * t51 + t44 * t52;
t108 = t10 / 0.2e1 + t15 / 0.2e1;
t45 = t63 * t121 - t64 * t91;
t46 = -t63 * t120 - t64 * t88;
t23 = Icges(7,5) * t46 + Icges(7,6) * t45 - t63 * t111;
t25 = Icges(7,4) * t46 + Icges(7,2) * t45 - t63 * t112;
t27 = Icges(7,1) * t46 + Icges(7,4) * t45 - t63 * t113;
t11 = -t92 * t23 + (-t25 * t88 + t27 * t91) * t89;
t16 = -t50 * t125 + t45 * t51 + t46 * t52;
t107 = -t16 / 0.2e1 - t11 / 0.2e1;
t79 = t90 ^ 2 + t93 ^ 2;
t106 = (m(4) + t135) * t79;
t105 = t90 * pkin(3) + t109;
t104 = rSges(6,1) * t92 - rSges(6,2) * t89;
t103 = -t46 * rSges(7,1) - t45 * rSges(7,2);
t98 = Icges(6,5) * t92 - Icges(6,6) * t89;
t97 = rSges(6,1) * t122 - rSges(6,2) * t123 - t63 * rSges(6,3);
t83 = t93 * qJ(2);
t96 = t93 * pkin(3) + t119 * t90 + t83;
t95 = t64 * pkin(7) + t96;
t94 = t64 * pkin(4) - pkin(7) * t63 + t105;
t77 = rSges(2,1) * t93 - t90 * rSges(2,2);
t76 = -t90 * rSges(2,1) - rSges(2,2) * t93;
t55 = -rSges(3,2) * t93 + t90 * rSges(3,3) + t116;
t54 = t93 * rSges(3,3) + t83 + (rSges(3,2) - pkin(1)) * t90;
t49 = t90 * rSges(4,1) + rSges(4,3) * t93 + t109;
t48 = t93 * rSges(4,1) + t83 + (-rSges(4,3) + t119) * t90;
t47 = t89 * t91 * t52;
t40 = rSges(5,1) * t64 + rSges(5,2) * t63 + t105;
t39 = t63 * rSges(5,1) - t64 * rSges(5,2) + t96;
t34 = -Icges(6,3) * t64 - t98 * t63;
t32 = t118 * t63;
t31 = t118 * t64;
t30 = -t89 * t126 - t92 * t50 + t47;
t29 = -rSges(7,3) * t125 - t103;
t21 = t94 + t97;
t20 = t124 + (pkin(4) + t104) * t63 + t95;
t19 = t64 * t97 + (t104 * t63 + t124) * t63;
t18 = -t53 * t125 + t29 * t92;
t17 = -t53 * t123 - t28 * t92;
t14 = t94 + t134;
t13 = (t127 + pkin(4) + (rSges(7,3) + pkin(8)) * t89) * t63 + t95 + t103;
t12 = (t28 * t63 + t29 * t64) * t89;
t9 = t134 * t64 + (-t29 + (pkin(8) * t89 + t127) * t63) * t63;
t8 = -t23 * t125 + t25 * t45 + t27 * t46;
t7 = -t22 * t125 + t24 * t45 + t26 * t46;
t6 = t23 * t123 + t25 * t43 + t27 * t44;
t5 = t22 * t123 + t24 * t43 + t26 * t44;
t4 = -t63 * t7 - t64 * t8;
t3 = -t5 * t63 - t6 * t64;
t2 = -t16 * t92 + (-t63 * t8 + t64 * t7) * t89;
t1 = -t15 * t92 + (t5 * t64 - t6 * t63) * t89;
t33 = [Icges(3,1) + Icges(4,2) + Icges(2,3) + Icges(5,3) + t47 + (Icges(6,2) * t92 - t50) * t92 + m(7) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t48 ^ 2 + t49 ^ 2) + m(5) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t54 ^ 2 + t55 ^ 2) + m(2) * (t76 ^ 2 + t77 ^ 2) + (Icges(6,1) * t89 + 0.2e1 * Icges(6,4) * t92 - t126) * t89; m(7) * (t90 * t13 - t14 * t93) + m(6) * (t90 * t20 - t21 * t93) + m(4) * (t90 * t48 - t49 * t93) + m(5) * (t90 * t39 - t40 * t93) + m(3) * (t90 * t54 - t55 * t93); m(3) * t79 + t106; m(7) * (t13 * t93 + t90 * t14) + m(6) * (t20 * t93 + t90 * t21) + m(4) * (t48 * t93 + t90 * t49) + m(5) * (t39 * t93 + t90 * t40); 0; t106; 0; 0; 0; t135; m(7) * (-t13 * t31 + t14 * t32) + (t133 / 0.2e1 + t132 / 0.2e1) * (Icges(6,6) * t92 + t138) + (-t20 * t128 + t137 * t64 + t107) * t64 + (t21 * t128 + t137 * t63 - t108) * t63; m(7) * (-t31 * t90 - t32 * t93) + (-t63 * t93 - t64 * t90) * t128; m(7) * (-t31 * t93 + t32 * t90) + (t63 * t90 - t64 * t93) * t128; m(6) * t19 + m(7) * t9; m(6) * (-t136 * t75 ^ 2 + t19 ^ 2) + m(7) * (t31 ^ 2 + t32 ^ 2 + t9 ^ 2) + (-t132 * t34 - t4) * t64 + (-t63 * t34 * t64 - t3 + t136 * (-Icges(6,3) * t63 + t98 * t64)) * t63; -t30 * t92 + m(7) * (t13 * t18 + t14 * t17) + (t107 * t63 + t108 * t64) * t89; m(7) * (-t17 * t93 + t18 * t90); m(7) * (t17 * t90 + t18 * t93); m(7) * t12; t1 * t131 - t64 * t2 / 0.2e1 + (-t10 * t63 - t11 * t64) * t129 + m(7) * (t12 * t9 + t17 * t32 - t18 * t31) + (t64 * t3 / 0.2e1 + t4 * t131) * t89; t92 ^ 2 * t30 + m(7) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + (t64 * t1 - t63 * t2 - t92 * (t10 * t64 - t11 * t63)) * t89;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t33(1) t33(2) t33(4) t33(7) t33(11) t33(16); t33(2) t33(3) t33(5) t33(8) t33(12) t33(17); t33(4) t33(5) t33(6) t33(9) t33(13) t33(18); t33(7) t33(8) t33(9) t33(10) t33(14) t33(19); t33(11) t33(12) t33(13) t33(14) t33(15) t33(20); t33(16) t33(17) t33(18) t33(19) t33(20) t33(21);];
Mq  = res;
