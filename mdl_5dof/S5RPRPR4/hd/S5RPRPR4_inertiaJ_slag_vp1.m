% Calculate joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:23
% EndTime: 2022-01-23 09:22:24
% DurationCPUTime: 0.82s
% Computational Cost: add. (1929->178), mult. (1340->250), div. (0->0), fcn. (1180->10), ass. (0->88)
t86 = qJ(3) + pkin(9);
t78 = sin(t86);
t80 = cos(t86);
t89 = sin(qJ(3));
t91 = cos(qJ(3));
t145 = Icges(4,5) * t91 + Icges(5,5) * t80 - Icges(4,6) * t89 - Icges(5,6) * t78;
t144 = Icges(4,3) + Icges(5,3);
t87 = qJ(1) + pkin(8);
t79 = sin(t87);
t143 = t79 * pkin(6);
t81 = cos(t87);
t142 = t79 * t81;
t82 = qJ(5) + t86;
t73 = sin(t82);
t74 = cos(t82);
t112 = rSges(6,1) * t74 - rSges(6,2) * t73;
t113 = rSges(5,1) * t80 - rSges(5,2) * t78;
t141 = t144 * t81 - t145 * t79;
t140 = t144 * t79 + t145 * t81;
t76 = t79 ^ 2;
t77 = t81 ^ 2;
t139 = t79 / 0.2e1;
t138 = -t81 / 0.2e1;
t137 = pkin(3) * t89;
t90 = sin(qJ(1));
t136 = t90 * pkin(1);
t75 = t91 * pkin(3) + pkin(2);
t135 = rSges(4,1) * t91;
t132 = rSges(4,2) * t89;
t129 = t81 * rSges(4,3);
t88 = -qJ(4) - pkin(6);
t59 = t81 * t75;
t72 = t81 * pkin(6);
t128 = t79 * (t72 + (-pkin(2) + t75) * t79) + t81 * (-t81 * pkin(2) - t143 + t59);
t94 = t79 * rSges(6,3) + t112 * t81;
t9 = t79 * (-t81 * rSges(6,3) + t112 * t79) + t81 * t94;
t56 = pkin(4) * t80 + t75;
t85 = pkin(7) - t88;
t127 = t81 * t56 + t85 * t79;
t126 = t79 * rSges(4,3) + t81 * t135;
t125 = t76 + t77;
t124 = Icges(4,4) * t89;
t123 = Icges(4,4) * t91;
t122 = Icges(5,4) * t78;
t121 = Icges(5,4) * t80;
t120 = Icges(6,4) * t73;
t119 = Icges(6,4) * t74;
t102 = Icges(6,1) * t74 - t120;
t48 = Icges(6,2) * t74 + t120;
t49 = Icges(6,1) * t73 + t119;
t105 = -t48 * t73 + t49 * t74;
t47 = Icges(6,5) * t73 + Icges(6,6) * t74;
t99 = -Icges(6,2) * t73 + t119;
t118 = (t105 * t81 + t74 * (Icges(6,6) * t79 + t81 * t99) + t73 * (Icges(6,5) * t79 + t102 * t81) + t79 * t47) * t139 + (t105 * t79 + t74 * (-Icges(6,6) * t81 + t79 * t99) + t73 * (-Icges(6,5) * t81 + t102 * t79) - t81 * t47) * t138;
t96 = Icges(6,5) * t74 - Icges(6,6) * t73;
t23 = -Icges(6,3) * t81 + t79 * t96;
t24 = Icges(6,3) * t79 + t81 * t96;
t117 = -t81 * (-t24 * t142 + t77 * t23) + t79 * (-t23 * t142 + t76 * t24);
t116 = -t78 * rSges(5,1) - t80 * rSges(5,2) - t137;
t114 = -t132 + t135;
t104 = Icges(4,1) * t91 - t124;
t103 = Icges(5,1) * t80 - t122;
t101 = -Icges(4,2) * t89 + t123;
t100 = -Icges(5,2) * t78 + t121;
t95 = t79 * rSges(5,3) + t113 * t81;
t50 = t73 * rSges(6,1) + t74 * rSges(6,2);
t93 = -pkin(4) * t78 - t137 - t50;
t92 = cos(qJ(1));
t84 = t92 * pkin(1);
t66 = t92 * rSges(2,1) - t90 * rSges(2,2);
t65 = -t90 * rSges(2,1) - t92 * rSges(2,2);
t64 = t89 * rSges(4,1) + t91 * rSges(4,2);
t44 = t81 * rSges(3,1) - t79 * rSges(3,2) + t84;
t43 = -t79 * rSges(3,1) - t81 * rSges(3,2) - t136;
t36 = t116 * t81;
t35 = t116 * t79;
t18 = t143 + t84 + (pkin(2) - t132) * t81 + t126;
t17 = t129 - t136 + t72 + (-pkin(2) - t114) * t79;
t16 = t93 * t81;
t15 = t93 * t79;
t14 = -t79 * t88 + t59 + t84 + t95;
t13 = -t136 + (rSges(5,3) - t88) * t81 + (-t113 - t75) * t79;
t12 = t81 * (-t81 * t132 + t126) + (t114 * t79 - t129) * t79;
t11 = t84 + t94 + t127;
t10 = -t136 + (rSges(6,3) + t85) * t81 + (-t112 - t56) * t79;
t4 = t81 * t95 + (-t81 * rSges(5,3) + t113 * t79) * t79 + t128;
t3 = t81 * (-t59 + t127) + (-t81 * t85 + (t56 - t75) * t79) * t79 + t9 + t128;
t1 = [t74 * t48 + t73 * t49 + t80 * (Icges(5,2) * t80 + t122) + t78 * (Icges(5,1) * t78 + t121) + t91 * (Icges(4,2) * t91 + t124) + t89 * (Icges(4,1) * t89 + t123) + Icges(2,3) + Icges(3,3) + m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + m(3) * (t43 ^ 2 + t44 ^ 2) + m(2) * (t65 ^ 2 + t66 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(6) * (t16 * t10 + t15 * t11) + m(5) * (t36 * t13 + t35 * t14) + m(4) * (-t17 * t81 - t18 * t79) * t64 + t118 + (t80 * (Icges(5,6) * t79 + t100 * t81) + t78 * (Icges(5,5) * t79 + t103 * t81) + t91 * (Icges(4,6) * t79 + t101 * t81) + t89 * (Icges(4,5) * t79 + t104 * t81)) * t139 + (t80 * (-Icges(5,6) * t81 + t100 * t79) + t78 * (-Icges(5,5) * t81 + t103 * t79) + t91 * (-Icges(4,6) * t81 + t101 * t79) + t89 * (-Icges(4,5) * t81 + t104 * t79)) * t138 + (Icges(4,5) * t89 + Icges(5,5) * t78 + Icges(4,6) * t91 + Icges(5,6) * t80) * (t76 / 0.2e1 + t77 / 0.2e1); m(4) * t12 + m(5) * t4 + m(6) * t3; m(6) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2 + t4 ^ 2) + m(4) * (t125 * t64 ^ 2 + t12 ^ 2) + t117 + t141 * t81 * t77 + (t140 * t76 + (t140 * t81 + t141 * t79) * t81) * t79; m(6) * (t79 * t10 - t81 * t11) + m(5) * (t79 * t13 - t81 * t14); 0; m(6) * (-t81 * t15 + t79 * t16) + m(5) * (-t81 * t35 + t79 * t36); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t125; m(6) * (-t10 * t81 - t11 * t79) * t50 + t118; m(6) * t9; m(6) * (t9 * t3 + (-t15 * t79 - t16 * t81) * t50) + t117; 0; m(6) * (t125 * t50 ^ 2 + t9 ^ 2) + t117;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
