% Calculate joint inertia matrix for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:37
% DurationCPUTime: 0.72s
% Computational Cost: add. (2183->141), mult. (3012->239), div. (0->0), fcn. (3250->8), ass. (0->84)
t80 = sin(pkin(7));
t77 = t80 ^ 2;
t81 = cos(pkin(7));
t78 = t81 ^ 2;
t123 = t77 + t78;
t79 = qJ(2) + qJ(3);
t75 = sin(t79);
t76 = cos(t79);
t88 = Icges(4,5) * t76 - Icges(4,6) * t75;
t51 = -Icges(4,3) * t81 + t80 * t88;
t52 = Icges(4,3) * t80 + t81 * t88;
t116 = t75 * t80;
t104 = Icges(5,3) * t75;
t84 = cos(qJ(4));
t110 = t81 * t84;
t82 = sin(qJ(4));
t113 = t80 * t82;
t65 = -t113 * t76 - t110;
t111 = t81 * t82;
t112 = t80 * t84;
t66 = t112 * t76 - t111;
t30 = Icges(5,5) * t66 + Icges(5,6) * t65 + t104 * t80;
t105 = Icges(5,6) * t75;
t32 = Icges(5,4) * t66 + Icges(5,2) * t65 + t105 * t80;
t106 = Icges(5,5) * t75;
t34 = Icges(5,1) * t66 + Icges(5,4) * t65 + t106 * t80;
t14 = t116 * t30 + t32 * t65 + t34 * t66;
t67 = -t111 * t76 + t112;
t68 = t110 * t76 + t113;
t31 = Icges(5,5) * t68 + Icges(5,6) * t67 + t104 * t81;
t33 = Icges(5,4) * t68 + Icges(5,2) * t67 + t105 * t81;
t35 = Icges(5,1) * t68 + Icges(5,4) * t67 + t106 * t81;
t15 = t116 * t31 + t33 * t65 + t35 * t66;
t8 = -t14 * t81 + t15 * t80;
t90 = Icges(4,4) * t76 - Icges(4,2) * t75;
t92 = Icges(4,1) * t76 - Icges(4,4) * t75;
t96 = -(Icges(4,6) * t80 + t81 * t90) * t75 + (Icges(4,5) * t80 + t81 * t92) * t76;
t97 = (-Icges(4,6) * t81 + t80 * t90) * t75 - (-Icges(4,5) * t81 + t80 * t92) * t76;
t121 = -t78 * t51 - (t96 * t80 + (-t52 + t97) * t81) * t80 - t8;
t83 = sin(qJ(2));
t120 = pkin(2) * t83;
t115 = t75 * t81;
t16 = t115 * t30 + t32 * t67 + t34 * t68;
t17 = t115 * t31 + t33 * t67 + t35 * t68;
t9 = -t16 * t81 + t17 * t80;
t117 = (t77 * t52 + t9 + (t97 * t81 + (-t51 + t96) * t80) * t81) * t80;
t45 = -Icges(5,3) * t76 + (Icges(5,5) * t84 - Icges(5,6) * t82) * t75;
t114 = t76 * t45;
t85 = cos(qJ(2));
t109 = t123 * t85 * pkin(2);
t25 = t123 * (rSges(4,1) * t76 - rSges(4,2) * t75);
t48 = -t76 * rSges(5,3) + (rSges(5,1) * t84 - rSges(5,2) * t82) * t75;
t108 = -pkin(3) * t75 + pkin(6) * t76 - t48;
t70 = rSges(4,1) * t75 + rSges(4,2) * t76;
t103 = -t70 - t120;
t36 = rSges(5,1) * t66 + rSges(5,2) * t65 + rSges(5,3) * t116;
t37 = rSges(5,1) * t68 + rSges(5,2) * t67 + rSges(5,3) * t115;
t20 = t80 * t36 + t81 * t37 + t123 * (pkin(3) * t76 + pkin(6) * t75);
t18 = -t76 * t30 + (-t32 * t82 + t34 * t84) * t75;
t19 = -t76 * t31 + (-t33 * t82 + t35 * t84) * t75;
t46 = -Icges(5,6) * t76 + (Icges(5,4) * t84 - Icges(5,2) * t82) * t75;
t47 = -Icges(5,5) * t76 + (Icges(5,1) * t84 - Icges(5,4) * t82) * t75;
t3 = -(t46 * t65 + t47 * t66) * t76 + (t15 * t81 + (t14 - t114) * t80) * t75;
t4 = -(t46 * t67 + t47 * t68) * t76 + (t16 * t80 + (t17 - t114) * t81) * t75;
t102 = t80 * t4 / 0.2e1 - t76 * (-t18 * t81 + t19 * t80) / 0.2e1 - t81 * t3 / 0.2e1 + t8 * t116 / 0.2e1 + t9 * t115 / 0.2e1;
t101 = t108 - t120;
t89 = Icges(3,5) * t85 - Icges(3,6) * t83;
t87 = t121 * t81 + t117;
t73 = rSges(3,1) * t83 + rSges(3,2) * t85;
t60 = Icges(3,3) * t80 + t81 * t89;
t59 = -Icges(3,3) * t81 + t80 * t89;
t50 = t103 * t81;
t49 = t103 * t80;
t40 = t108 * t81;
t39 = t108 * t80;
t38 = t123 * (rSges(3,1) * t85 - rSges(3,2) * t83);
t27 = t101 * t81;
t26 = t101 * t80;
t24 = -t115 * t48 - t37 * t76;
t23 = t116 * t48 + t36 * t76;
t22 = t25 + t109;
t21 = (t36 * t81 - t37 * t80) * t75;
t11 = t20 + t109;
t1 = [m(2) + m(3) + m(4) + m(5); m(3) * t38 + m(4) * t22 + m(5) * t11; t80 * t77 * t60 + m(3) * (t123 * t73 ^ 2 + t38 ^ 2) + m(4) * (t22 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t11 ^ 2 + t26 ^ 2 + t27 ^ 2) + t117 + (-t78 * t59 + (-t80 * t59 + t81 * t60) * t80 + t121) * t81; m(4) * t25 + m(5) * t20; m(4) * (t25 * t22 + (-t49 * t80 - t50 * t81) * t70) + m(5) * (t11 * t20 + t26 * t39 + t27 * t40) + t87; m(4) * (t123 * t70 ^ 2 + t25 ^ 2) + m(5) * (t20 ^ 2 + t39 ^ 2 + t40 ^ 2) + t87; m(5) * t21; m(5) * (t11 * t21 + t23 * t27 + t24 * t26) + t102; m(5) * (t20 * t21 + t23 * t40 + t24 * t39) + t102; m(5) * (t21 ^ 2 + t23 ^ 2 + t24 ^ 2) + t4 * t115 + t3 * t116 - t76 * (t76 ^ 2 * t45 + (t19 * t81 + t18 * t80 - (-t46 * t82 + t47 * t84) * t76) * t75);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
