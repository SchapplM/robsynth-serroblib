% Calculate joint inertia matrix for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:15
% EndTime: 2019-12-31 16:50:16
% DurationCPUTime: 0.69s
% Computational Cost: add. (1650->156), mult. (2001->250), div. (0->0), fcn. (2134->8), ass. (0->89)
t73 = sin(qJ(3));
t113 = Icges(4,5) * t73;
t112 = t113 / 0.2e1;
t71 = qJ(1) + pkin(7);
t69 = cos(t71);
t108 = t69 ^ 2;
t68 = sin(t71);
t109 = t68 ^ 2;
t111 = t108 + t109;
t76 = cos(qJ(3));
t100 = t69 * t76;
t101 = t69 * t73;
t75 = cos(qJ(4));
t72 = sin(qJ(4));
t98 = t72 * t76;
t44 = t68 * t75 - t69 * t98;
t97 = t75 * t76;
t45 = t68 * t72 + t69 * t97;
t28 = t45 * rSges(5,1) + t44 * rSges(5,2) + rSges(5,3) * t101;
t110 = pkin(3) * t100 + pkin(6) * t101 + t28;
t107 = t68 / 0.2e1;
t106 = -t76 / 0.2e1;
t105 = pkin(3) * t76;
t74 = sin(qJ(1));
t104 = t74 * pkin(1);
t103 = t68 * t73;
t102 = t69 * rSges(4,3);
t47 = -Icges(5,6) * t76 + (Icges(5,4) * t75 - Icges(5,2) * t72) * t73;
t99 = t72 * t47;
t49 = -t76 * rSges(5,3) + (rSges(5,1) * t75 - rSges(5,2) * t72) * t73;
t96 = -pkin(3) * t73 + pkin(6) * t76 - t49;
t93 = Icges(4,4) * t76;
t92 = Icges(5,5) * t73;
t91 = Icges(5,6) * t73;
t90 = Icges(5,3) * t73;
t77 = cos(qJ(1));
t70 = t77 * pkin(1);
t89 = t69 * pkin(2) + t68 * pkin(5) + t70;
t42 = -t68 * t98 - t69 * t75;
t43 = t68 * t97 - t69 * t72;
t46 = -Icges(5,3) * t76 + (Icges(5,5) * t75 - Icges(5,6) * t72) * t73;
t48 = -Icges(5,5) * t76 + (Icges(5,1) * t75 - Icges(5,4) * t72) * t73;
t13 = t103 * t46 + t42 * t47 + t43 * t48;
t21 = Icges(5,5) * t43 + Icges(5,6) * t42 + t68 * t90;
t23 = Icges(5,4) * t43 + Icges(5,2) * t42 + t68 * t91;
t25 = Icges(5,1) * t43 + Icges(5,4) * t42 + t68 * t92;
t9 = -t76 * t21 + (-t23 * t72 + t25 * t75) * t73;
t88 = t9 / 0.2e1 + t13 / 0.2e1;
t22 = Icges(5,5) * t45 + Icges(5,6) * t44 + t69 * t90;
t24 = Icges(5,4) * t45 + Icges(5,2) * t44 + t69 * t91;
t26 = Icges(5,1) * t45 + Icges(5,4) * t44 + t69 * t92;
t10 = -t76 * t22 + (-t24 * t72 + t26 * t75) * t73;
t14 = t101 * t46 + t44 * t47 + t45 * t48;
t87 = t10 / 0.2e1 + t14 / 0.2e1;
t86 = t69 * pkin(5) - t104;
t85 = rSges(4,1) * t76 - rSges(4,2) * t73;
t84 = -t43 * rSges(5,1) - t42 * rSges(5,2);
t80 = -Icges(4,2) * t73 + t93;
t79 = Icges(4,5) * t76 - Icges(4,6) * t73;
t78 = rSges(4,1) * t100 - rSges(4,2) * t101 + t68 * rSges(4,3);
t58 = rSges(2,1) * t77 - rSges(2,2) * t74;
t57 = -rSges(2,1) * t74 - rSges(2,2) * t77;
t56 = rSges(4,1) * t73 + rSges(4,2) * t76;
t51 = rSges(3,1) * t69 - rSges(3,2) * t68 + t70;
t50 = -rSges(3,1) * t68 - rSges(3,2) * t69 - t104;
t39 = t73 * t75 * t48;
t33 = -Icges(4,3) * t69 + t68 * t79;
t32 = t96 * t69;
t31 = t96 * t68;
t30 = t78 + t89;
t29 = t102 + (-pkin(2) - t85) * t68 + t86;
t27 = rSges(5,3) * t103 - t84;
t20 = t69 * t78 + (t68 * t85 - t102) * t68;
t19 = -t76 * t46 - t73 * t99 + t39;
t18 = -t101 * t49 - t28 * t76;
t17 = t103 * t49 + t27 * t76;
t16 = t89 + t110;
t15 = (-t105 - pkin(2) + (-rSges(5,3) - pkin(6)) * t73) * t68 + t84 + t86;
t12 = (t27 * t69 - t28 * t68) * t73;
t11 = t110 * t69 + (t27 + (pkin(6) * t73 + t105) * t68) * t68;
t8 = t101 * t22 + t24 * t44 + t26 * t45;
t7 = t101 * t21 + t23 * t44 + t25 * t45;
t6 = t103 * t22 + t24 * t42 + t26 * t43;
t5 = t103 * t21 + t23 * t42 + t25 * t43;
t4 = t68 * t8 - t69 * t7;
t3 = -t5 * t69 + t6 * t68;
t2 = -t14 * t76 + (t68 * t7 + t69 * t8) * t73;
t1 = -t13 * t76 + (t5 * t68 + t6 * t69) * t73;
t34 = [Icges(2,3) + Icges(3,3) + t39 + (Icges(4,4) * t73 + Icges(4,2) * t76 - t46) * t76 + (Icges(4,1) * t73 + t93 - t99) * t73 + m(2) * (t57 ^ 2 + t58 ^ 2) + m(3) * (t50 ^ 2 + t51 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2); 0; m(3) + m(4) + m(5); ((-Icges(4,6) * t69 + t68 * t80) * t106 + t69 * t112 - t88) * t69 + (t76 * (Icges(4,6) * t68 + t80 * t69) / 0.2e1 + t68 * t112 + t87) * t68 + m(5) * (t15 * t32 + t16 * t31) + m(4) * (-t29 * t69 - t30 * t68) * t56 + (t109 / 0.2e1 + t108 / 0.2e1) * (Icges(4,6) * t76 + t113); m(4) * t20 + m(5) * t11; m(4) * (t111 * t56 ^ 2 + t20 ^ 2) + m(5) * (t11 ^ 2 + t31 ^ 2 + t32 ^ 2) + (-t108 * t33 - t3) * t69 + (-t68 * t33 * t69 + t4 + t111 * (Icges(4,3) * t68 + t69 * t79)) * t68; -t19 * t76 + m(5) * (t15 * t17 + t16 * t18) + (t68 * t88 + t69 * t87) * t73; m(5) * t12; m(5) * (t11 * t12 + t17 * t32 + t18 * t31) + t2 * t107 - t69 * t1 / 0.2e1 + (t10 * t68 - t9 * t69) * t106 + (t69 * t4 / 0.2e1 + t3 * t107) * t73; m(5) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + t76 ^ 2 * t19 + (t69 * t2 + t68 * t1 - t76 * (t10 * t69 + t68 * t9)) * t73;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t34(1), t34(2), t34(4), t34(7); t34(2), t34(3), t34(5), t34(8); t34(4), t34(5), t34(6), t34(9); t34(7), t34(8), t34(9), t34(10);];
Mq = res;
