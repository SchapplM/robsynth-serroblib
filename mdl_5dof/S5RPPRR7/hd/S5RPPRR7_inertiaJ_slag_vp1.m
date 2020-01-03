% Calculate joint inertia matrix for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:33
% DurationCPUTime: 0.84s
% Computational Cost: add. (1842->172), mult. (2205->272), div. (0->0), fcn. (2302->8), ass. (0->93)
t82 = cos(qJ(4));
t123 = Icges(5,5) * t82;
t79 = sin(qJ(4));
t122 = Icges(5,6) * t79;
t121 = t123 / 0.2e1 - t122 / 0.2e1;
t77 = qJ(1) + pkin(8);
t74 = sin(t77);
t72 = t74 ^ 2;
t75 = cos(t77);
t73 = t75 ^ 2;
t104 = t72 + t73;
t120 = (rSges(5,1) * t79 + rSges(5,2) * t82) * t75;
t119 = t75 / 0.2e1;
t116 = -pkin(2) - pkin(6);
t115 = pkin(4) * t79;
t80 = sin(qJ(1));
t114 = t80 * pkin(1);
t113 = t74 * t79;
t112 = t74 * t82;
t111 = t75 * t82;
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t50 = Icges(6,6) * t79 + (Icges(6,4) * t81 - Icges(6,2) * t78) * t82;
t110 = t78 * t50;
t109 = t78 * t79;
t108 = t79 * t81;
t49 = Icges(6,3) * t79 + (Icges(6,5) * t81 - Icges(6,6) * t78) * t82;
t51 = Icges(6,5) * t79 + (Icges(6,1) * t81 - Icges(6,4) * t78) * t82;
t107 = t82 * t81 * t51 + t79 * t49;
t45 = -t74 * t109 + t75 * t81;
t46 = t74 * t108 + t75 * t78;
t106 = t46 * rSges(6,1) + t45 * rSges(6,2);
t52 = t79 * rSges(6,3) + (rSges(6,1) * t81 - rSges(6,2) * t78) * t82;
t105 = t82 * pkin(4) + t79 * pkin(7) + t52;
t101 = Icges(6,5) * t82;
t100 = Icges(6,6) * t82;
t99 = Icges(6,3) * t82;
t98 = rSges(5,1) * t113 + rSges(5,2) * t112 + t75 * rSges(5,3);
t83 = cos(qJ(1));
t76 = t83 * pkin(1);
t97 = t75 * pkin(2) + t74 * qJ(3) + t76;
t13 = -t49 * t112 + t45 * t50 + t46 * t51;
t21 = Icges(6,5) * t46 + Icges(6,6) * t45 - t74 * t99;
t23 = Icges(6,4) * t46 + Icges(6,2) * t45 - t74 * t100;
t25 = Icges(6,1) * t46 + Icges(6,4) * t45 - t74 * t101;
t9 = t79 * t21 + (-t23 * t78 + t25 * t81) * t82;
t96 = -t9 / 0.2e1 - t13 / 0.2e1;
t47 = t75 * t109 + t74 * t81;
t48 = -t75 * t108 + t74 * t78;
t22 = Icges(6,5) * t48 + Icges(6,6) * t47 + t75 * t99;
t24 = Icges(6,4) * t48 + Icges(6,2) * t47 + t75 * t100;
t26 = Icges(6,1) * t48 + Icges(6,4) * t47 + t75 * t101;
t10 = t79 * t22 + (-t24 * t78 + t26 * t81) * t82;
t14 = t49 * t111 + t47 * t50 + t48 * t51;
t95 = t14 / 0.2e1 + t10 / 0.2e1;
t94 = (-rSges(6,3) - pkin(7)) * t82;
t93 = t75 * qJ(3) - t114;
t92 = t75 * pkin(6) + t97;
t90 = -t48 * rSges(6,1) - t47 * rSges(6,2);
t85 = Icges(5,5) * t79 + Icges(5,6) * t82;
t29 = t120 + (-rSges(5,3) + t116) * t74 + t93;
t30 = t92 + t98;
t84 = m(5) * (t74 * t29 - t75 * t30);
t66 = pkin(4) * t113;
t62 = t83 * rSges(2,1) - t80 * rSges(2,2);
t61 = t82 * rSges(5,1) - t79 * rSges(5,2);
t60 = -t80 * rSges(2,1) - t83 * rSges(2,2);
t54 = t75 * rSges(3,1) - t74 * rSges(3,2) + t76;
t53 = -t74 * rSges(3,1) - t75 * rSges(3,2) - t114;
t35 = Icges(5,3) * t75 + t85 * t74;
t34 = -t75 * rSges(4,2) + t74 * rSges(4,3) + t97;
t33 = t75 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t74 + t93;
t32 = t105 * t75;
t31 = t105 * t74;
t28 = rSges(6,3) * t111 - t90;
t27 = -rSges(6,3) * t112 + t106;
t20 = -t74 * t98 + (t74 * rSges(5,3) - t120) * t75;
t19 = (-t82 * t110 + t107) * t79;
t18 = t52 * t111 - t79 * t28;
t17 = t52 * t112 + t79 * t27;
t16 = t74 * t94 + t106 + t66 + t92;
t15 = t116 * t74 + (t94 + t115) * t75 + t90 + t93;
t12 = (-t27 * t75 - t28 * t74) * t82;
t11 = (pkin(7) * t112 - t27 - t66) * t74 + (t28 + (pkin(7) * t82 - t115) * t75) * t75;
t8 = t22 * t111 + t47 * t24 + t48 * t26;
t7 = t21 * t111 + t47 * t23 + t48 * t25;
t6 = -t22 * t112 + t45 * t24 + t46 * t26;
t5 = -t21 * t112 + t45 * t23 + t46 * t25;
t4 = t7 * t75 + t8 * t74;
t3 = t5 * t75 + t6 * t74;
t2 = t14 * t79 + (-t7 * t74 + t75 * t8) * t82;
t1 = t13 * t79 + (-t5 * t74 + t6 * t75) * t82;
t36 = [Icges(4,1) + Icges(2,3) + Icges(3,3) + (Icges(5,1) * t82 - t110) * t82 + m(6) * (t15 ^ 2 + t16 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(3) * (t53 ^ 2 + t54 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(2) * (t60 ^ 2 + t62 ^ 2) + t107 + (-0.2e1 * Icges(5,4) * t82 + Icges(5,2) * t79) * t79; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t74 * t15 - t75 * t16) + t84 + m(4) * (t74 * t33 - t75 * t34); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t104; m(6) * (t31 * t15 - t32 * t16) + t61 * t84 + (t72 / 0.2e1 + t73 / 0.2e1) * (-t122 + t123) + (t121 * t75 - t96) * t75 + (t121 * t74 + t95) * t74; m(5) * t20 + m(6) * t11; m(6) * (t31 * t74 + t32 * t75) + m(5) * t104 * t61; m(5) * (t104 * t61 ^ 2 + t20 ^ 2) + m(6) * (t11 ^ 2 + t31 ^ 2 + t32 ^ 2) + (t73 * t35 + t3) * t75 + (t74 * t35 * t75 + t4 + t104 * (Icges(5,3) * t74 - t85 * t75)) * t74; t19 + m(6) * (t18 * t15 + t17 * t16) + (t96 * t74 + t95 * t75) * t82; m(6) * t12; m(6) * (-t17 * t75 + t18 * t74); m(6) * (t12 * t11 - t17 * t32 + t18 * t31) + t1 * t119 + t74 * t2 / 0.2e1 + t79 * (t10 * t74 + t9 * t75) / 0.2e1 + (-t74 * t3 / 0.2e1 + t4 * t119) * t82; m(6) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + t79 * t19 + (-t74 * t1 + t75 * t2 + t79 * (t10 * t75 - t74 * t9)) * t82;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t36(1), t36(2), t36(4), t36(7), t36(11); t36(2), t36(3), t36(5), t36(8), t36(12); t36(4), t36(5), t36(6), t36(9), t36(13); t36(7), t36(8), t36(9), t36(10), t36(14); t36(11), t36(12), t36(13), t36(14), t36(15);];
Mq = res;
