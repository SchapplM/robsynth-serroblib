% Calculate joint inertia matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:33
% EndTime: 2019-12-31 18:05:35
% DurationCPUTime: 0.86s
% Computational Cost: add. (977->183), mult. (2311->287), div. (0->0), fcn. (2382->6), ass. (0->95)
t80 = cos(qJ(4));
t119 = Icges(5,5) * t80;
t118 = t119 / 0.2e1;
t78 = sin(qJ(1));
t74 = t78 ^ 2;
t81 = cos(qJ(1));
t75 = t81 ^ 2;
t66 = t74 + t75;
t77 = sin(qJ(4));
t117 = t77 / 0.2e1;
t116 = -t78 / 0.2e1;
t56 = m(5) * t66;
t115 = pkin(4) * t77;
t114 = -rSges(5,3) - pkin(6);
t113 = rSges(6,3) + pkin(7);
t76 = sin(qJ(5));
t79 = cos(qJ(5));
t40 = Icges(6,6) * t77 + (Icges(6,4) * t79 - Icges(6,2) * t76) * t80;
t112 = t76 * t40;
t111 = t77 * t81;
t110 = t78 * t76;
t109 = t78 * t79;
t108 = t78 * t80;
t107 = t80 * t81;
t106 = t81 * t76;
t105 = t81 * t79;
t104 = -pkin(1) - qJ(3);
t37 = Icges(6,3) * t77 + (Icges(6,5) * t79 - Icges(6,6) * t76) * t80;
t43 = Icges(6,5) * t77 + (Icges(6,1) * t79 - Icges(6,4) * t76) * t80;
t103 = t80 * t79 * t43 + t77 * t37;
t46 = t77 * rSges(6,3) + (rSges(6,1) * t79 - rSges(6,2) * t76) * t80;
t102 = t80 * pkin(4) + t77 * pkin(7) + t46;
t53 = -t77 * t106 - t109;
t54 = t77 * t105 - t110;
t101 = t54 * rSges(6,1) + t53 * rSges(6,2);
t100 = rSges(5,1) * t111 + rSges(5,2) * t107;
t99 = t81 * pkin(1) + t78 * qJ(2);
t98 = Icges(5,4) * t77;
t96 = Icges(6,5) * t80;
t95 = Icges(6,6) * t80;
t94 = Icges(6,3) * t80;
t93 = t56 + (m(4) + m(6)) * t66;
t92 = t81 * qJ(3) + t99;
t51 = -t77 * t110 + t105;
t52 = t77 * t109 + t106;
t12 = -t37 * t108 + t51 * t40 + t52 * t43;
t21 = Icges(6,5) * t52 + Icges(6,6) * t51 - t78 * t94;
t23 = Icges(6,4) * t52 + Icges(6,2) * t51 - t78 * t95;
t25 = Icges(6,1) * t52 + Icges(6,4) * t51 - t78 * t96;
t9 = t77 * t21 + (-t23 * t76 + t25 * t79) * t80;
t91 = -t9 / 0.2e1 - t12 / 0.2e1;
t22 = Icges(6,5) * t54 + Icges(6,6) * t53 - t81 * t94;
t24 = Icges(6,4) * t54 + Icges(6,2) * t53 - t81 * t95;
t26 = Icges(6,1) * t54 + Icges(6,4) * t53 - t81 * t96;
t10 = t77 * t22 + (-t24 * t76 + t26 * t79) * t80;
t13 = -t37 * t107 + t53 * t40 + t54 * t43;
t90 = -t13 / 0.2e1 - t10 / 0.2e1;
t89 = -rSges(5,1) * t77 - rSges(5,2) * t80;
t88 = -t52 * rSges(6,1) - t51 * rSges(6,2);
t84 = Icges(5,2) * t80 + t98;
t83 = Icges(5,5) * t77 + Icges(5,6) * t80;
t72 = t81 * qJ(2);
t29 = t72 + t114 * t81 + (t89 + t104) * t78;
t30 = t114 * t78 + t100 + t92;
t82 = m(5) * (t81 * t29 + t78 * t30);
t69 = pkin(4) * t111;
t64 = t81 * rSges(2,1) - t78 * rSges(2,2);
t63 = t80 * rSges(5,1) - t77 * rSges(5,2);
t62 = -t78 * rSges(2,1) - t81 * rSges(2,2);
t48 = -t81 * rSges(3,2) + t78 * rSges(3,3) + t99;
t47 = t81 * rSges(3,3) + t72 + (rSges(3,2) - pkin(1)) * t78;
t38 = Icges(5,3) * t81 + t83 * t78;
t36 = t78 * rSges(4,2) + t81 * rSges(4,3) + t92;
t35 = t81 * rSges(4,2) + t72 + (-rSges(4,3) + t104) * t78;
t32 = t102 * t81;
t31 = t102 * t78;
t28 = -rSges(6,3) * t107 + t101;
t27 = -rSges(6,3) * t108 - t88;
t20 = -t81 * t100 + t89 * t74;
t19 = t46 * t107 + t77 * t28;
t18 = -t46 * t108 - t77 * t27;
t17 = -t78 * pkin(6) - t113 * t107 + t101 + t69 + t92;
t16 = -t81 * pkin(6) + t72 + (t113 * t80 + t104 - t115) * t78 + t88;
t15 = (-t80 * t112 + t103) * t77;
t14 = (-t27 * t81 + t28 * t78) * t80;
t11 = (pkin(7) * t107 - t28 - t69) * t81 + (-t27 + (pkin(7) * t80 - t115) * t78) * t78;
t8 = -t22 * t107 + t53 * t24 + t54 * t26;
t7 = -t21 * t107 + t53 * t23 + t54 * t25;
t6 = -t22 * t108 + t51 * t24 + t52 * t26;
t5 = -t21 * t108 + t51 * t23 + t52 * t25;
t4 = t7 * t81 - t8 * t78;
t3 = t5 * t81 - t6 * t78;
t2 = t13 * t77 + (-t7 * t78 - t8 * t81) * t80;
t1 = t12 * t77 + (-t5 * t78 - t6 * t81) * t80;
t33 = [-t77 * (Icges(5,4) * t80 - Icges(5,2) * t77) + Icges(3,1) + Icges(4,1) + Icges(2,3) + (Icges(5,1) * t80 - t112 - t98) * t80 + m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(3) * (t47 ^ 2 + t48 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(2) * (t62 ^ 2 + t64 ^ 2) + t103; m(6) * (t78 * t16 - t81 * t17) + m(5) * (t78 * t29 - t81 * t30) + m(3) * (t78 * t47 - t81 * t48) + m(4) * (t78 * t35 - t81 * t36); m(3) * t66 + t93; m(6) * (t81 * t16 + t78 * t17) + t82 + m(4) * (t81 * t35 + t78 * t36); 0; t93; (-t77 * (Icges(5,6) * t81 + t84 * t78) / 0.2e1 + t81 * t118 - t91) * t81 + ((-Icges(5,6) * t78 + t84 * t81) * t117 + t78 * t118 + t90) * t78 + m(6) * (t32 * t16 + t31 * t17) + t63 * t82 + (t74 / 0.2e1 + t75 / 0.2e1) * (-Icges(5,6) * t77 + t119); m(6) * (-t31 * t81 + t32 * t78); m(6) * (t31 * t78 + t32 * t81) + t63 * t56; m(5) * (t66 * t63 ^ 2 + t20 ^ 2) + m(6) * (t11 ^ 2 + t31 ^ 2 + t32 ^ 2) + (t75 * t38 + t3) * t81 + (t78 * t38 * t81 - t4 - t66 * (-Icges(5,3) * t78 + t83 * t81)) * t78; m(6) * (t18 * t16 + t19 * t17) + t15 + (t91 * t78 + t90 * t81) * t80; m(6) * (t18 * t78 - t19 * t81); m(6) * (t18 * t81 + t19 * t78); m(6) * (t14 * t11 + t18 * t32 + t19 * t31) + t2 * t116 + t81 * t1 / 0.2e1 + (-t10 * t78 + t9 * t81) * t117 + (-t81 * t4 / 0.2e1 + t3 * t116) * t80; m(6) * (t14 ^ 2 + t18 ^ 2 + t19 ^ 2) + t77 * t15 + (-t81 * t2 - t78 * t1 + t77 * (-t10 * t81 - t78 * t9)) * t80;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t33(1), t33(2), t33(4), t33(7), t33(11); t33(2), t33(3), t33(5), t33(8), t33(12); t33(4), t33(5), t33(6), t33(9), t33(13); t33(7), t33(8), t33(9), t33(10), t33(14); t33(11), t33(12), t33(13), t33(14), t33(15);];
Mq = res;
