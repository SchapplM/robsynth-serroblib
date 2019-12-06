% Calculate joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:37
% EndTime: 2019-12-05 17:35:40
% DurationCPUTime: 0.87s
% Computational Cost: add. (1661->206), mult. (2045->291), div. (0->0), fcn. (2125->8), ass. (0->101)
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t48 = -Icges(6,3) * t76 + (Icges(6,5) * t80 - Icges(6,6) * t78) * t75;
t49 = -Icges(5,3) * t76 + (Icges(5,5) * t80 - Icges(5,6) * t78) * t75;
t125 = -t48 - t49;
t50 = -Icges(6,6) * t76 + (Icges(6,4) * t80 - Icges(6,2) * t78) * t75;
t51 = -Icges(5,6) * t76 + (Icges(5,4) * t80 - Icges(5,2) * t78) * t75;
t124 = (-t50 - t51) * t78;
t52 = -Icges(6,5) * t76 + (Icges(6,1) * t80 - Icges(6,4) * t78) * t75;
t53 = -Icges(5,5) * t76 + (Icges(5,1) * t80 - Icges(5,4) * t78) * t75;
t123 = (t52 + t53) * t75 * t80;
t74 = qJ(1) + pkin(7);
t72 = cos(t74);
t112 = t72 * t78;
t71 = sin(t74);
t115 = t71 * t75;
t109 = t76 * t78;
t56 = t71 * t109 + t72 * t80;
t108 = t76 * t80;
t57 = -t71 * t108 + t112;
t77 = -qJ(5) - pkin(6);
t122 = t57 * rSges(6,1) + t56 * rSges(6,2) + pkin(4) * t112 + t77 * t115;
t68 = t80 * pkin(4) + pkin(3);
t117 = pkin(3) - t68;
t121 = t117 * t76;
t120 = t76 ^ 2;
t79 = sin(qJ(1));
t119 = t79 * pkin(1);
t81 = cos(qJ(1));
t118 = t81 * pkin(1);
t116 = pkin(6) + t77;
t114 = t71 * t78;
t113 = t72 * t75;
t107 = t75 * t124 + t125 * t76 + t123;
t89 = Icges(6,3) * t75;
t21 = Icges(6,5) * t57 + Icges(6,6) * t56 - t71 * t89;
t90 = Icges(5,3) * t75;
t23 = Icges(5,5) * t57 + Icges(5,6) * t56 - t71 * t90;
t106 = -t23 - t21;
t58 = -t72 * t109 + t71 * t80;
t59 = t72 * t108 + t114;
t22 = Icges(6,5) * t59 + Icges(6,6) * t58 + t72 * t89;
t24 = Icges(5,5) * t59 + Icges(5,6) * t58 + t72 * t90;
t105 = t24 + t22;
t91 = Icges(6,6) * t75;
t26 = Icges(6,4) * t59 + Icges(6,2) * t58 + t72 * t91;
t92 = Icges(5,6) * t75;
t28 = Icges(5,4) * t59 + Icges(5,2) * t58 + t72 * t92;
t104 = t26 + t28;
t25 = Icges(6,4) * t57 + Icges(6,2) * t56 - t71 * t91;
t27 = Icges(5,4) * t57 + Icges(5,2) * t56 - t71 * t92;
t103 = t27 + t25;
t93 = Icges(6,5) * t75;
t30 = Icges(6,1) * t59 + Icges(6,4) * t58 + t72 * t93;
t94 = Icges(5,5) * t75;
t32 = Icges(5,1) * t59 + Icges(5,4) * t58 + t72 * t94;
t102 = t30 + t32;
t29 = Icges(6,1) * t57 + Icges(6,4) * t56 - t71 * t93;
t31 = Icges(5,1) * t57 + Icges(5,4) * t56 - t71 * t94;
t101 = t31 + t29;
t100 = -(pkin(6) * t75 + t121) * t71 + rSges(6,3) * t115 - t122;
t84 = -t59 * rSges(6,1) - t58 * rSges(6,2);
t99 = pkin(4) * t114 + (-t116 * t75 - t121) * t72 + rSges(6,3) * t113 - t84;
t97 = t57 * rSges(5,1) + t56 * rSges(5,2);
t95 = t71 ^ 2 + t72 ^ 2;
t88 = t72 * qJ(3) - t119;
t87 = -t68 * t76 - pkin(2);
t86 = ((t116 - rSges(6,3)) * t76 + (rSges(6,1) * t80 - rSges(6,2) * t78 - t117) * t75) * t75;
t85 = -t59 * rSges(5,1) - t58 * rSges(5,2);
t83 = -rSges(4,1) * t76 + rSges(4,2) * t75 - pkin(2);
t82 = -pkin(3) * t76 - pkin(2) + (-rSges(5,3) - pkin(6)) * t75;
t65 = -t81 * rSges(2,1) + t79 * rSges(2,2);
t64 = -t79 * rSges(2,1) - t81 * rSges(2,2);
t61 = -t72 * rSges(3,1) + t71 * rSges(3,2) - t118;
t60 = -t71 * rSges(3,1) - t72 * rSges(3,2) - t119;
t55 = -t76 * rSges(5,3) + (rSges(5,1) * t80 - rSges(5,2) * t78) * t75;
t40 = -t118 + (-rSges(4,3) - qJ(3)) * t71 + t83 * t72;
t39 = t72 * rSges(4,3) + t83 * t71 + t88;
t38 = rSges(5,3) * t113 - t85;
t36 = -rSges(5,3) * t115 + t97;
t18 = t55 * t113 + t76 * t38;
t17 = t55 * t115 - t76 * t36;
t16 = -t71 * qJ(3) + t82 * t72 - t118 + t85;
t15 = t82 * t71 + t88 + t97;
t14 = -t118 + (-pkin(4) * t78 - qJ(3)) * t71 + ((-rSges(6,3) + t77) * t75 + t87) * t72 + t84;
t13 = (-rSges(6,3) * t75 + t87) * t71 + t88 + t122;
t12 = t49 * t113 + t58 * t51 + t59 * t53;
t11 = t48 * t113 + t58 * t50 + t59 * t52;
t10 = -t49 * t115 + t56 * t51 + t57 * t53;
t9 = -t48 * t115 + t56 * t50 + t57 * t52;
t8 = (-t36 * t72 - t38 * t71) * t75;
t7 = t72 * t86 + t99 * t76;
t6 = t100 * t76 + t71 * t86;
t5 = -t76 * t24 + (-t28 * t78 + t32 * t80) * t75;
t4 = -t76 * t23 + (-t27 * t78 + t31 * t80) * t75;
t3 = -t76 * t22 + (-t26 * t78 + t30 * t80) * t75;
t2 = -t76 * t21 + (-t25 * t78 + t29 * t80) * t75;
t1 = (t100 * t72 - t99 * t71) * t75;
t19 = [Icges(2,3) + Icges(3,3) + (Icges(4,2) * t76 + t125) * t76 + (Icges(4,1) * t75 + 0.2e1 * Icges(4,4) * t76 + t124) * t75 + m(2) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2) + t123; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t71 * t39 + t72 * t40) + m(5) * (t71 * t15 + t72 * t16) + m(6) * (t71 * t13 + t72 * t14); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t95; -t107 * t76 + m(5) * (t17 * t15 + t18 * t16) + m(6) * (t6 * t13 + t7 * t14) + ((t5 / 0.2e1 + t3 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t72 + (-t4 / 0.2e1 - t2 / 0.2e1 - t10 / 0.2e1 - t9 / 0.2e1) * t71) * t75; m(5) * t8 + m(6) * t1; m(5) * (t17 * t71 + t18 * t72) + m(6) * (t6 * t71 + t7 * t72); m(5) * (t17 ^ 2 + t18 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t6 ^ 2 + t7 ^ 2) + t107 * t120 + (((t102 * t59 + t104 * t58 + t105 * t113) * t113 + (-t11 - t12 - t3 - t5) * t76) * t72 + ((t101 * t57 + t103 * t56 + t106 * t115) * t115 + (t4 + t10 + t2 + t9) * t76 + ((t105 * t71 + t106 * t72) * t75 - t101 * t59 - t103 * t58 - t102 * t57 - t104 * t56) * t113) * t71) * t75; m(6) * (t13 * t72 - t14 * t71) * t75; -m(6) * t76; 0; m(6) * (-t76 * t1 + (t6 * t72 - t7 * t71) * t75); m(6) * (t95 * t75 ^ 2 + t120);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
