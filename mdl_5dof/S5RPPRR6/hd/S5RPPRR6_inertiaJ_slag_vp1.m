% Calculate joint inertia matrix for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:43
% DurationCPUTime: 0.84s
% Computational Cost: add. (2585->182), mult. (2211->285), div. (0->0), fcn. (2314->10), ass. (0->98)
t77 = pkin(9) + qJ(4);
t72 = sin(t77);
t121 = Icges(5,5) * t72;
t120 = t121 / 0.2e1;
t78 = qJ(1) + pkin(8);
t73 = sin(t78);
t70 = t73 ^ 2;
t75 = cos(t78);
t71 = t75 ^ 2;
t104 = t70 + t71;
t74 = cos(t77);
t110 = t74 * t75;
t113 = t72 * t75;
t82 = sin(qJ(5));
t109 = t75 * t82;
t84 = cos(qJ(5));
t111 = t73 * t84;
t50 = -t74 * t109 + t111;
t108 = t75 * t84;
t112 = t73 * t82;
t51 = t74 * t108 + t112;
t28 = t51 * rSges(6,1) + t50 * rSges(6,2) + rSges(6,3) * t113;
t119 = pkin(4) * t110 + pkin(7) * t113 + t28;
t118 = t73 / 0.2e1;
t117 = -t74 / 0.2e1;
t116 = pkin(4) * t74;
t83 = sin(qJ(1));
t115 = t83 * pkin(1);
t114 = t72 * t73;
t43 = -Icges(6,6) * t74 + (Icges(6,4) * t84 - Icges(6,2) * t82) * t72;
t107 = t82 * t43;
t45 = -t74 * rSges(6,3) + (rSges(6,1) * t84 - rSges(6,2) * t82) * t72;
t106 = -t72 * pkin(4) + t74 * pkin(7) - t45;
t102 = Icges(5,4) * t74;
t101 = Icges(6,5) * t72;
t100 = Icges(6,6) * t72;
t99 = Icges(6,3) * t72;
t98 = rSges(4,3) + qJ(3);
t42 = -Icges(6,3) * t74 + (Icges(6,5) * t84 - Icges(6,6) * t82) * t72;
t44 = -Icges(6,5) * t74 + (Icges(6,1) * t84 - Icges(6,4) * t82) * t72;
t48 = -t74 * t112 - t108;
t49 = t74 * t111 - t109;
t13 = t42 * t114 + t48 * t43 + t49 * t44;
t21 = Icges(6,5) * t49 + Icges(6,6) * t48 + t73 * t99;
t23 = Icges(6,4) * t49 + Icges(6,2) * t48 + t73 * t100;
t25 = Icges(6,1) * t49 + Icges(6,4) * t48 + t73 * t101;
t9 = -t74 * t21 + (-t23 * t82 + t25 * t84) * t72;
t97 = t9 / 0.2e1 + t13 / 0.2e1;
t22 = Icges(6,5) * t51 + Icges(6,6) * t50 + t75 * t99;
t24 = Icges(6,4) * t51 + Icges(6,2) * t50 + t75 * t100;
t26 = Icges(6,1) * t51 + Icges(6,4) * t50 + t75 * t101;
t10 = -t74 * t22 + (-t24 * t82 + t26 * t84) * t72;
t14 = t42 * t113 + t50 * t43 + t51 * t44;
t96 = t10 / 0.2e1 + t14 / 0.2e1;
t80 = cos(pkin(9));
t69 = t80 * pkin(3) + pkin(2);
t85 = cos(qJ(1));
t76 = t85 * pkin(1);
t81 = -pkin(6) - qJ(3);
t95 = t75 * t69 - t73 * t81 + t76;
t94 = rSges(5,1) * t74 - rSges(5,2) * t72;
t93 = -t49 * rSges(6,1) - t48 * rSges(6,2);
t89 = -Icges(5,2) * t72 + t102;
t88 = Icges(5,5) * t74 - Icges(5,6) * t72;
t87 = rSges(5,1) * t110 - rSges(5,2) * t113 + t73 * rSges(5,3);
t79 = sin(pkin(9));
t86 = rSges(4,1) * t80 - rSges(4,2) * t79 + pkin(2);
t67 = t85 * rSges(2,1) - t83 * rSges(2,2);
t66 = -t83 * rSges(2,1) - t85 * rSges(2,2);
t58 = t72 * rSges(5,1) + t74 * rSges(5,2);
t53 = t75 * rSges(3,1) - t73 * rSges(3,2) + t76;
t52 = -t73 * rSges(3,1) - t75 * rSges(3,2) - t115;
t36 = -Icges(5,3) * t75 + t88 * t73;
t35 = t72 * t84 * t44;
t34 = t98 * t73 + t86 * t75 + t76;
t33 = -t86 * t73 + t98 * t75 - t115;
t32 = t106 * t75;
t31 = t106 * t73;
t30 = t87 + t95;
t29 = -t115 + (rSges(5,3) - t81) * t75 + (-t69 - t94) * t73;
t27 = rSges(6,3) * t114 - t93;
t20 = t75 * t87 + (-t75 * rSges(5,3) + t94 * t73) * t73;
t19 = -t72 * t107 - t74 * t42 + t35;
t18 = -t45 * t113 - t74 * t28;
t17 = t45 * t114 + t74 * t27;
t16 = t95 + t119;
t15 = -t115 - t75 * t81 + (-t116 - t69 + (-rSges(6,3) - pkin(7)) * t72) * t73 + t93;
t12 = (t27 * t75 - t28 * t73) * t72;
t11 = t119 * t75 + (t27 + (pkin(7) * t72 + t116) * t73) * t73;
t8 = t22 * t113 + t50 * t24 + t51 * t26;
t7 = t21 * t113 + t50 * t23 + t51 * t25;
t6 = t22 * t114 + t48 * t24 + t49 * t26;
t5 = t21 * t114 + t48 * t23 + t49 * t25;
t4 = -t7 * t75 + t8 * t73;
t3 = -t5 * t75 + t6 * t73;
t2 = -t14 * t74 + (t7 * t73 + t75 * t8) * t72;
t1 = -t13 * t74 + (t5 * t73 + t6 * t75) * t72;
t37 = [Icges(4,2) * t80 ^ 2 + Icges(2,3) + Icges(3,3) + t35 + (Icges(4,1) * t79 + 0.2e1 * Icges(4,4) * t80) * t79 + (Icges(5,4) * t72 + Icges(5,2) * t74 - t42) * t74 + (Icges(5,1) * t72 + t102 - t107) * t72 + m(6) * (t15 ^ 2 + t16 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(2) * (t66 ^ 2 + t67 ^ 2) + m(3) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(6) * (t73 * t15 - t75 * t16) + m(5) * (t73 * t29 - t75 * t30) + m(4) * (t73 * t33 - t75 * t34); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t104; ((-Icges(5,6) * t75 + t89 * t73) * t117 + t75 * t120 - t97) * t75 + (t74 * (Icges(5,6) * t73 + t89 * t75) / 0.2e1 + t73 * t120 + t96) * t73 + m(6) * (t32 * t15 + t31 * t16) + m(5) * (-t29 * t75 - t30 * t73) * t58 + (t70 / 0.2e1 + t71 / 0.2e1) * (Icges(5,6) * t74 + t121); m(5) * t20 + m(6) * t11; m(6) * (-t31 * t75 + t32 * t73); m(5) * (t104 * t58 ^ 2 + t20 ^ 2) + m(6) * (t11 ^ 2 + t31 ^ 2 + t32 ^ 2) + (-t71 * t36 - t3) * t75 + (-t73 * t36 * t75 + t4 + t104 * (Icges(5,3) * t73 + t88 * t75)) * t73; m(6) * (t17 * t15 + t18 * t16) - t19 * t74 + (t97 * t73 + t96 * t75) * t72; m(6) * t12; m(6) * (t17 * t73 - t18 * t75); m(6) * (t12 * t11 + t17 * t32 + t18 * t31) + t2 * t118 - t75 * t1 / 0.2e1 + (t10 * t73 - t9 * t75) * t117 + (t75 * t4 / 0.2e1 + t3 * t118) * t72; m(6) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + t74 ^ 2 * t19 + (t75 * t2 + t73 * t1 - t74 * (t10 * t75 + t73 * t9)) * t72;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t37(1), t37(2), t37(4), t37(7), t37(11); t37(2), t37(3), t37(5), t37(8), t37(12); t37(4), t37(5), t37(6), t37(9), t37(13); t37(7), t37(8), t37(9), t37(10), t37(14); t37(11), t37(12), t37(13), t37(14), t37(15);];
Mq = res;
