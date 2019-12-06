% Calculate joint inertia matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:53
% EndTime: 2019-12-05 15:00:56
% DurationCPUTime: 0.75s
% Computational Cost: add. (2104->164), mult. (2132->272), div. (0->0), fcn. (2258->8), ass. (0->93)
t73 = sin(pkin(7));
t68 = t73 ^ 2;
t75 = cos(pkin(7));
t69 = t75 ^ 2;
t93 = t68 + t69;
t71 = pkin(8) + qJ(3);
t67 = cos(t71);
t109 = t67 ^ 2;
t108 = t73 / 0.2e1;
t107 = -m(5) - m(6);
t65 = sin(t71);
t106 = t65 * t73;
t105 = t65 * t75;
t70 = pkin(9) + qJ(5);
t64 = sin(t70);
t66 = cos(t70);
t36 = -Icges(6,3) * t67 + (Icges(6,5) * t66 - Icges(6,6) * t64) * t65;
t104 = t67 * t36;
t103 = t73 * t64;
t102 = t73 * t66;
t72 = sin(pkin(9));
t101 = t73 * t72;
t74 = cos(pkin(9));
t100 = t73 * t74;
t99 = t75 * t64;
t98 = t75 * t66;
t97 = t75 * t72;
t96 = t75 * t74;
t59 = t65 * pkin(3) - t67 * qJ(4);
t95 = t67 * rSges(5,3) - (rSges(5,1) * t74 - rSges(5,2) * t72) * t65 - t59;
t83 = pkin(3) * t67 + qJ(4) * t65;
t94 = t93 * t83;
t92 = Icges(5,5) * t65;
t91 = Icges(6,5) * t65;
t90 = Icges(5,6) * t65;
t89 = Icges(6,6) * t65;
t88 = Icges(5,3) * t65;
t87 = Icges(6,3) * t65;
t86 = m(5) / 0.2e1 + m(6) / 0.2e1;
t39 = -t67 * rSges(6,3) + (rSges(6,1) * t66 - rSges(6,2) * t64) * t65;
t62 = t74 * pkin(4) + pkin(3);
t76 = -pkin(6) - qJ(4);
t85 = -(qJ(4) + t76) * t67 - (-pkin(3) + t62) * t65 - t39 - t59;
t78 = Icges(4,5) * t67 - Icges(4,6) * t65;
t77 = t62 * t67 - t65 * t76 - t83;
t60 = t65 * rSges(4,1) + t67 * rSges(4,2);
t57 = t67 * t96 + t101;
t56 = -t67 * t97 + t100;
t55 = t67 * t100 - t97;
t54 = -t67 * t101 - t96;
t53 = t67 * t98 + t103;
t52 = -t67 * t99 + t102;
t51 = t67 * t102 - t99;
t50 = -t67 * t103 - t98;
t42 = Icges(4,3) * t73 + t75 * t78;
t41 = -Icges(4,3) * t75 + t73 * t78;
t38 = -Icges(6,5) * t67 + (Icges(6,1) * t66 - Icges(6,4) * t64) * t65;
t37 = -Icges(6,6) * t67 + (Icges(6,4) * t66 - Icges(6,2) * t64) * t65;
t34 = t95 * t75;
t33 = t95 * t73;
t32 = Icges(5,1) * t57 + Icges(5,4) * t56 + t75 * t92;
t31 = Icges(5,1) * t55 + Icges(5,4) * t54 + t73 * t92;
t30 = Icges(5,4) * t57 + Icges(5,2) * t56 + t75 * t90;
t29 = Icges(5,4) * t55 + Icges(5,2) * t54 + t73 * t90;
t28 = Icges(5,5) * t57 + Icges(5,6) * t56 + t75 * t88;
t27 = Icges(5,5) * t55 + Icges(5,6) * t54 + t73 * t88;
t26 = t93 * (rSges(4,1) * t67 - rSges(4,2) * t65);
t25 = t53 * rSges(6,1) + t52 * rSges(6,2) + rSges(6,3) * t105;
t24 = t51 * rSges(6,1) + t50 * rSges(6,2) + rSges(6,3) * t106;
t23 = Icges(6,1) * t53 + Icges(6,4) * t52 + t75 * t91;
t22 = Icges(6,1) * t51 + Icges(6,4) * t50 + t73 * t91;
t21 = Icges(6,4) * t53 + Icges(6,2) * t52 + t75 * t89;
t20 = Icges(6,4) * t51 + Icges(6,2) * t50 + t73 * t89;
t19 = Icges(6,5) * t53 + Icges(6,6) * t52 + t75 * t87;
t18 = Icges(6,5) * t51 + Icges(6,6) * t50 + t73 * t87;
t17 = t85 * t75;
t16 = t85 * t73;
t15 = -t39 * t105 - t67 * t25;
t14 = t39 * t106 + t67 * t24;
t13 = (t24 * t75 - t25 * t73) * t65;
t12 = t73 * (t55 * rSges(5,1) + t54 * rSges(5,2) + rSges(5,3) * t106) + t75 * (t57 * rSges(5,1) + t56 * rSges(5,2) + rSges(5,3) * t105) + t94;
t11 = -t67 * t19 + (-t21 * t64 + t23 * t66) * t65;
t10 = -t67 * t18 + (-t20 * t64 + t22 * t66) * t65;
t9 = t19 * t105 + t52 * t21 + t53 * t23;
t8 = t18 * t105 + t52 * t20 + t53 * t22;
t7 = t19 * t106 + t50 * t21 + t51 * t23;
t6 = t18 * t106 + t50 * t20 + t51 * t22;
t5 = (t75 * t77 + t25) * t75 + (t73 * t77 + t24) * t73 + t94;
t4 = t9 * t73 - t8 * t75;
t3 = -t6 * t75 + t7 * t73;
t2 = -(t52 * t37 + t53 * t38) * t67 + (t8 * t73 + (t9 - t104) * t75) * t65;
t1 = -(t50 * t37 + t51 * t38) * t67 + (t7 * t75 + (t6 - t104) * t73) * t65;
t35 = [m(2) + m(3) + m(4) - t107; 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t86) * t93; m(4) * t26 + m(5) * t12 + m(6) * t5; m(5) * (-t33 * t75 + t34 * t73) + m(6) * (-t16 * t75 + t17 * t73); m(6) * (t16 ^ 2 + t17 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(4) * (t93 * t60 ^ 2 + t26 ^ 2) + (-t69 * t41 - t3 + (t27 * t106 + t54 * t29 + t55 * t31) * t75) * t75 + (t4 + t68 * t42 + (t28 * t105 + t56 * t30 + t57 * t32) * t73 + (-t27 * t105 - t28 * t106 - t56 * t29 - t54 * t30 - t57 * t31 - t55 * t32 - t73 * t41 + t75 * t42) * t75) * t73; t107 * t67; 0; m(6) * (-t67 * t5 + (t16 * t73 + t17 * t75) * t65) + m(5) * (-t67 * t12 + (t33 * t73 + t34 * t75) * t65); 0.2e1 * t86 * (t93 * t65 ^ 2 + t109); m(6) * t13; m(6) * (t14 * t73 - t15 * t75); -t75 * t1 / 0.2e1 + m(6) * (t13 * t5 + t14 * t17 + t15 * t16) - t67 * (-t10 * t75 + t11 * t73) / 0.2e1 + t2 * t108 + (t75 * t4 / 0.2e1 + t3 * t108) * t65; m(6) * (-t13 * t67 + (t14 * t75 + t15 * t73) * t65); m(6) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) + t2 * t105 + t1 * t106 - t67 * (t109 * t36 + (t11 * t75 + t10 * t73 - (-t37 * t64 + t38 * t66) * t67) * t65);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t35(1), t35(2), t35(4), t35(7), t35(11); t35(2), t35(3), t35(5), t35(8), t35(12); t35(4), t35(5), t35(6), t35(9), t35(13); t35(7), t35(8), t35(9), t35(10), t35(14); t35(11), t35(12), t35(13), t35(14), t35(15);];
Mq = res;
