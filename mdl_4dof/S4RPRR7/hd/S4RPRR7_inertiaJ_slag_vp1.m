% Calculate joint inertia matrix for
% S4RPRR7
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR7_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:53:44
% DurationCPUTime: 0.81s
% Computational Cost: add. (1599->165), mult. (2099->273), div. (0->0), fcn. (2214->8), ass. (0->91)
t70 = pkin(7) + qJ(3);
t67 = sin(t70);
t114 = Icges(4,5) * t67;
t113 = t114 / 0.2e1;
t77 = sin(qJ(1));
t71 = t77 ^ 2;
t79 = cos(qJ(1));
t72 = t79 ^ 2;
t98 = t71 + t72;
t68 = cos(t70);
t106 = t68 * t79;
t107 = t67 * t79;
t76 = sin(qJ(4));
t102 = t79 * t76;
t78 = cos(qJ(4));
t103 = t77 * t78;
t50 = -t68 * t102 + t103;
t101 = t79 * t78;
t104 = t77 * t76;
t51 = t68 * t101 + t104;
t28 = t51 * rSges(5,1) + t50 * rSges(5,2) + rSges(5,3) * t107;
t112 = pkin(3) * t106 + pkin(6) * t107 + t28;
t111 = -t68 / 0.2e1;
t110 = t77 / 0.2e1;
t109 = pkin(3) * t68;
t108 = t67 * t77;
t37 = -Icges(5,6) * t68 + (Icges(5,4) * t78 - Icges(5,2) * t76) * t67;
t105 = t76 * t37;
t39 = -t68 * rSges(5,3) + (rSges(5,1) * t78 - rSges(5,2) * t76) * t67;
t100 = -t67 * pkin(3) + t68 * pkin(6) - t39;
t96 = Icges(4,4) * t68;
t95 = Icges(5,5) * t67;
t94 = Icges(5,6) * t67;
t93 = Icges(5,3) * t67;
t92 = rSges(3,3) + qJ(2);
t36 = -Icges(5,3) * t68 + (Icges(5,5) * t78 - Icges(5,6) * t76) * t67;
t38 = -Icges(5,5) * t68 + (Icges(5,1) * t78 - Icges(5,4) * t76) * t67;
t48 = -t68 * t104 - t101;
t49 = t68 * t103 - t102;
t12 = t36 * t108 + t48 * t37 + t49 * t38;
t21 = Icges(5,5) * t49 + Icges(5,6) * t48 + t77 * t93;
t23 = Icges(5,4) * t49 + Icges(5,2) * t48 + t77 * t94;
t25 = Icges(5,1) * t49 + Icges(5,4) * t48 + t77 * t95;
t9 = -t68 * t21 + (-t23 * t76 + t25 * t78) * t67;
t91 = t9 / 0.2e1 + t12 / 0.2e1;
t22 = Icges(5,5) * t51 + Icges(5,6) * t50 + t79 * t93;
t24 = Icges(5,4) * t51 + Icges(5,2) * t50 + t79 * t94;
t26 = Icges(5,1) * t51 + Icges(5,4) * t50 + t79 * t95;
t10 = -t68 * t22 + (-t24 * t76 + t26 * t78) * t67;
t13 = t36 * t107 + t50 * t37 + t51 * t38;
t90 = t10 / 0.2e1 + t13 / 0.2e1;
t74 = cos(pkin(7));
t66 = t74 * pkin(2) + pkin(1);
t75 = -pkin(5) - qJ(2);
t89 = t79 * t66 - t77 * t75;
t88 = rSges(4,1) * t68 - rSges(4,2) * t67;
t87 = -t49 * rSges(5,1) - t48 * rSges(5,2);
t83 = -Icges(4,2) * t67 + t96;
t82 = Icges(4,5) * t68 - Icges(4,6) * t67;
t81 = rSges(4,1) * t106 - rSges(4,2) * t107 + t77 * rSges(4,3);
t73 = sin(pkin(7));
t80 = rSges(3,1) * t74 - rSges(3,2) * t73 + pkin(1);
t59 = t79 * rSges(2,1) - t77 * rSges(2,2);
t58 = -t77 * rSges(2,1) - t79 * rSges(2,2);
t56 = t67 * rSges(4,1) + t68 * rSges(4,2);
t40 = -Icges(4,3) * t79 + t82 * t77;
t35 = t92 * t77 + t80 * t79;
t34 = -t80 * t77 + t92 * t79;
t33 = t67 * t78 * t38;
t32 = t81 + t89;
t31 = (rSges(4,3) - t75) * t79 + (-t66 - t88) * t77;
t30 = t100 * t79;
t29 = t100 * t77;
t27 = rSges(5,3) * t108 - t87;
t20 = t79 * t81 + (-t79 * rSges(4,3) + t88 * t77) * t77;
t19 = t89 + t112;
t18 = -t79 * t75 + (-t109 - t66 + (-rSges(5,3) - pkin(6)) * t67) * t77 + t87;
t17 = -t39 * t107 - t68 * t28;
t16 = t39 * t108 + t68 * t27;
t15 = -t67 * t105 - t68 * t36 + t33;
t14 = (t27 * t79 - t28 * t77) * t67;
t11 = t112 * t79 + (t27 + (pkin(6) * t67 + t109) * t77) * t77;
t8 = t22 * t107 + t50 * t24 + t51 * t26;
t7 = t21 * t107 + t50 * t23 + t51 * t25;
t6 = t22 * t108 + t48 * t24 + t49 * t26;
t5 = t21 * t108 + t48 * t23 + t49 * t25;
t4 = -t7 * t79 + t8 * t77;
t3 = -t5 * t79 + t6 * t77;
t2 = -t13 * t68 + (t7 * t77 + t79 * t8) * t67;
t1 = -t12 * t68 + (t5 * t77 + t6 * t79) * t67;
t41 = [Icges(3,2) * t74 ^ 2 + Icges(2,3) + t33 + (Icges(3,1) * t73 + 0.2e1 * Icges(3,4) * t74) * t73 + (Icges(4,4) * t67 + Icges(4,2) * t68 - t36) * t68 + (Icges(4,1) * t67 - t105 + t96) * t67 + m(2) * (t58 ^ 2 + t59 ^ 2) + m(3) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t31 ^ 2 + t32 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2); m(3) * (t77 * t34 - t79 * t35) + m(4) * (t77 * t31 - t79 * t32) + m(5) * (t77 * t18 - t79 * t19); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t98; ((-Icges(4,6) * t79 + t83 * t77) * t111 + t79 * t113 - t91) * t79 + (t68 * (Icges(4,6) * t77 + t83 * t79) / 0.2e1 + t77 * t113 + t90) * t77 + m(5) * (t30 * t18 + t29 * t19) + m(4) * (-t31 * t79 - t32 * t77) * t56 + (t71 / 0.2e1 + t72 / 0.2e1) * (Icges(4,6) * t68 + t114); m(5) * (-t29 * t79 + t30 * t77); m(4) * (t98 * t56 ^ 2 + t20 ^ 2) + m(5) * (t11 ^ 2 + t29 ^ 2 + t30 ^ 2) + (-t72 * t40 - t3) * t79 + (-t77 * t40 * t79 + t4 + t98 * (Icges(4,3) * t77 + t82 * t79)) * t77; -t15 * t68 + m(5) * (t16 * t18 + t17 * t19) + (t91 * t77 + t90 * t79) * t67; m(5) * (t16 * t77 - t17 * t79); m(5) * (t14 * t11 + t16 * t30 + t17 * t29) + t2 * t110 - t79 * t1 / 0.2e1 + (t10 * t77 - t9 * t79) * t111 + (t79 * t4 / 0.2e1 + t3 * t110) * t67; m(5) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + t68 ^ 2 * t15 + (t79 * t2 + t77 * t1 - t68 * (t10 * t79 + t77 * t9)) * t67;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t41(1), t41(2), t41(4), t41(7); t41(2), t41(3), t41(5), t41(8); t41(4), t41(5), t41(6), t41(9); t41(7), t41(8), t41(9), t41(10);];
Mq = res;
