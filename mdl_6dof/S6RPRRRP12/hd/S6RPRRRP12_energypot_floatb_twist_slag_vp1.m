% Calculate potential energy for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP12_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:43:52
% EndTime: 2019-03-09 06:43:53
% DurationCPUTime: 0.49s
% Computational Cost: add. (621->131), mult. (1532->170), div. (0->0), fcn. (1949->14), ass. (0->70)
t54 = cos(pkin(6));
t49 = sin(pkin(12));
t60 = cos(qJ(1));
t85 = t60 * t49;
t52 = cos(pkin(12));
t58 = sin(qJ(1));
t86 = t58 * t52;
t37 = -t54 * t86 - t85;
t50 = sin(pkin(7));
t53 = cos(pkin(7));
t51 = sin(pkin(6));
t90 = t51 * t58;
t71 = -t37 * t50 + t53 * t90;
t91 = t51 * t52;
t70 = -t50 * t91 + t54 * t53;
t99 = rSges(7,1) + pkin(5);
t98 = rSges(7,2) + pkin(11);
t97 = rSges(5,3) + pkin(10);
t96 = rSges(6,3) + pkin(11);
t95 = cos(qJ(3));
t94 = cos(qJ(4));
t92 = t49 * t51;
t89 = t51 * t60;
t87 = t58 * t49;
t84 = t60 * t52;
t83 = t51 * qJ(2);
t82 = rSges(7,3) + qJ(6);
t81 = pkin(8) + r_base(3);
t80 = t58 * pkin(1) + r_base(2);
t77 = t50 * t95;
t76 = t53 * t95;
t75 = t54 * qJ(2) + t81;
t74 = t60 * pkin(1) + t58 * t83 + r_base(1);
t73 = t51 * t77;
t35 = t54 * t84 - t87;
t72 = -t35 * t50 - t53 * t89;
t38 = -t54 * t87 + t84;
t69 = t38 * pkin(2) + t71 * pkin(9) + t74;
t68 = pkin(2) * t92 + t70 * pkin(9) + t75;
t57 = sin(qJ(3));
t23 = t38 * t95 + (t37 * t53 + t50 * t90) * t57;
t67 = t23 * pkin(3) + t69;
t29 = t54 * t50 * t57 + (t52 * t53 * t57 + t95 * t49) * t51;
t66 = t29 * pkin(3) + t68;
t56 = sin(qJ(4));
t12 = t23 * t94 + t71 * t56;
t22 = -t37 * t76 + t38 * t57 - t58 * t73;
t65 = t12 * pkin(4) + t22 * pkin(10) + t67;
t19 = t29 * t94 + t70 * t56;
t28 = -t54 * t77 + t57 * t92 - t76 * t91;
t64 = t19 * pkin(4) + t28 * pkin(10) + t66;
t36 = t54 * t85 + t86;
t63 = t36 * pkin(2) + t72 * pkin(9) - t60 * t83 + t80;
t21 = t36 * t95 + (t35 * t53 - t50 * t89) * t57;
t62 = t21 * pkin(3) + t63;
t10 = t21 * t94 + t72 * t56;
t20 = -t35 * t76 + t36 * t57 + t60 * t73;
t61 = t10 * pkin(4) + t20 * pkin(10) + t62;
t59 = cos(qJ(5));
t55 = sin(qJ(5));
t18 = t29 * t56 - t70 * t94;
t11 = t23 * t56 - t71 * t94;
t9 = t21 * t56 - t72 * t94;
t6 = t19 * t59 + t28 * t55;
t5 = t19 * t55 - t28 * t59;
t4 = t12 * t59 + t22 * t55;
t3 = t12 * t55 - t22 * t59;
t2 = t10 * t59 + t20 * t55;
t1 = t10 * t55 - t20 * t59;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t60 * rSges(2,1) - t58 * rSges(2,2) + r_base(1)) + g(2) * (t58 * rSges(2,1) + t60 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t81)) - m(3) * (g(1) * (t38 * rSges(3,1) + t37 * rSges(3,2) + t74) + g(2) * (t36 * rSges(3,1) + t35 * rSges(3,2) + t80) + g(3) * (t54 * rSges(3,3) + t75) + (g(1) * rSges(3,3) * t58 + g(3) * (rSges(3,1) * t49 + rSges(3,2) * t52) + g(2) * (-rSges(3,3) - qJ(2)) * t60) * t51) - m(4) * (g(1) * (t23 * rSges(4,1) - t22 * rSges(4,2) + t71 * rSges(4,3) + t69) + g(2) * (t21 * rSges(4,1) - t20 * rSges(4,2) + t72 * rSges(4,3) + t63) + g(3) * (t29 * rSges(4,1) - t28 * rSges(4,2) + t70 * rSges(4,3) + t68)) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t97 * t22 + t67) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t97 * t20 + t62) + g(3) * (t19 * rSges(5,1) - t18 * rSges(5,2) + t97 * t28 + t66)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t96 * t11 + t65) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t96 * t9 + t61) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t96 * t18 + t64)) - m(7) * (g(1) * (t98 * t11 + t82 * t3 + t99 * t4 + t65) + g(2) * (t82 * t1 + t99 * t2 + t98 * t9 + t61) + g(3) * (t98 * t18 + t82 * t5 + t99 * t6 + t64));
U  = t7;
