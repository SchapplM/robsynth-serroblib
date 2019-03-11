% Calculate potential energy for
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:31
% EndTime: 2019-03-10 02:59:32
% DurationCPUTime: 0.50s
% Computational Cost: add. (621->131), mult. (1532->169), div. (0->0), fcn. (1949->14), ass. (0->69)
t52 = cos(pkin(6));
t56 = sin(qJ(2));
t60 = cos(qJ(1));
t84 = t60 * t56;
t57 = sin(qJ(1));
t59 = cos(qJ(2));
t85 = t57 * t59;
t37 = -t52 * t85 - t84;
t49 = sin(pkin(7));
t51 = cos(pkin(7));
t50 = sin(pkin(6));
t90 = t50 * t57;
t71 = -t37 * t49 + t51 * t90;
t89 = t50 * t59;
t70 = -t49 * t89 + t51 * t52;
t98 = rSges(7,1) + pkin(5);
t97 = rSges(7,2) + pkin(12);
t96 = rSges(5,3) + pkin(11);
t95 = rSges(6,3) + pkin(12);
t94 = cos(qJ(3));
t93 = cos(qJ(4));
t91 = t50 * t56;
t88 = t50 * t60;
t86 = t57 * t56;
t83 = t60 * t59;
t82 = rSges(7,3) + qJ(6);
t81 = pkin(8) + r_base(3);
t80 = pkin(1) * t57 + r_base(2);
t77 = t49 * t94;
t76 = t51 * t94;
t75 = pkin(9) * t52 + t81;
t74 = pkin(1) * t60 + pkin(9) * t90 + r_base(1);
t73 = t50 * t77;
t35 = t52 * t83 - t86;
t72 = -t35 * t49 - t51 * t88;
t38 = -t52 * t86 + t83;
t69 = t38 * pkin(2) + pkin(10) * t71 + t74;
t68 = pkin(2) * t91 + pkin(10) * t70 + t75;
t55 = sin(qJ(3));
t23 = t38 * t94 + (t37 * t51 + t49 * t90) * t55;
t67 = pkin(3) * t23 + t69;
t29 = t52 * t49 * t55 + (t51 * t55 * t59 + t56 * t94) * t50;
t66 = pkin(3) * t29 + t68;
t54 = sin(qJ(4));
t12 = t23 * t93 + t54 * t71;
t22 = -t37 * t76 + t38 * t55 - t57 * t73;
t65 = pkin(4) * t12 + t22 * pkin(11) + t67;
t19 = t29 * t93 + t54 * t70;
t28 = -t52 * t77 + t55 * t91 - t76 * t89;
t64 = pkin(4) * t19 + t28 * pkin(11) + t66;
t36 = t52 * t84 + t85;
t63 = pkin(2) * t36 - pkin(9) * t88 + pkin(10) * t72 + t80;
t21 = t36 * t94 + (t35 * t51 - t49 * t88) * t55;
t62 = pkin(3) * t21 + t63;
t10 = t21 * t93 + t54 * t72;
t20 = -t35 * t76 + t36 * t55 + t60 * t73;
t61 = pkin(4) * t10 + t20 * pkin(11) + t62;
t58 = cos(qJ(5));
t53 = sin(qJ(5));
t18 = t29 * t54 - t70 * t93;
t11 = t23 * t54 - t71 * t93;
t9 = t21 * t54 - t72 * t93;
t6 = t19 * t58 + t28 * t53;
t5 = t19 * t53 - t28 * t58;
t4 = t12 * t58 + t22 * t53;
t3 = t12 * t53 - t22 * t58;
t2 = t10 * t58 + t20 * t53;
t1 = t10 * t53 - t20 * t58;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t60 - rSges(2,2) * t57 + r_base(1)) + g(2) * (rSges(2,1) * t57 + rSges(2,2) * t60 + r_base(2)) + g(3) * (rSges(2,3) + t81)) - m(3) * (g(1) * (rSges(3,1) * t38 + rSges(3,2) * t37 + t74) + g(2) * (t36 * rSges(3,1) + t35 * rSges(3,2) + t80) + g(3) * (t52 * rSges(3,3) + t75) + (g(1) * rSges(3,3) * t57 + g(3) * (rSges(3,1) * t56 + rSges(3,2) * t59) + g(2) * (-rSges(3,3) - pkin(9)) * t60) * t50) - m(4) * (g(1) * (t23 * rSges(4,1) - t22 * rSges(4,2) + rSges(4,3) * t71 + t69) + g(2) * (t21 * rSges(4,1) - t20 * rSges(4,2) + rSges(4,3) * t72 + t63) + g(3) * (t29 * rSges(4,1) - t28 * rSges(4,2) + rSges(4,3) * t70 + t68)) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t22 * t96 + t67) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + t20 * t96 + t62) + g(3) * (t19 * rSges(5,1) - t18 * rSges(5,2) + t28 * t96 + t66)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t11 * t95 + t65) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t9 * t95 + t61) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t18 * t95 + t64)) - m(7) * (g(1) * (t11 * t97 + t3 * t82 + t4 * t98 + t65) + g(2) * (t1 * t82 + t2 * t98 + t9 * t97 + t61) + g(3) * (t18 * t97 + t5 * t82 + t6 * t98 + t64));
U  = t7;
