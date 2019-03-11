% Calculate potential energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:25
% EndTime: 2019-03-09 18:56:26
% DurationCPUTime: 0.82s
% Computational Cost: add. (619->155), mult. (1464->210), div. (0->0), fcn. (1852->16), ass. (0->74)
t48 = sin(pkin(7));
t51 = cos(pkin(7));
t55 = sin(qJ(3));
t60 = cos(qJ(3));
t64 = (rSges(4,1) * t55 + rSges(4,2) * t60) * t48 + t51 * pkin(10);
t57 = sin(qJ(1));
t95 = g(1) * t57;
t62 = cos(qJ(1));
t94 = g(2) * t62;
t92 = -rSges(6,3) - pkin(11);
t91 = pkin(12) + rSges(7,3);
t52 = cos(pkin(6));
t61 = cos(qJ(2));
t79 = t61 * t62;
t56 = sin(qJ(2));
t82 = t57 * t56;
t35 = t52 * t79 - t82;
t90 = t35 * t48;
t81 = t57 * t61;
t83 = t56 * t62;
t37 = -t52 * t81 - t83;
t89 = t37 * t48;
t49 = sin(pkin(6));
t88 = t49 * t57;
t87 = t49 * t62;
t86 = t51 * t55;
t85 = t51 * t60;
t84 = t51 * t61;
t80 = t61 * t48;
t78 = pkin(10) + qJ(4);
t77 = pkin(8) + r_base(3);
t76 = t57 * pkin(1) + r_base(2);
t75 = t52 * pkin(9) + t77;
t74 = t62 * pkin(1) + pkin(9) * t88 + r_base(1);
t47 = sin(pkin(13));
t50 = cos(pkin(13));
t72 = t47 * t60 + t50 * t55;
t40 = -t47 * t55 + t50 * t60;
t32 = pkin(3) * t48 * t55 + t78 * t51;
t33 = pkin(3) * t86 - t78 * t48;
t43 = pkin(3) * t60 + pkin(2);
t71 = t52 * t32 + t75 + (t33 * t61 + t43 * t56) * t49;
t38 = -t52 * t82 + t79;
t70 = t32 * t88 + t37 * t33 + t38 * t43 + t74;
t69 = t40 * t48;
t29 = t72 * t48;
t31 = t72 * t51;
t15 = t29 * t52 + (t31 * t61 + t40 * t56) * t49;
t68 = t15 * pkin(4) + t71;
t12 = t29 * t88 + t31 * t37 + t38 * t40;
t67 = t12 * pkin(4) + t70;
t66 = t49 * t69;
t36 = t52 * t83 + t81;
t65 = t35 * t33 + t36 * t43 + (-pkin(9) - t32) * t87 + t76;
t10 = -t29 * t87 + t35 * t31 + t36 * t40;
t63 = t10 * pkin(4) + t65;
t59 = cos(qJ(5));
t58 = cos(qJ(6));
t54 = sin(qJ(5));
t53 = sin(qJ(6));
t34 = -t49 * t80 + t51 * t52;
t30 = t40 * t51;
t21 = t51 * t88 - t89;
t20 = -t51 * t87 - t90;
t14 = (t30 * t61 - t56 * t72) * t49 + t52 * t69;
t11 = t37 * t30 - t38 * t72 + t57 * t66;
t9 = t30 * t35 - t36 * t72 - t62 * t66;
t8 = t15 * t59 + t34 * t54;
t7 = t15 * t54 - t34 * t59;
t4 = t12 * t59 + t21 * t54;
t3 = t12 * t54 - t21 * t59;
t2 = t10 * t59 + t20 * t54;
t1 = t10 * t54 - t20 * t59;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t62 - t57 * rSges(2,2) + r_base(1)) + g(2) * (t57 * rSges(2,1) + rSges(2,2) * t62 + r_base(2)) + g(3) * (rSges(2,3) + t77)) - m(3) * (g(1) * (rSges(3,1) * t38 + rSges(3,2) * t37 + t74) + g(2) * (t36 * rSges(3,1) + t35 * rSges(3,2) + t76) + g(3) * (rSges(3,3) * t52 + t75) + (rSges(3,3) * t95 + g(3) * (rSges(3,1) * t56 + rSges(3,2) * t61) + (-rSges(3,3) - pkin(9)) * t94) * t49) - m(4) * (g(1) * (t38 * pkin(2) - pkin(10) * t89 + (t37 * t86 + t38 * t60) * rSges(4,1) + (t37 * t85 - t38 * t55) * rSges(4,2) + t21 * rSges(4,3) + t74) + g(2) * (t36 * pkin(2) - pkin(10) * t90 + (t35 * t86 + t36 * t60) * rSges(4,1) + (t35 * t85 - t36 * t55) * rSges(4,2) + t20 * rSges(4,3) + t76) + (t64 * t95 + (-pkin(9) - t64) * t94) * t49 + (t34 * rSges(4,3) + t75 + t64 * t52 + (t56 * pkin(2) - pkin(10) * t80 + (t55 * t84 + t56 * t60) * rSges(4,1) + (-t55 * t56 + t60 * t84) * rSges(4,2)) * t49) * g(3)) - m(5) * (g(1) * (rSges(5,1) * t12 + rSges(5,2) * t11 + rSges(5,3) * t21 + t70) + g(2) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t20 * rSges(5,3) + t65) + g(3) * (rSges(5,1) * t15 + rSges(5,2) * t14 + rSges(5,3) * t34 + t71)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t92 * t11 + t67) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t92 * t9 + t63) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + t92 * t14 + t68)) - m(7) * (g(1) * (t4 * pkin(5) - t11 * pkin(11) + (-t11 * t53 + t4 * t58) * rSges(7,1) + (-t11 * t58 - t4 * t53) * rSges(7,2) + t91 * t3 + t67) + g(2) * (t2 * pkin(5) - t9 * pkin(11) + (t2 * t58 - t53 * t9) * rSges(7,1) + (-t2 * t53 - t58 * t9) * rSges(7,2) + t91 * t1 + t63) + g(3) * (t8 * pkin(5) - t14 * pkin(11) + (-t14 * t53 + t58 * t8) * rSges(7,1) + (-t14 * t58 - t53 * t8) * rSges(7,2) + t91 * t7 + t68));
U  = t5;
