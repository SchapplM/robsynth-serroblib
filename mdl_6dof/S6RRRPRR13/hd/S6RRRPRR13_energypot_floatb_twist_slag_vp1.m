% Calculate potential energy for
% S6RRRPRR13
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR13_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:51:49
% EndTime: 2019-03-09 19:51:49
% DurationCPUTime: 0.56s
% Computational Cost: add. (547->143), mult. (1191->188), div. (0->0), fcn. (1486->16), ass. (0->65)
t55 = cos(pkin(6));
t60 = sin(qJ(1));
t62 = cos(qJ(2));
t81 = t60 * t62;
t59 = sin(qJ(2));
t63 = cos(qJ(1));
t83 = t59 * t63;
t34 = -t55 * t81 - t83;
t51 = sin(pkin(7));
t54 = cos(pkin(7));
t52 = sin(pkin(6));
t87 = t52 * t60;
t24 = -t34 * t51 + t54 * t87;
t86 = t52 * t62;
t31 = -t51 * t86 + t54 * t55;
t94 = pkin(12) + rSges(7,3);
t93 = cos(qJ(3));
t80 = t62 * t63;
t82 = t60 * t59;
t32 = t55 * t80 - t82;
t85 = t52 * t63;
t23 = -t32 * t51 - t54 * t85;
t50 = sin(pkin(13));
t92 = t23 * t50;
t91 = t24 * t50;
t90 = t31 * t50;
t88 = t52 * t59;
t79 = qJ(4) + rSges(5,3);
t78 = pkin(8) + r_base(3);
t77 = t60 * pkin(1) + r_base(2);
t74 = t51 * t93;
t73 = t54 * t93;
t72 = t55 * pkin(9) + t78;
t71 = t63 * pkin(1) + pkin(9) * t87 + r_base(1);
t70 = t52 * t74;
t35 = -t55 * t82 + t80;
t69 = t35 * pkin(2) + t24 * pkin(10) + t71;
t68 = pkin(2) * t88 + t31 * pkin(10) + t72;
t58 = sin(qJ(3));
t13 = -t34 * t73 + t35 * t58 - t60 * t70;
t14 = t35 * t93 + (t34 * t54 + t51 * t87) * t58;
t53 = cos(pkin(13));
t43 = pkin(4) * t53 + pkin(3);
t56 = -pkin(11) - qJ(4);
t67 = pkin(4) * t91 - t13 * t56 + t14 * t43 + t69;
t21 = -t55 * t74 + t58 * t88 - t73 * t86;
t22 = t55 * t51 * t58 + (t54 * t58 * t62 + t93 * t59) * t52;
t66 = pkin(4) * t90 - t21 * t56 + t22 * t43 + t68;
t33 = t55 * t83 + t81;
t65 = t33 * pkin(2) - pkin(9) * t85 + t23 * pkin(10) + t77;
t11 = -t32 * t73 + t33 * t58 + t63 * t70;
t12 = t33 * t93 + (t32 * t54 - t51 * t85) * t58;
t64 = pkin(4) * t92 - t11 * t56 + t12 * t43 + t65;
t61 = cos(qJ(6));
t57 = sin(qJ(6));
t49 = pkin(13) + qJ(5);
t45 = cos(t49);
t44 = sin(t49);
t8 = t22 * t45 + t31 * t44;
t7 = t22 * t44 - t31 * t45;
t4 = t14 * t45 + t24 * t44;
t3 = t14 * t44 - t24 * t45;
t2 = t12 * t45 + t23 * t44;
t1 = t12 * t44 - t23 * t45;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t63 - t60 * rSges(2,2) + r_base(1)) + g(2) * (t60 * rSges(2,1) + rSges(2,2) * t63 + r_base(2)) + g(3) * (rSges(2,3) + t78)) - m(3) * (g(1) * (rSges(3,1) * t35 + rSges(3,2) * t34 + t71) + g(2) * (t33 * rSges(3,1) + t32 * rSges(3,2) + t77) + g(3) * (t55 * rSges(3,3) + t72) + (g(1) * rSges(3,3) * t60 + g(3) * (rSges(3,1) * t59 + rSges(3,2) * t62) + g(2) * (-rSges(3,3) - pkin(9)) * t63) * t52) - m(4) * (g(1) * (rSges(4,1) * t14 - rSges(4,2) * t13 + rSges(4,3) * t24 + t69) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t23 * rSges(4,3) + t65) + g(3) * (rSges(4,1) * t22 - rSges(4,2) * t21 + rSges(4,3) * t31 + t68)) - m(5) * (g(1) * (t14 * pkin(3) + (t14 * t53 + t91) * rSges(5,1) + (-t14 * t50 + t24 * t53) * rSges(5,2) + t79 * t13 + t69) + g(2) * (t12 * pkin(3) + (t12 * t53 + t92) * rSges(5,1) + (-t12 * t50 + t23 * t53) * rSges(5,2) + t79 * t11 + t65) + g(3) * (t22 * pkin(3) + (t22 * t53 + t90) * rSges(5,1) + (-t22 * t50 + t31 * t53) * rSges(5,2) + t79 * t21 + t68)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t13 + t67) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t11 * rSges(6,3) + t64) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + rSges(6,3) * t21 + t66)) - m(7) * (g(1) * (t4 * pkin(5) + (t13 * t57 + t4 * t61) * rSges(7,1) + (t13 * t61 - t4 * t57) * rSges(7,2) + t94 * t3 + t67) + g(2) * (t2 * pkin(5) + (t11 * t57 + t2 * t61) * rSges(7,1) + (t11 * t61 - t2 * t57) * rSges(7,2) + t94 * t1 + t64) + g(3) * (t8 * pkin(5) + (t21 * t57 + t61 * t8) * rSges(7,1) + (t21 * t61 - t57 * t8) * rSges(7,2) + t94 * t7 + t66));
U  = t5;
