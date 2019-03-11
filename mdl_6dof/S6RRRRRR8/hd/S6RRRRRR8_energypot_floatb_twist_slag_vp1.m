% Calculate potential energy for
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:11
% EndTime: 2019-03-10 04:48:11
% DurationCPUTime: 0.56s
% Computational Cost: add. (547->143), mult. (1191->188), div. (0->0), fcn. (1486->16), ass. (0->65)
t53 = cos(pkin(6));
t58 = sin(qJ(1));
t61 = cos(qJ(2));
t80 = t58 * t61;
t57 = sin(qJ(2));
t62 = cos(qJ(1));
t81 = t57 * t62;
t34 = -t53 * t80 - t81;
t50 = sin(pkin(7));
t52 = cos(pkin(7));
t51 = sin(pkin(6));
t86 = t51 * t58;
t24 = -t34 * t50 + t52 * t86;
t85 = t51 * t61;
t31 = -t50 * t85 + t52 * t53;
t94 = pkin(11) + rSges(5,3);
t93 = pkin(13) + rSges(7,3);
t92 = cos(qJ(3));
t79 = t61 * t62;
t82 = t57 * t58;
t32 = t53 * t79 - t82;
t84 = t51 * t62;
t23 = -t32 * t50 - t52 * t84;
t55 = sin(qJ(4));
t91 = t23 * t55;
t90 = t24 * t55;
t89 = t31 * t55;
t87 = t51 * t57;
t78 = pkin(8) + r_base(3);
t77 = t58 * pkin(1) + r_base(2);
t74 = t50 * t92;
t73 = t52 * t92;
t72 = t53 * pkin(9) + t78;
t71 = t62 * pkin(1) + pkin(9) * t86 + r_base(1);
t70 = t51 * t74;
t35 = -t53 * t82 + t79;
t69 = t35 * pkin(2) + t24 * pkin(10) + t71;
t68 = pkin(2) * t87 + t31 * pkin(10) + t72;
t56 = sin(qJ(3));
t13 = -t34 * t73 + t35 * t56 - t58 * t70;
t14 = t35 * t92 + (t34 * t52 + t50 * t86) * t56;
t60 = cos(qJ(4));
t43 = pkin(4) * t60 + pkin(3);
t63 = -pkin(12) - pkin(11);
t67 = pkin(4) * t90 - t13 * t63 + t14 * t43 + t69;
t21 = -t53 * t74 + t56 * t87 - t73 * t85;
t22 = t53 * t50 * t56 + (t52 * t56 * t61 + t92 * t57) * t51;
t66 = pkin(4) * t89 - t21 * t63 + t22 * t43 + t68;
t33 = t53 * t81 + t80;
t65 = t33 * pkin(2) - pkin(9) * t84 + t23 * pkin(10) + t77;
t11 = -t32 * t73 + t33 * t56 + t62 * t70;
t12 = t33 * t92 + (t32 * t52 - t50 * t84) * t56;
t64 = pkin(4) * t91 - t11 * t63 + t12 * t43 + t65;
t59 = cos(qJ(6));
t54 = sin(qJ(6));
t49 = qJ(4) + qJ(5);
t45 = cos(t49);
t44 = sin(t49);
t8 = t22 * t45 + t31 * t44;
t7 = t22 * t44 - t31 * t45;
t4 = t14 * t45 + t24 * t44;
t3 = t14 * t44 - t24 * t45;
t2 = t12 * t45 + t23 * t44;
t1 = t12 * t44 - t23 * t45;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t62 - rSges(2,2) * t58 + r_base(1)) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t62 + r_base(2)) + g(3) * (rSges(2,3) + t78)) - m(3) * (g(1) * (rSges(3,1) * t35 + rSges(3,2) * t34 + t71) + g(2) * (t33 * rSges(3,1) + t32 * rSges(3,2) + t77) + g(3) * (t53 * rSges(3,3) + t72) + (g(1) * rSges(3,3) * t58 + g(3) * (rSges(3,1) * t57 + rSges(3,2) * t61) + g(2) * (-rSges(3,3) - pkin(9)) * t62) * t51) - m(4) * (g(1) * (rSges(4,1) * t14 - rSges(4,2) * t13 + t24 * rSges(4,3) + t69) + g(2) * (rSges(4,1) * t12 - rSges(4,2) * t11 + rSges(4,3) * t23 + t65) + g(3) * (rSges(4,1) * t22 - rSges(4,2) * t21 + rSges(4,3) * t31 + t68)) - m(5) * (g(1) * (t14 * pkin(3) + (t14 * t60 + t90) * rSges(5,1) + (-t14 * t55 + t24 * t60) * rSges(5,2) + t94 * t13 + t69) + g(2) * (t12 * pkin(3) + (t12 * t60 + t91) * rSges(5,1) + (-t12 * t55 + t23 * t60) * rSges(5,2) + t94 * t11 + t65) + g(3) * (t22 * pkin(3) + (t22 * t60 + t89) * rSges(5,1) + (-t22 * t55 + t31 * t60) * rSges(5,2) + t94 * t21 + t68)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t13 + t67) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t11 + t64) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + rSges(6,3) * t21 + t66)) - m(7) * (g(1) * (t4 * pkin(5) + (t13 * t54 + t4 * t59) * rSges(7,1) + (t13 * t59 - t4 * t54) * rSges(7,2) + t93 * t3 + t67) + g(2) * (t2 * pkin(5) + (t11 * t54 + t2 * t59) * rSges(7,1) + (t11 * t59 - t2 * t54) * rSges(7,2) + t93 * t1 + t64) + g(3) * (t8 * pkin(5) + (t21 * t54 + t59 * t8) * rSges(7,1) + (t21 * t59 - t54 * t8) * rSges(7,2) + t93 * t7 + t66));
U  = t5;
