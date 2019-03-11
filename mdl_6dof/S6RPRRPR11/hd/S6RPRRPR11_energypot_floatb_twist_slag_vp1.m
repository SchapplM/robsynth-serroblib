% Calculate potential energy for
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR11_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:38:58
% EndTime: 2019-03-09 05:38:59
% DurationCPUTime: 0.55s
% Computational Cost: add. (584->136), mult. (1388->177), div. (0->0), fcn. (1753->16), ass. (0->63)
t42 = sin(pkin(12));
t48 = cos(pkin(6));
t53 = cos(qJ(1));
t46 = cos(pkin(12));
t52 = sin(qJ(1));
t77 = t52 * t46;
t25 = -t42 * t53 - t48 * t77;
t43 = sin(pkin(7));
t47 = cos(pkin(7));
t44 = sin(pkin(6));
t82 = t44 * t52;
t63 = -t25 * t43 + t47 * t82;
t83 = t44 * t46;
t62 = -t43 * t83 + t47 * t48;
t88 = rSges(5,3) + pkin(10);
t87 = cos(qJ(3));
t86 = cos(qJ(4));
t84 = t42 * t44;
t81 = t44 * t53;
t79 = t48 * t53;
t78 = t52 * t42;
t76 = pkin(11) + qJ(5) + rSges(7,3);
t75 = t44 * qJ(2);
t74 = qJ(5) + rSges(6,3);
t73 = pkin(8) + r_base(3);
t72 = t52 * pkin(1) + r_base(2);
t69 = t43 * t87;
t68 = t47 * t87;
t67 = t48 * qJ(2) + t73;
t66 = t53 * pkin(1) + t52 * t75 + r_base(1);
t65 = t44 * t69;
t23 = t46 * t79 - t78;
t64 = -t23 * t43 - t47 * t81;
t40 = pkin(13) + qJ(6);
t35 = sin(t40);
t36 = cos(t40);
t45 = cos(pkin(13));
t61 = t36 * rSges(7,1) - t35 * rSges(7,2) + pkin(5) * t45 + pkin(4);
t41 = sin(pkin(13));
t60 = t35 * rSges(7,1) + t36 * rSges(7,2) + t41 * pkin(5) + pkin(10);
t26 = t46 * t53 - t48 * t78;
t59 = t26 * pkin(2) + t63 * pkin(9) + t66;
t58 = pkin(2) * t84 + t62 * pkin(9) + t67;
t51 = sin(qJ(3));
t12 = t26 * t87 + (t25 * t47 + t43 * t82) * t51;
t57 = t12 * pkin(3) + t59;
t17 = t48 * t43 * t51 + (t46 * t47 * t51 + t87 * t42) * t44;
t56 = t17 * pkin(3) + t58;
t24 = t42 * t79 + t77;
t55 = t24 * pkin(2) + t64 * pkin(9) - t53 * t75 + t72;
t10 = t24 * t87 + (t23 * t47 - t43 * t81) * t51;
t54 = t10 * pkin(3) + t55;
t50 = sin(qJ(4));
t16 = -t48 * t69 + t51 * t84 - t68 * t83;
t11 = -t25 * t68 + t26 * t51 - t52 * t65;
t9 = -t23 * t68 + t24 * t51 + t53 * t65;
t8 = t17 * t86 + t62 * t50;
t7 = t17 * t50 - t62 * t86;
t4 = t12 * t86 + t63 * t50;
t3 = t12 * t50 - t63 * t86;
t2 = t10 * t86 + t64 * t50;
t1 = t10 * t50 - t64 * t86;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t53 - t52 * rSges(2,2) + r_base(1)) + g(2) * (t52 * rSges(2,1) + rSges(2,2) * t53 + r_base(2)) + g(3) * (rSges(2,3) + t73)) - m(3) * (g(1) * (rSges(3,1) * t26 + rSges(3,2) * t25 + t66) + g(2) * (t24 * rSges(3,1) + t23 * rSges(3,2) + t72) + g(3) * (rSges(3,3) * t48 + t67) + (g(1) * rSges(3,3) * t52 + g(3) * (rSges(3,1) * t42 + rSges(3,2) * t46) + g(2) * (-rSges(3,3) - qJ(2)) * t53) * t44) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t63 * rSges(4,3) + t59) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t64 * rSges(4,3) + t55) + g(3) * (t17 * rSges(4,1) - t16 * rSges(4,2) + t62 * rSges(4,3) + t58)) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + t88 * t11 + t57) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t88 * t9 + t54) + g(3) * (rSges(5,1) * t8 - rSges(5,2) * t7 + t88 * t16 + t56)) - m(6) * (g(1) * (t4 * pkin(4) + t11 * pkin(10) + (t11 * t41 + t4 * t45) * rSges(6,1) + (t11 * t45 - t4 * t41) * rSges(6,2) + t74 * t3 + t57) + g(2) * (t2 * pkin(4) + t9 * pkin(10) + (t2 * t45 + t41 * t9) * rSges(6,1) + (-t2 * t41 + t45 * t9) * rSges(6,2) + t74 * t1 + t54) + g(3) * (t8 * pkin(4) + t16 * pkin(10) + (t16 * t41 + t45 * t8) * rSges(6,1) + (t16 * t45 - t41 * t8) * rSges(6,2) + t74 * t7 + t56)) - m(7) * (g(1) * (t60 * t11 + t76 * t3 + t61 * t4 + t57) + g(2) * (t76 * t1 + t61 * t2 + t60 * t9 + t54) + g(3) * (t60 * t16 + t61 * t8 + t76 * t7 + t56));
U  = t5;
