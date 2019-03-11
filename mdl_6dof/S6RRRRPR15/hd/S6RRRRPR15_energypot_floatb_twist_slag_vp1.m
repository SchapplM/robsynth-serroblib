% Calculate potential energy for
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR15_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:05
% EndTime: 2019-03-10 00:34:05
% DurationCPUTime: 0.50s
% Computational Cost: add. (553->126), mult. (1342->158), div. (0->0), fcn. (1691->14), ass. (0->64)
t43 = cos(pkin(6));
t47 = sin(qJ(2));
t51 = cos(qJ(1));
t77 = t51 * t47;
t48 = sin(qJ(1));
t50 = cos(qJ(2));
t78 = t48 * t50;
t28 = -t43 * t78 - t77;
t40 = sin(pkin(7));
t42 = cos(pkin(7));
t41 = sin(pkin(6));
t83 = t41 * t48;
t64 = -t28 * t40 + t42 * t83;
t82 = t41 * t50;
t63 = -t40 * t82 + t43 * t42;
t90 = rSges(6,1) + pkin(11);
t89 = rSges(5,3) + pkin(11);
t88 = pkin(12) + rSges(7,3);
t87 = cos(qJ(3));
t86 = cos(qJ(4));
t84 = t41 * t47;
t81 = t41 * t51;
t79 = t48 * t47;
t76 = t51 * t50;
t75 = rSges(6,3) + qJ(5);
t74 = pkin(8) + r_base(3);
t73 = t48 * pkin(1) + r_base(2);
t70 = t40 * t87;
t69 = t42 * t87;
t68 = t43 * pkin(9) + t74;
t67 = t51 * pkin(1) + pkin(9) * t83 + r_base(1);
t66 = t41 * t70;
t26 = t43 * t76 - t79;
t65 = -t26 * t40 - t42 * t81;
t44 = sin(qJ(6));
t49 = cos(qJ(6));
t62 = rSges(7,1) * t44 + rSges(7,2) * t49 + qJ(5);
t61 = rSges(7,1) * t49 - t44 * rSges(7,2) + pkin(5) + pkin(11);
t29 = -t43 * t79 + t76;
t60 = t29 * pkin(2) + t64 * pkin(10) + t67;
t59 = pkin(2) * t84 + t63 * pkin(10) + t68;
t46 = sin(qJ(3));
t15 = t29 * t87 + (t28 * t42 + t40 * t83) * t46;
t58 = t15 * pkin(3) + t60;
t45 = sin(qJ(4));
t6 = t15 * t86 + t64 * t45;
t57 = t6 * pkin(4) + t58;
t20 = t43 * t40 * t46 + (t42 * t46 * t50 + t87 * t47) * t41;
t56 = t20 * pkin(3) + t59;
t11 = t20 * t86 + t63 * t45;
t55 = t11 * pkin(4) + t56;
t27 = t43 * t77 + t78;
t54 = t27 * pkin(2) - pkin(9) * t81 + t65 * pkin(10) + t73;
t13 = t27 * t87 + (t26 * t42 - t40 * t81) * t46;
t53 = t13 * pkin(3) + t54;
t4 = t13 * t86 + t65 * t45;
t52 = t4 * pkin(4) + t53;
t19 = -t43 * t70 + t46 * t84 - t69 * t82;
t14 = -t28 * t69 + t29 * t46 - t48 * t66;
t12 = -t26 * t69 + t27 * t46 + t51 * t66;
t10 = t20 * t45 - t63 * t86;
t5 = t15 * t45 - t64 * t86;
t3 = t13 * t45 - t65 * t86;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t51 - t48 * rSges(2,2) + r_base(1)) + g(2) * (t48 * rSges(2,1) + rSges(2,2) * t51 + r_base(2)) + g(3) * (rSges(2,3) + t74)) - m(3) * (g(1) * (rSges(3,1) * t29 + rSges(3,2) * t28 + t67) + g(2) * (t27 * rSges(3,1) + t26 * rSges(3,2) + t73) + g(3) * (t43 * rSges(3,3) + t68) + (g(1) * rSges(3,3) * t48 + g(3) * (rSges(3,1) * t47 + rSges(3,2) * t50) + g(2) * (-rSges(3,3) - pkin(9)) * t51) * t41) - m(4) * (g(1) * (t15 * rSges(4,1) - t14 * rSges(4,2) + t64 * rSges(4,3) + t60) + g(2) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t65 * rSges(4,3) + t54) + g(3) * (t20 * rSges(4,1) - t19 * rSges(4,2) + t63 * rSges(4,3) + t59)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t89 * t14 + t58) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t89 * t12 + t53) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) + t89 * t19 + t56)) - m(6) * (g(1) * (-t6 * rSges(6,2) + t90 * t14 + t75 * t5 + t57) + g(2) * (-t4 * rSges(6,2) + t90 * t12 + t75 * t3 + t52) + g(3) * (-t11 * rSges(6,2) + t75 * t10 + t90 * t19 + t55)) - m(7) * (g(1) * (t61 * t14 + t62 * t5 + t88 * t6 + t57) + g(2) * (t61 * t12 + t62 * t3 + t88 * t4 + t52) + g(3) * (t62 * t10 + t88 * t11 + t61 * t19 + t55));
U  = t1;
