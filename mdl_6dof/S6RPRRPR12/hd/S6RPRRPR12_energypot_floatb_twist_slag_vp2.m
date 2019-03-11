% Calculate potential energy for
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:07
% EndTime: 2019-03-09 05:48:08
% DurationCPUTime: 0.75s
% Computational Cost: add. (553->115), mult. (1356->138), div. (0->0), fcn. (1691->14), ass. (0->65)
t40 = sin(pkin(12));
t45 = cos(pkin(6));
t51 = cos(qJ(1));
t43 = cos(pkin(12));
t49 = sin(qJ(1));
t74 = t49 * t43;
t28 = -t40 * t51 - t45 * t74;
t41 = sin(pkin(7));
t44 = cos(pkin(7));
t42 = sin(pkin(6));
t79 = t42 * t49;
t63 = -t28 * t41 + t44 * t79;
t80 = t42 * t43;
t62 = -t41 * t80 + t44 * t45;
t92 = -m(1) - m(2);
t91 = -m(7) * pkin(11) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t46 = sin(qJ(6));
t50 = cos(qJ(6));
t90 = -t46 * mrSges(7,1) - t50 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t89 = -m(7) * (pkin(5) + pkin(10)) - t50 * mrSges(7,1) + t46 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t75 = t49 * t40;
t29 = t43 * t51 - t45 * t75;
t48 = sin(qJ(3));
t84 = cos(qJ(3));
t69 = t41 * t84;
t65 = t42 * t69;
t68 = t44 * t84;
t14 = -t28 * t68 + t29 * t48 - t49 * t65;
t87 = pkin(10) * t14;
t81 = t40 * t42;
t19 = -t45 * t69 + t48 * t81 - t68 * t80;
t86 = pkin(10) * t19;
t76 = t45 * t51;
t26 = t43 * t76 - t75;
t27 = t40 * t76 + t74;
t12 = -t26 * t68 + t27 * t48 + t51 * t65;
t85 = t12 * pkin(10);
t83 = cos(qJ(4));
t78 = t42 * t51;
t73 = qJ(2) * t42;
t72 = pkin(8) + r_base(3);
t67 = t45 * qJ(2) + t72;
t66 = t51 * pkin(1) + t49 * t73 + r_base(1);
t64 = -t26 * t41 - t44 * t78;
t61 = t49 * pkin(1) - t51 * t73 + r_base(2);
t60 = t29 * pkin(2) + t63 * pkin(9) + t66;
t59 = pkin(2) * t81 + t62 * pkin(9) + t67;
t15 = t29 * t84 + (t28 * t44 + t41 * t79) * t48;
t58 = t15 * pkin(3) + t60;
t20 = t45 * t41 * t48 + (t43 * t44 * t48 + t40 * t84) * t42;
t57 = t20 * pkin(3) + t59;
t47 = sin(qJ(4));
t5 = t15 * t47 - t63 * t83;
t6 = t15 * t83 + t47 * t63;
t56 = t6 * pkin(4) + qJ(5) * t5 + t58;
t10 = t20 * t47 - t62 * t83;
t11 = t20 * t83 + t47 * t62;
t55 = t11 * pkin(4) + qJ(5) * t10 + t57;
t54 = t27 * pkin(2) + pkin(9) * t64 + t61;
t13 = t27 * t84 + (t26 * t44 - t41 * t78) * t48;
t53 = t13 * pkin(3) + t54;
t3 = t13 * t47 - t64 * t83;
t4 = t13 * t83 + t47 * t64;
t52 = t4 * pkin(4) + t3 * qJ(5) + t53;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t72 - mrSges(2,3) - m(3) * t67 - t45 * mrSges(3,3) - (t40 * mrSges(3,1) + t43 * mrSges(3,2)) * t42 - m(4) * t59 - t20 * mrSges(4,1) - t62 * mrSges(4,3) - m(5) * (t57 + t86) - m(6) * (t55 + t86) - m(7) * t55 + t90 * t10 + t89 * t19 + t91 * t11) * g(3) + (-m(7) * t52 - m(5) * (t53 + t85) - m(6) * (t52 + t85) + mrSges(3,3) * t78 - t64 * mrSges(4,3) - m(3) * t61 - m(4) * t54 - t49 * mrSges(2,1) - t51 * mrSges(2,2) - t26 * mrSges(3,2) - t27 * mrSges(3,1) - t13 * mrSges(4,1) - mrSges(1,2) + t92 * r_base(2) + t91 * t4 + t90 * t3 + t89 * t12) * g(2) + (-m(7) * t56 - m(5) * (t58 + t87) - m(6) * (t56 + t87) - m(3) * t66 - t63 * mrSges(4,3) - m(4) * t60 - mrSges(3,3) * t79 + t49 * mrSges(2,2) - t51 * mrSges(2,1) - t29 * mrSges(3,1) - t28 * mrSges(3,2) - t15 * mrSges(4,1) - mrSges(1,1) + t92 * r_base(1) + t91 * t6 + t90 * t5 + t89 * t14) * g(1);
U  = t1;
