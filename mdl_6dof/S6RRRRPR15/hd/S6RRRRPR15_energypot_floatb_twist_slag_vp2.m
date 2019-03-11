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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR15_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:05
% EndTime: 2019-03-10 00:34:06
% DurationCPUTime: 0.79s
% Computational Cost: add. (553->115), mult. (1356->136), div. (0->0), fcn. (1691->14), ass. (0->65)
t43 = cos(pkin(6));
t48 = sin(qJ(1));
t50 = cos(qJ(2));
t74 = t48 * t50;
t47 = sin(qJ(2));
t51 = cos(qJ(1));
t76 = t47 * t51;
t28 = -t43 * t74 - t76;
t40 = sin(pkin(7));
t42 = cos(pkin(7));
t41 = sin(pkin(6));
t80 = t41 * t48;
t63 = -t28 * t40 + t42 * t80;
t79 = t41 * t50;
t62 = -t40 * t79 + t42 * t43;
t92 = -m(1) - m(2);
t91 = -m(7) * pkin(12) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t44 = sin(qJ(6));
t49 = cos(qJ(6));
t90 = -t44 * mrSges(7,1) - t49 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t89 = -m(7) * (pkin(5) + pkin(11)) - t49 * mrSges(7,1) + t44 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t73 = t50 * t51;
t75 = t48 * t47;
t29 = -t43 * t75 + t73;
t46 = sin(qJ(3));
t84 = cos(qJ(3));
t69 = t40 * t84;
t65 = t41 * t69;
t68 = t42 * t84;
t14 = -t28 * t68 + t29 * t46 - t48 * t65;
t87 = pkin(11) * t14;
t81 = t41 * t47;
t19 = -t43 * t69 + t46 * t81 - t68 * t79;
t86 = pkin(11) * t19;
t26 = t43 * t73 - t75;
t27 = t43 * t76 + t74;
t12 = -t26 * t68 + t27 * t46 + t51 * t65;
t85 = t12 * pkin(11);
t83 = cos(qJ(4));
t78 = t41 * t51;
t72 = pkin(8) + r_base(3);
t67 = t43 * pkin(9) + t72;
t66 = t51 * pkin(1) + pkin(9) * t80 + r_base(1);
t64 = -t26 * t40 - t42 * t78;
t61 = t48 * pkin(1) - pkin(9) * t78 + r_base(2);
t60 = t29 * pkin(2) + t63 * pkin(10) + t66;
t59 = pkin(2) * t81 + t62 * pkin(10) + t67;
t15 = t29 * t84 + (t28 * t42 + t40 * t80) * t46;
t58 = t15 * pkin(3) + t60;
t20 = t43 * t40 * t46 + (t42 * t46 * t50 + t47 * t84) * t41;
t57 = t20 * pkin(3) + t59;
t45 = sin(qJ(4));
t5 = t15 * t45 - t63 * t83;
t6 = t15 * t83 + t45 * t63;
t56 = t6 * pkin(4) + qJ(5) * t5 + t58;
t10 = t20 * t45 - t62 * t83;
t11 = t20 * t83 + t45 * t62;
t55 = t11 * pkin(4) + qJ(5) * t10 + t57;
t54 = t27 * pkin(2) + pkin(10) * t64 + t61;
t13 = t27 * t84 + (t26 * t42 - t40 * t78) * t46;
t53 = t13 * pkin(3) + t54;
t3 = t13 * t45 - t64 * t83;
t4 = t13 * t83 + t45 * t64;
t52 = t4 * pkin(4) + t3 * qJ(5) + t53;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t72 - mrSges(2,3) - m(3) * t67 - t43 * mrSges(3,3) - (t47 * mrSges(3,1) + t50 * mrSges(3,2)) * t41 - m(4) * t59 - t20 * mrSges(4,1) - t62 * mrSges(4,3) - m(5) * (t57 + t86) - m(6) * (t55 + t86) - m(7) * t55 + t90 * t10 + t89 * t19 + t91 * t11) * g(3) + (-m(5) * (t53 + t85) - m(6) * (t52 + t85) - m(7) * t52 + mrSges(3,3) * t78 - m(3) * t61 - t64 * mrSges(4,3) - m(4) * t54 - t51 * mrSges(2,2) - t48 * mrSges(2,1) - t26 * mrSges(3,2) - t27 * mrSges(3,1) - t13 * mrSges(4,1) - mrSges(1,2) + t92 * r_base(2) + t91 * t4 + t90 * t3 + t89 * t12) * g(2) + (-m(5) * (t58 + t87) - m(6) * (t56 + t87) - m(7) * t56 - m(3) * t66 - m(4) * t60 - t63 * mrSges(4,3) - mrSges(3,3) * t80 - t51 * mrSges(2,1) + t48 * mrSges(2,2) - t28 * mrSges(3,2) - t29 * mrSges(3,1) - t15 * mrSges(4,1) - mrSges(1,1) + t92 * r_base(1) + t91 * t6 + t90 * t5 + t89 * t14) * g(1);
U  = t1;
