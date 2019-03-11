% Calculate potential energy for
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:48:55
% EndTime: 2019-03-08 18:48:56
% DurationCPUTime: 0.78s
% Computational Cost: add. (553->115), mult. (1356->140), div. (0->0), fcn. (1691->14), ass. (0->65)
t42 = sin(pkin(7));
t46 = cos(pkin(7));
t47 = cos(pkin(6));
t43 = sin(pkin(6));
t44 = cos(pkin(12));
t78 = t43 * t44;
t62 = -t42 * t78 + t46 * t47;
t40 = sin(pkin(12));
t45 = cos(pkin(11));
t41 = sin(pkin(11));
t79 = t41 * t47;
t28 = -t40 * t45 - t44 * t79;
t76 = t43 * t46;
t63 = -t28 * t42 + t41 * t76;
t92 = -m(1) - m(2);
t91 = -m(7) * pkin(10) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t48 = sin(qJ(6));
t51 = cos(qJ(6));
t90 = -t48 * mrSges(7,1) - t51 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t89 = -m(7) * (pkin(5) + pkin(9)) - t51 * mrSges(7,1) + t48 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t75 = t45 * t47;
t26 = -t40 * t41 + t44 * t75;
t27 = t40 * t75 + t41 * t44;
t50 = sin(qJ(3));
t84 = cos(qJ(3));
t69 = t42 * t84;
t65 = t43 * t69;
t68 = t46 * t84;
t10 = -t26 * t68 + t27 * t50 + t45 * t65;
t87 = pkin(9) * t10;
t29 = -t40 * t79 + t44 * t45;
t12 = -t28 * t68 + t29 * t50 - t41 * t65;
t86 = pkin(9) * t12;
t81 = t40 * t43;
t19 = -t47 * t69 + t50 * t81 - t68 * t78;
t85 = pkin(9) * t19;
t83 = cos(qJ(4));
t80 = t41 * t43;
t77 = t43 * t45;
t73 = qJ(2) * t43;
t70 = qJ(1) + r_base(3);
t67 = t45 * pkin(1) + t41 * t73 + r_base(1);
t66 = t47 * qJ(2) + t70;
t64 = -t26 * t42 - t45 * t76;
t61 = t41 * pkin(1) - t45 * t73 + r_base(2);
t60 = t29 * pkin(2) + t63 * pkin(8) + t67;
t13 = t29 * t84 + (t28 * t46 + t42 * t80) * t50;
t59 = t13 * pkin(3) + t60;
t58 = pkin(2) * t81 + t62 * pkin(8) + t66;
t20 = t47 * t42 * t50 + (t44 * t46 * t50 + t40 * t84) * t43;
t57 = t20 * pkin(3) + t58;
t49 = sin(qJ(4));
t5 = t13 * t49 - t63 * t83;
t6 = t13 * t83 + t49 * t63;
t56 = t6 * pkin(4) + qJ(5) * t5 + t59;
t14 = t20 * t49 - t62 * t83;
t15 = t20 * t83 + t49 * t62;
t55 = t15 * pkin(4) + qJ(5) * t14 + t57;
t54 = t27 * pkin(2) + pkin(8) * t64 + t61;
t11 = t27 * t84 + (t26 * t46 - t42 * t77) * t50;
t53 = t11 * pkin(3) + t54;
t3 = t11 * t49 - t64 * t83;
t4 = t11 * t83 + t49 * t64;
t52 = t4 * pkin(4) + qJ(5) * t3 + t53;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t70 - mrSges(2,3) - m(3) * t66 - t47 * mrSges(3,3) - (t40 * mrSges(3,1) + t44 * mrSges(3,2)) * t43 - m(4) * t58 - t20 * mrSges(4,1) - t62 * mrSges(4,3) - m(5) * (t57 + t85) - m(6) * (t55 + t85) - m(7) * t55 + t90 * t14 + t89 * t19 + t91 * t15) * g(3) + (-m(5) * (t53 + t87) - m(6) * (t52 + t87) - m(7) * t52 + mrSges(3,3) * t77 - t64 * mrSges(4,3) - m(3) * t61 - m(4) * t54 - t45 * mrSges(2,2) - t41 * mrSges(2,1) - t26 * mrSges(3,2) - t27 * mrSges(3,1) - t11 * mrSges(4,1) - mrSges(1,2) + t92 * r_base(2) + t91 * t4 + t90 * t3 + t89 * t10) * g(2) + (-m(5) * (t59 + t86) - m(6) * (t56 + t86) - m(7) * t56 - m(3) * t67 - m(4) * t60 - t63 * mrSges(4,3) - mrSges(3,3) * t80 - t45 * mrSges(2,1) + t41 * mrSges(2,2) - t29 * mrSges(3,1) - t28 * mrSges(3,2) - t13 * mrSges(4,1) - mrSges(1,1) + t92 * r_base(1) + t91 * t6 + t90 * t5 + t89 * t12) * g(1);
U  = t1;
