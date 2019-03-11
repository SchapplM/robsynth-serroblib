% Calculate Gravitation load on the joints for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:08
% EndTime: 2019-03-09 04:35:10
% DurationCPUTime: 0.78s
% Computational Cost: add. (434->103), mult. (532->118), div. (0->0), fcn. (510->8), ass. (0->48)
t79 = mrSges(5,1) - mrSges(6,2);
t74 = -mrSges(6,3) - mrSges(7,2);
t72 = mrSges(5,2) + t74;
t78 = mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t77 = -t72 * t27 + t79 * t30;
t76 = m(6) + m(7);
t75 = mrSges(3,2) - mrSges(4,3);
t73 = -m(5) - t76;
t28 = sin(qJ(3));
t31 = cos(qJ(3));
t71 = t31 * mrSges(4,1) + t78 * t28;
t70 = m(7) * qJ(6) + mrSges(7,3);
t69 = pkin(8) * t73;
t44 = m(7) * (-pkin(4) - qJ(6)) - mrSges(7,3);
t53 = qJ(5) * t27;
t47 = -pkin(3) - t53;
t68 = (-m(7) * pkin(5) - mrSges(7,1) - t78) * t31 + (-m(7) * t47 - t44 * t30 - m(6) * (-pkin(4) * t30 + t47) + m(5) * pkin(3) + mrSges(4,1) + t77) * t28;
t67 = -t28 * mrSges(7,1) - t71;
t66 = t70 + t79;
t29 = sin(qJ(1));
t65 = pkin(1) * t29;
t64 = g(3) * t28;
t23 = t28 * pkin(8);
t24 = t31 * pkin(3);
t32 = cos(qJ(1));
t25 = t32 * pkin(1);
t26 = qJ(1) + pkin(9);
t22 = cos(t26);
t63 = t22 * t28;
t62 = t22 * t31;
t60 = t27 * t31;
t56 = t28 * t30;
t55 = t30 * t31;
t54 = t24 + t23;
t21 = sin(t26);
t52 = t22 * pkin(2) + t21 * pkin(7) + t25;
t51 = -pkin(2) - t24;
t50 = t22 * pkin(7) - t65;
t46 = pkin(4) * t55 + t31 * t53 + t54;
t45 = pkin(3) * t62 + pkin(8) * t63 + t52;
t7 = -t21 * t30 + t22 * t60;
t8 = t21 * t27 + t22 * t55;
t36 = t8 * pkin(4) + qJ(5) * t7 + t45;
t6 = t21 * t55 - t22 * t27;
t5 = t21 * t60 + t22 * t30;
t1 = [(-t32 * mrSges(2,1) + t29 * mrSges(2,2) - m(3) * t25 - m(4) * t52 - m(5) * t45 - m(6) * t36 - m(7) * (pkin(5) * t63 + t36) + t75 * t21 - t66 * t8 + t72 * t7 + (-mrSges(3,1) + t67) * t22) * g(2) + (m(3) * t65 + t29 * mrSges(2,1) + t32 * mrSges(2,2) + (-m(4) - m(5)) * t50 - t76 * (-t6 * pkin(4) - qJ(5) * t5 + t50) + t75 * t22 + t66 * t6 - t72 * t5 + (mrSges(3,1) + m(4) * pkin(2) - m(7) * t51 - (m(7) * (-pkin(5) - pkin(8)) - mrSges(7,1)) * t28 + (-m(5) - m(6)) * (t51 - t23) + t71) * t21) * g(1) (-m(3) - m(4) + t73) * g(3) (-m(5) * t54 - m(6) * t46 - m(7) * (pkin(5) * t28 + t46) + (-t70 * t30 - t77) * t31 + t67) * g(3) + (t68 * t22 + t62 * t69) * g(1) + (t31 * t69 + t68) * g(2) * t21 -(-mrSges(5,1) * t27 - mrSges(5,2) * t30) * t64 + ((t74 * t30 + (m(6) * pkin(4) - mrSges(6,2) - t44) * t27) * t28 - t76 * qJ(5) * t56) * g(3) + (-t76 * (-t5 * pkin(4) + qJ(5) * t6) + t72 * t6 + t66 * t5) * g(2) + (-t76 * (-t7 * pkin(4) + qJ(5) * t8) + t72 * t8 + t66 * t7) * g(1), t76 * (-g(1) * t7 - g(2) * t5 - t27 * t64) (-g(1) * t8 - g(2) * t6 - g(3) * t56) * m(7)];
taug  = t1(:);
