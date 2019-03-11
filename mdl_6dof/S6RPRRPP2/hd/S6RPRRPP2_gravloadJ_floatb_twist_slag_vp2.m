% Calculate Gravitation load on the joints for
% S6RPRRPP2
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:30
% EndTime: 2019-03-09 04:31:31
% DurationCPUTime: 0.82s
% Computational Cost: add. (430->101), mult. (527->113), div. (0->0), fcn. (503->8), ass. (0->47)
t82 = mrSges(5,1) + mrSges(6,1);
t73 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t81 = mrSges(6,2) + mrSges(5,3);
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t80 = -t73 * t27 + t82 * t30;
t26 = qJ(1) + pkin(9);
t21 = sin(t26);
t79 = g(2) * t21;
t28 = sin(qJ(3));
t63 = g(3) * t28;
t78 = m(6) + m(7);
t77 = mrSges(3,2) - mrSges(4,3);
t22 = cos(t26);
t75 = g(1) * t22 + t79;
t74 = -m(5) - t78;
t31 = cos(qJ(3));
t45 = m(7) * (-pkin(4) - pkin(5)) - mrSges(7,1);
t54 = qJ(5) * t27;
t47 = -pkin(3) - t54;
t72 = (m(7) * qJ(6) + mrSges(7,3) - t81) * t31 + (-m(7) * t47 - t45 * t30 - m(6) * (-pkin(4) * t30 + t47) + m(5) * pkin(3) + t80) * t28;
t71 = t31 * mrSges(4,1) + (-mrSges(4,2) + t81) * t28;
t70 = m(7) * pkin(5) + mrSges(7,1);
t69 = pkin(8) * t74;
t68 = t28 * mrSges(7,3) - t71;
t67 = t70 + t82;
t29 = sin(qJ(1));
t66 = pkin(1) * t29;
t23 = t28 * pkin(8);
t24 = t31 * pkin(3);
t32 = cos(qJ(1));
t25 = t32 * pkin(1);
t62 = t22 * t31;
t60 = t27 * t31;
t56 = t30 * t31;
t55 = t24 + t23;
t53 = qJ(6) * t28;
t52 = t22 * pkin(2) + t21 * pkin(7) + t25;
t51 = -pkin(2) - t24;
t50 = t22 * pkin(7) - t66;
t46 = pkin(4) * t56 + t31 * t54 + t55;
t44 = pkin(3) * t62 + t22 * t23 + t52;
t8 = t21 * t27 + t22 * t56;
t7 = -t21 * t30 + t22 * t60;
t6 = t21 * t56 - t22 * t27;
t5 = t21 * t60 + t22 * t30;
t1 = [(-m(3) * t25 - m(4) * t52 - m(5) * t44 - t32 * mrSges(2,1) + t29 * mrSges(2,2) - t78 * (t8 * pkin(4) + qJ(5) * t7 + t44) + t77 * t21 - t67 * t8 + t73 * t7 + (m(7) * t53 - mrSges(3,1) + t68) * t22) * g(2) + (m(3) * t66 + t29 * mrSges(2,1) + t32 * mrSges(2,2) + (-m(4) - m(5)) * t50 - t78 * (-t6 * pkin(4) - qJ(5) * t5 + t50) + t77 * t22 + t67 * t6 - t73 * t5 + (mrSges(3,1) + m(4) * pkin(2) - m(7) * t51 - (m(7) * (-pkin(8) + qJ(6)) + mrSges(7,3)) * t28 + (-m(5) - m(6)) * (t51 - t23) + t71) * t21) * g(1) (-m(3) - m(4) + t74) * g(3), t75 * (mrSges(4,1) * t28 + mrSges(4,2) * t31) + (t72 * t22 + t62 * t69) * g(1) + (-m(5) * t55 - m(6) * t46 - m(7) * (t46 - t53) + (-t70 * t30 - t80) * t31 + t68) * g(3) + (t31 * t69 + t72) * t79 (-t78 * (-t5 * pkin(4) + qJ(5) * t6) + t73 * t6 + t67 * t5) * g(2) + (-t78 * (-t7 * pkin(4) + qJ(5) * t8) + t73 * t8 + t67 * t7) * g(1) + ((-t78 * qJ(5) + t73) * t30 + (m(6) * pkin(4) - t45 + t82) * t27) * t63, t78 * (-g(1) * t7 - g(2) * t5 - t27 * t63) (-g(3) * t31 + t75 * t28) * m(7)];
taug  = t1(:);
