% Calculate Gravitation load on the joints for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:05
% EndTime: 2019-03-09 08:06:07
% DurationCPUTime: 0.88s
% Computational Cost: add. (410->104), mult. (564->115), div. (0->0), fcn. (543->10), ass. (0->57)
t96 = mrSges(5,1) + mrSges(6,1);
t84 = mrSges(5,2) - mrSges(6,3);
t85 = m(6) + m(7);
t60 = m(5) + t85;
t95 = m(4) + t60;
t94 = mrSges(6,2) + mrSges(5,3);
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t83 = g(1) * t33 + g(2) * t30;
t25 = sin(pkin(10));
t26 = cos(pkin(10));
t28 = sin(qJ(6));
t31 = cos(qJ(6));
t45 = t25 * t28 + t26 * t31;
t46 = t25 * t31 - t26 * t28;
t93 = t45 * mrSges(7,1) + t46 * mrSges(7,2) - t84 * t25 + t96 * t26;
t92 = m(7) * pkin(8) + mrSges(7,3) - t94;
t29 = sin(qJ(2));
t32 = cos(qJ(2));
t51 = t32 * mrSges(3,1) - t29 * mrSges(3,2);
t89 = -m(3) * pkin(1) - mrSges(2,1) - t51;
t24 = qJ(2) + pkin(9);
t21 = sin(t24);
t16 = t21 * qJ(4);
t22 = cos(t24);
t17 = t22 * pkin(3);
t88 = t17 + t16;
t87 = pkin(4) * t26 + qJ(5) * t25;
t82 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t49 = t22 * mrSges(4,1) - t21 * mrSges(4,2);
t80 = t94 * t21 + t49;
t78 = m(7) * pkin(5) + t96;
t74 = pkin(5) * t26;
t71 = g(3) * t21;
t23 = t32 * pkin(2);
t67 = t26 * t33;
t27 = -qJ(3) - pkin(7);
t66 = t27 * t33;
t65 = t30 * t25;
t64 = t30 * t26;
t63 = t33 * t25;
t59 = t23 + t88;
t20 = t23 + pkin(1);
t56 = -t20 - t17;
t55 = t33 * t20 - t30 * t27;
t54 = t22 * t87 + t59;
t5 = t22 * t65 + t67;
t6 = t22 * t64 - t63;
t53 = t28 * t6 - t31 * t5;
t52 = -t28 * t5 - t31 * t6;
t44 = t33 * t88 + t55;
t43 = -pkin(3) - t87;
t8 = t22 * t67 + t65;
t7 = t22 * t63 - t64;
t2 = t28 * t7 + t31 * t8;
t1 = -t28 * t8 + t31 * t7;
t3 = [(-m(4) * t55 - m(5) * t44 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t85 * (t8 * pkin(4) + t7 * qJ(5) + t44) - t78 * t8 + t84 * t7 + t82 * t30 + (t92 * t21 - t49 + t89) * t33) * g(2) + (m(5) * t66 - t52 * mrSges(7,1) - t53 * mrSges(7,2) - t85 * (-t6 * pkin(4) - t5 * qJ(5) - t66) + t78 * t6 - t84 * t5 + (m(4) * t27 + t82) * t33 + (m(4) * t20 - m(7) * t56 - (m(7) * (pkin(8) - qJ(4)) + mrSges(7,3)) * t21 + (-m(5) - m(6)) * (t56 - t16) + t80 - t89) * t30) * g(1) (-t51 - m(4) * t23 - m(5) * t59 - m(6) * t54 - m(7) * (-pkin(8) * t21 + t54) + t21 * mrSges(7,3) + (-m(7) * t74 - t93) * t22 - t80) * g(3) + (mrSges(3,2) * t32 + (mrSges(4,1) - m(7) * (t43 - t74) - m(6) * t43 + m(5) * pkin(3) + t93) * t21 + (-qJ(4) * t60 + mrSges(4,2) + t92) * t22 + (t95 * pkin(2) + mrSges(3,1)) * t29) * t83 (-g(1) * t30 + g(2) * t33) * t95 (t22 * g(3) - t83 * t21) * t60, t85 * (-g(1) * t7 - g(2) * t5 - t25 * t71) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-mrSges(7,1) * t53 + mrSges(7,2) * t52) - (mrSges(7,1) * t46 - mrSges(7,2) * t45) * t71];
taug  = t3(:);
