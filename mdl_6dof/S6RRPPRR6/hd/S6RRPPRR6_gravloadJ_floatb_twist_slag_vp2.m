% Calculate Gravitation load on the joints for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:06
% EndTime: 2019-03-09 09:13:08
% DurationCPUTime: 1.09s
% Computational Cost: add. (425->134), mult. (670->163), div. (0->0), fcn. (671->10), ass. (0->67)
t40 = sin(qJ(2));
t36 = sin(pkin(10));
t43 = cos(qJ(2));
t84 = t36 * t43;
t37 = cos(pkin(10));
t28 = pkin(4) * t37 + pkin(3);
t86 = -pkin(2) - t28;
t104 = pkin(4) * t84 + t40 * t86;
t39 = sin(qJ(6));
t42 = cos(qJ(6));
t103 = m(7) * pkin(5) + t42 * mrSges(7,1) - t39 * mrSges(7,2) + mrSges(6,1);
t41 = sin(qJ(1));
t75 = qJ(3) * t43;
t24 = t41 * t75;
t102 = t104 * t41 + t24;
t44 = cos(qJ(1));
t25 = t44 * t75;
t101 = t104 * t44 + t25;
t99 = (mrSges(3,1) + mrSges(4,1)) * t43 + (-mrSges(3,2) + mrSges(4,3)) * t40;
t98 = -mrSges(2,1) - t99;
t69 = m(7) * pkin(9) + mrSges(7,3);
t97 = -mrSges(6,2) + t69;
t96 = m(5) * qJ(4) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t34 = t44 * pkin(7);
t38 = -pkin(8) - qJ(4);
t95 = (-pkin(1) + t86 * t43 + (-pkin(4) * t36 - qJ(3)) * t40) * t41 + t44 * t38 + t34;
t73 = pkin(10) + qJ(5);
t63 = sin(t73);
t64 = cos(t73);
t12 = t40 * t64 - t43 * t63;
t94 = g(1) * t44 + g(2) * t41;
t3 = t12 * t41;
t11 = t40 * t63 + t43 * t64;
t4 = t11 * t41;
t92 = -t4 * mrSges(6,2) + t103 * t3;
t49 = t44 * t63;
t50 = t44 * t64;
t5 = -t40 * t49 - t43 * t50;
t6 = -t40 * t50 + t43 * t49;
t91 = t5 * mrSges(6,2) - t103 * t6;
t90 = -t12 * mrSges(6,2) - t103 * t11;
t89 = -pkin(2) - pkin(3);
t33 = t43 * pkin(2);
t85 = t36 * t40;
t81 = t43 * t44;
t31 = t40 * qJ(3);
t77 = t33 + t31;
t76 = t44 * pkin(1) + t41 * pkin(7);
t74 = m(5) + m(6) + m(7);
t23 = pkin(4) * t85;
t71 = t40 * t89;
t67 = -pkin(1) - t31;
t66 = t43 * t28 + t23 + t77;
t65 = pkin(2) * t81 + t44 * t31 + t76;
t60 = -t39 * t44 - t4 * t42;
t59 = t4 * t39 - t42 * t44;
t54 = -t37 * t40 + t84;
t53 = t37 * t43 + t85;
t48 = t44 * t23 + t28 * t81 + t41 * t38 + t65;
t46 = t43 * mrSges(4,3) + (-m(4) * pkin(2) - mrSges(4,1)) * t40;
t10 = t53 * t44;
t9 = t54 * t44;
t8 = t53 * t41;
t7 = t54 * t41;
t2 = -t39 * t41 - t42 * t5;
t1 = t39 * t5 - t41 * t42;
t13 = [(-m(3) * t76 - m(4) * t65 - m(5) * (pkin(3) * t81 + t65) - t10 * mrSges(5,1) + t9 * mrSges(5,2) - m(6) * t48 + t5 * mrSges(6,1) - m(7) * (-pkin(5) * t5 + t48) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t97 * t6 + t98 * t44 + t96 * t41) * g(2) + (t8 * mrSges(5,1) - t7 * mrSges(5,2) - m(6) * t95 + t4 * mrSges(6,1) - t60 * mrSges(7,1) - t59 * mrSges(7,2) - (-pkin(5) * t4 + t95) * m(7) - t97 * t3 + (-m(3) - m(4) - m(5)) * t34 + (m(3) * pkin(1) - m(4) * (t67 - t33) - m(5) * (t89 * t43 + t67) - t98) * t41 + t96 * t44) * g(1), t94 * (mrSges(3,1) * t40 + mrSges(3,2) * t43) + (-m(4) * t24 - t46 * t41 - m(5) * (t41 * t71 + t24) - t7 * mrSges(5,1) - t8 * mrSges(5,2) - m(6) * t102 - m(7) * (-pkin(9) * t4 + t102) + t4 * mrSges(7,3) + t92) * g(2) + (-m(4) * t25 - t46 * t44 - m(5) * (t44 * t71 + t25) - t9 * mrSges(5,1) - t10 * mrSges(5,2) - m(6) * t101 - m(7) * (t5 * pkin(9) + t101) - t5 * mrSges(7,3) + t91) * g(1) + (-m(4) * t77 - m(5) * (pkin(3) * t43 + t77) - t53 * mrSges(5,1) + t54 * mrSges(5,2) - m(6) * t66 - m(7) * (-pkin(9) * t12 + t66) + t12 * mrSges(7,3) + t90 - t99) * g(3) (t43 * g(3) - t40 * t94) * (m(4) + t74) (g(1) * t41 - g(2) * t44) * t74 (-t69 * t12 - t90) * g(3) + (-t69 * t4 - t92) * g(2) + (t69 * t5 - t91) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t59 * mrSges(7,1) + t60 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t39 - mrSges(7,2) * t42) * t12];
taug  = t13(:);
