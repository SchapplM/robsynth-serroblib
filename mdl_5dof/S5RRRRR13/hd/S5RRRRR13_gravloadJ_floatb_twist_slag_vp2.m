% Calculate Gravitation load on the joints for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:23
% EndTime: 2024-09-27 17:30:24
% DurationCPUTime: 0.26s
% Computational Cost: add. (584->83), mult. (368->98), div. (0->0), fcn. (326->16), ass. (0->54)
t82 = -m(6) * pkin(4) - mrSges(5,1);
t81 = -mrSges(5,3) - mrSges(6,3);
t49 = qJ(4) + qJ(5);
t66 = pkin(5) - t49;
t64 = sin(t66);
t42 = pkin(5) + t49;
t69 = sin(t42) / 0.2e1;
t21 = t69 - t64 / 0.2e1;
t50 = qJ(1) + qJ(2);
t47 = qJ(3) + t50;
t39 = sin(t47);
t40 = cos(t47);
t45 = cos(t49);
t63 = -t40 * t21 - t39 * t45;
t33 = cos(t42) / 0.2e1;
t37 = cos(t66);
t22 = t37 / 0.2e1 + t33;
t43 = sin(t49);
t9 = -t40 * t22 + t39 * t43;
t79 = -t9 * mrSges(6,1) + t63 * mrSges(6,2);
t10 = -t39 * t22 - t40 * t43;
t62 = t39 * t21 - t40 * t45;
t78 = t10 * mrSges(6,1) + t62 * mrSges(6,2);
t46 = cos(t50);
t38 = pkin(2) * t46;
t55 = cos(qJ(4));
t77 = t55 * pkin(4);
t51 = sin(pkin(5));
t76 = t39 * t51;
t52 = cos(pkin(5));
t53 = sin(qJ(4));
t74 = t52 * t53;
t73 = t52 * t55;
t72 = (t69 + t64 / 0.2e1) * mrSges(6,1) + (t33 - t37 / 0.2e1) * mrSges(6,2);
t71 = t40 * pkin(3) + pkin(9) * t76;
t70 = m(4) + m(5) + m(6);
t68 = t38 + t71;
t23 = pkin(4) * t74 - t51 * (pkin(9) + pkin(10));
t41 = pkin(3) + t77;
t67 = -t39 * t23 + t40 * t41;
t65 = t38 + t67;
t61 = -t39 * t53 + t40 * t73;
t17 = -t39 * t73 - t40 * t53;
t18 = -t39 * t74 + t40 * t55;
t60 = -t40 * mrSges(4,1) - t18 * mrSges(5,1) + t62 * mrSges(6,1) + t39 * mrSges(4,2) - t17 * mrSges(5,2) - t10 * mrSges(6,2) + t81 * t76;
t44 = sin(t50);
t59 = -t46 * mrSges(3,1) + t44 * mrSges(3,2) + t60;
t16 = -t39 * t55 - t40 * t74;
t58 = -t16 * mrSges(5,1) - t63 * mrSges(6,1) + t61 * mrSges(5,2) - t9 * mrSges(6,2) + (m(5) * pkin(3) + m(6) * t41 + mrSges(4,1)) * t39 + (m(6) * t23 + mrSges(4,2) + (-m(5) * pkin(9) + t81) * t51) * t40;
t57 = t46 * mrSges(3,2) + (t70 * pkin(2) + mrSges(3,1)) * t44 + t58;
t56 = cos(qJ(1));
t54 = sin(qJ(1));
t48 = t56 * pkin(1);
t1 = [(-m(4) * (t38 + t48) - m(5) * (t48 + t68) - m(6) * (t48 + t65) + t54 * mrSges(2,2) + t59 + (-m(3) * pkin(1) - mrSges(2,1)) * t56) * g(2) + (t56 * mrSges(2,2) + t57 + (mrSges(2,1) + (m(3) + t70) * pkin(1)) * t54) * g(1), (-m(4) * t38 - m(5) * t68 - m(6) * t65 + t59) * g(2) + t57 * g(1), (-m(5) * t71 - m(6) * t67 + t60) * g(2) + t58 * g(1), (-t72 + (-m(6) * t77 - mrSges(5,1) * t55 + mrSges(5,2) * t53) * t51) * g(3) + (-t16 * mrSges(5,2) + t82 * t61 - t79) * g(2) + (t18 * mrSges(5,2) + t82 * t17 - t78) * g(1), -g(1) * t78 - g(2) * t79 - g(3) * t72];
taug = t1(:);
