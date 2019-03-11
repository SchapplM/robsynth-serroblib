% Calculate Gravitation load on the joints for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:57
% EndTime: 2019-03-09 03:54:58
% DurationCPUTime: 0.48s
% Computational Cost: add. (321->93), mult. (340->108), div. (0->0), fcn. (277->10), ass. (0->55)
t82 = -m(4) - m(5);
t81 = -m(6) - m(7);
t29 = qJ(3) + pkin(10);
t24 = qJ(5) + t29;
t20 = sin(t24);
t21 = cos(t24);
t50 = -pkin(5) * t20 + t21 * pkin(9);
t80 = m(7) * t50;
t36 = cos(qJ(1));
t31 = sin(qJ(6));
t63 = mrSges(7,2) * t31;
t53 = t21 * t63;
t64 = mrSges(6,2) * t20;
t79 = (-t53 - t64) * t36;
t33 = sin(qJ(1));
t34 = cos(qJ(6));
t65 = mrSges(7,1) * t34;
t52 = t21 * t65;
t61 = t21 * t33;
t62 = t20 * t33;
t78 = -mrSges(6,1) * t61 - mrSges(7,3) * t62 - (t52 - t53) * t33;
t77 = mrSges(6,1) * t21 + t20 * mrSges(7,3) + t52;
t76 = -g(1) * t33 + g(2) * t36;
t22 = sin(t29);
t23 = cos(t29);
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t45 = t20 * mrSges(6,1) + t21 * mrSges(6,2);
t71 = pkin(3) * t32;
t75 = -m(5) * t71 - t32 * mrSges(4,1) - t22 * mrSges(5,1) - t35 * mrSges(4,2) - t23 * mrSges(5,2) - t45;
t74 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t17 = t21 * mrSges(7,3);
t73 = mrSges(2,2) - mrSges(3,3) + t17 + t80 + t75;
t70 = pkin(3) * t35;
t10 = pkin(4) * t23 + t70;
t72 = m(6) * t10;
t60 = t31 * t36;
t59 = t33 * t31;
t58 = t33 * t34;
t57 = t34 * t36;
t30 = -qJ(4) - pkin(7);
t56 = pkin(5) * t61 + pkin(9) * t62;
t55 = t36 * pkin(1) + t33 * qJ(2);
t54 = -m(5) + t81;
t26 = t36 * qJ(2);
t51 = -t33 * pkin(1) + t26;
t9 = pkin(4) * t22 + t71;
t48 = -pkin(5) * t21 - pkin(9) * t20;
t42 = t17 + (t63 - t65) * t20;
t28 = -pkin(8) + t30;
t4 = t20 * t57 - t59;
t3 = t20 * t60 + t58;
t2 = t20 * t58 + t60;
t1 = -t20 * t59 + t57;
t5 = [(-t2 * mrSges(7,1) - t1 * mrSges(7,2) + t81 * (-t28 * t36 + t33 * t9 + t55) + (-m(3) + t82) * t55 + (-m(4) * pkin(7) + m(5) * t30 - t74) * t36 + t73 * t33) * g(2) + (-m(3) * t51 - t4 * mrSges(7,1) + t3 * mrSges(7,2) + t81 * (t33 * t28 + t36 * t9 + t51) + t82 * t26 + (-m(4) * (-pkin(1) - pkin(7)) - m(5) * (-pkin(1) + t30) + t74) * t33 + t73 * t36) * g(1), t76 * (m(3) + m(4) - t54) ((t72 - m(7) * (-t10 + t48) + t77) * t36 + t79) * g(2) + (-(-t64 + t72) * t33 - m(7) * (t10 * t33 + t56) + t78) * g(1) + (m(6) * t9 - m(7) * (t50 - t9) - t42 - t75) * g(3) + (m(5) * t70 + mrSges(4,1) * t35 + mrSges(5,1) * t23 - mrSges(4,2) * t32 - mrSges(5,2) * t22) * t76 (g(1) * t36 + g(2) * t33) * t54 (t45 - t42 - t80) * g(3) + ((-m(7) * t48 + t77) * t36 + t79) * g(2) + (-m(7) * t56 + mrSges(6,2) * t62 + t78) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t3 + mrSges(7,2) * t4) - g(3) * (-mrSges(7,1) * t31 - mrSges(7,2) * t34) * t21];
taug  = t5(:);
