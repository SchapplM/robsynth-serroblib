% Calculate Gravitation load on the joints for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:27
% EndTime: 2019-03-09 03:44:29
% DurationCPUTime: 0.67s
% Computational Cost: add. (375->102), mult. (404->118), div. (0->0), fcn. (355->10), ass. (0->54)
t88 = mrSges(4,1) - mrSges(5,2);
t87 = -mrSges(4,2) + mrSges(5,3);
t86 = -mrSges(6,3) - mrSges(7,3);
t30 = qJ(5) + qJ(6);
t24 = sin(t30);
t25 = cos(t30);
t31 = sin(qJ(5));
t34 = cos(qJ(5));
t69 = pkin(5) * t31;
t85 = -m(7) * t69 - mrSges(6,1) * t31 - t24 * mrSges(7,1) - t34 * mrSges(6,2) - t25 * mrSges(7,2);
t37 = -pkin(9) - pkin(8);
t84 = -m(6) * (-pkin(3) - pkin(8)) - m(7) * (-pkin(3) + t37) - t86;
t29 = qJ(1) + pkin(10);
t22 = sin(t29);
t23 = cos(t29);
t77 = g(1) * t23 + g(2) * t22;
t52 = m(5) + m(6) + m(7);
t81 = -m(4) - t52;
t79 = -m(7) * pkin(5) - mrSges(6,1);
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t78 = t87 * t32 + t88 * t35;
t75 = t86 * t35 - t78;
t68 = pkin(5) * t34;
t73 = -m(6) * pkin(4) - m(7) * (pkin(4) + t68) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t60 = t25 * t32;
t5 = -t22 * t24 + t23 * t60;
t61 = t24 * t32;
t6 = t22 * t25 + t23 * t61;
t71 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t22 * t60 + t23 * t24;
t8 = -t22 * t61 + t23 * t25;
t70 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t65 = g(3) * t35;
t33 = sin(qJ(1));
t64 = t33 * pkin(1);
t27 = t35 * pkin(3);
t36 = cos(qJ(1));
t28 = t36 * pkin(1);
t63 = mrSges(7,1) * t25;
t62 = t23 * t35;
t59 = t31 * t32;
t58 = t32 * t34;
t55 = t35 * t37;
t26 = t32 * qJ(4);
t54 = t27 + t26;
t51 = t23 * pkin(2) + t22 * pkin(7) + t28;
t49 = -pkin(2) - t26;
t9 = -t22 * t31 + t23 * t58;
t11 = t22 * t58 + t23 * t31;
t17 = t35 * t24 * mrSges(7,2);
t12 = -t22 * t59 + t23 * t34;
t10 = t22 * t34 + t23 * t59;
t1 = [(-m(6) * pkin(8) * t62 - m(3) * t28 - m(4) * t51 - t36 * mrSges(2,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) + t33 * mrSges(2,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) - t52 * (pkin(3) * t62 + t23 * t26 + t51) + t73 * t22 + (-mrSges(3,1) - m(7) * (pkin(5) * t59 - t55) + t75) * t23) * g(2) + (m(3) * t64 + t33 * mrSges(2,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t36 * mrSges(2,2) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + t81 * (t23 * pkin(7) - t64) + t73 * t23 + (mrSges(3,1) + m(4) * pkin(2) - m(5) * (t49 - t27) - m(6) * t49 - m(7) * (-pkin(2) + (-qJ(4) - t69) * t32) + t84 * t35 + t78) * t22) * g(1) (-m(3) + t81) * g(3) (-m(5) * t54 - m(6) * (pkin(8) * t35 + t54) - m(7) * (t54 - t55) + t85 * t32 + t75) * g(3) + ((m(5) * pkin(3) + t84 + t88) * t32 + (-qJ(4) * t52 + t85 - t87) * t35) * t77 (-t77 * t32 + t65) * t52 -(-mrSges(6,1) * t34 + mrSges(6,2) * t31) * t65 - g(3) * (t17 + (-m(7) * t68 - t63) * t35) + (-t12 * mrSges(6,2) + t79 * t11 - t70) * g(2) + (t10 * mrSges(6,2) + t79 * t9 - t71) * g(1), -g(1) * t71 - g(2) * t70 - g(3) * (-t35 * t63 + t17)];
taug  = t1(:);
