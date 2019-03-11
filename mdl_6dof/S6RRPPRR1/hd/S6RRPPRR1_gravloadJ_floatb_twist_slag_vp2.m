% Calculate Gravitation load on the joints for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:58
% EndTime: 2019-03-09 08:46:00
% DurationCPUTime: 0.92s
% Computational Cost: add. (470->115), mult. (599->136), div. (0->0), fcn. (591->10), ass. (0->59)
t91 = -m(6) - m(7);
t29 = sin(qJ(6));
t32 = cos(qJ(6));
t90 = m(7) * pkin(5) + t32 * mrSges(7,1) - t29 * mrSges(7,2) + mrSges(6,1);
t27 = qJ(2) + pkin(10);
t24 = sin(t27);
t25 = cos(t27);
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t89 = -t33 * mrSges(3,1) + t30 * mrSges(3,2) - (mrSges(4,1) + mrSges(5,1)) * t25 - (-mrSges(4,2) + mrSges(5,3)) * t24;
t68 = sin(qJ(5));
t69 = cos(qJ(5));
t7 = t24 * t68 + t25 * t69;
t8 = t24 * t69 - t25 * t68;
t87 = -t8 * mrSges(6,2) - t90 * t7;
t34 = cos(qJ(1));
t57 = t34 * t68;
t58 = t34 * t69;
t5 = -t24 * t57 - t25 * t58;
t6 = -t24 * t58 + t25 * t57;
t86 = t5 * mrSges(6,2) - t90 * t6;
t31 = sin(qJ(1));
t3 = t8 * t31;
t4 = t7 * t31;
t85 = -t4 * mrSges(6,2) + t90 * t3;
t84 = g(1) * t34 + g(2) * t31;
t61 = m(7) * pkin(9) + mrSges(7,3);
t83 = -m(3) * pkin(1) - mrSges(2,1) + t89;
t82 = -mrSges(6,2) + t61;
t28 = -qJ(3) - pkin(7);
t81 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3) + t91 * (-pkin(8) - t28);
t66 = qJ(4) * t25;
t13 = t31 * t66;
t73 = pkin(2) * t30;
t74 = -pkin(3) - pkin(4);
t42 = t74 * t24 - t73;
t77 = t42 * t31 + t13;
t20 = t24 * qJ(4);
t67 = t25 * t34;
t76 = pkin(3) * t67 + t34 * t20;
t15 = t34 * t66;
t75 = t42 * t34 + t15;
t22 = t25 * pkin(3);
t26 = t33 * pkin(2);
t65 = m(5) - t91;
t64 = t22 + t20 + t26;
t23 = t26 + pkin(1);
t18 = t34 * t23;
t56 = -t28 * t31 + t18;
t55 = -t23 - t20;
t54 = pkin(4) * t67 + t18 + t76;
t53 = t25 * pkin(4) + t64;
t48 = t4 * t29 - t32 * t34;
t47 = -t29 * t34 - t4 * t32;
t36 = (t74 * t25 + t55) * t31;
t35 = m(5) * (-pkin(3) * t24 - t73) - t24 * mrSges(5,1) + t25 * mrSges(5,3);
t2 = -t29 * t31 - t32 * t5;
t1 = t29 * t5 - t31 * t32;
t9 = [(-m(4) * t56 - m(5) * (t56 + t76) - m(6) * t54 + t5 * mrSges(6,1) - m(7) * (-pkin(5) * t5 + t54) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t82 * t6 + t83 * t34 + t81 * t31) * g(2) + (-m(6) * t36 + t4 * mrSges(6,1) - t47 * mrSges(7,1) - t48 * mrSges(7,2) - (-t4 * pkin(5) + t36) * m(7) - t82 * t3 + (m(4) * t23 - m(5) * (t55 - t22) - t83) * t31 + ((m(4) + m(5)) * t28 + t81) * t34) * g(1) (-m(5) * t13 - t35 * t31 - m(6) * t77 + t4 * mrSges(7,3) - (-pkin(9) * t4 + t77) * m(7) + t85) * g(2) + (-m(5) * t15 - t35 * t34 - m(6) * t75 - t5 * mrSges(7,3) - (t5 * pkin(9) + t75) * m(7) + t86) * g(1) + (-m(4) * t26 - m(5) * t64 - m(6) * t53 - m(7) * (-t8 * pkin(9) + t53) + t8 * mrSges(7,3) + t87 + t89) * g(3) + (m(4) * t73 + mrSges(3,1) * t30 + mrSges(4,1) * t24 + mrSges(3,2) * t33 + mrSges(4,2) * t25) * t84 (-g(1) * t31 + g(2) * t34) * (m(4) + t65) (t25 * g(3) - t84 * t24) * t65 (-t61 * t8 - t87) * g(3) + (-t61 * t4 - t85) * g(2) + (t61 * t5 - t86) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t48 * mrSges(7,1) + t47 * mrSges(7,2)) - g(3) * (-t29 * mrSges(7,1) - t32 * mrSges(7,2)) * t8];
taug  = t9(:);
