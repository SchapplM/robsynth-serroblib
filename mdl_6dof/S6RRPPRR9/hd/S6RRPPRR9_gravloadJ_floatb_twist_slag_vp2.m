% Calculate Gravitation load on the joints for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:25
% EndTime: 2019-03-09 09:28:27
% DurationCPUTime: 1.02s
% Computational Cost: add. (426->109), mult. (1007->145), div. (0->0), fcn. (1134->10), ass. (0->57)
t82 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t84 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t94 = m(6) + m(7);
t100 = t94 * (qJ(3) - pkin(9)) + t84;
t90 = m(5) + t94;
t35 = sin(qJ(6));
t39 = cos(qJ(6));
t91 = m(7) * pkin(5) + t39 * mrSges(7,1) - t35 * mrSges(7,2) + mrSges(6,1);
t36 = sin(qJ(5));
t40 = cos(qJ(5));
t89 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t99 = -t91 * t36 - t82 * t40 - t89;
t95 = -m(4) - m(5);
t88 = -mrSges(4,1) - mrSges(5,1) - mrSges(3,3);
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t67 = cos(pkin(6));
t60 = t42 * t67;
t18 = t37 * t38 - t41 * t60;
t19 = t37 * t60 + t38 * t41;
t34 = sin(pkin(6));
t73 = t34 * t42;
t53 = -t19 * t36 + t40 * t73;
t86 = -t18 * t39 + t35 * t53;
t85 = t18 * t35 + t39 * t53;
t54 = -t35 * mrSges(7,1) - t39 * mrSges(7,2);
t81 = -t54 - t100;
t80 = t90 * qJ(4) - t99;
t77 = t34 * t37;
t76 = t34 * t38;
t75 = t34 * t40;
t74 = t34 * t41;
t71 = pkin(2) * t74 + qJ(3) * t77;
t70 = t42 * pkin(1) + pkin(8) * t76;
t69 = qJ(3) * t18;
t61 = t38 * t67;
t20 = t42 * t37 + t41 * t61;
t68 = qJ(3) * t20;
t21 = -t37 * t61 + t41 * t42;
t66 = t21 * pkin(2) + t70;
t62 = -pkin(1) * t38 + pkin(8) * t73;
t57 = -t19 * pkin(2) + t62;
t52 = t19 * t40 + t36 * t73;
t51 = pkin(3) * t76 + qJ(4) * t21 + t66;
t49 = pkin(4) * t76 + t51;
t47 = pkin(3) * t73 - qJ(4) * t19 + t57;
t45 = pkin(4) * t73 + t47;
t17 = t36 * t77 + t67 * t40;
t14 = t20 * pkin(2);
t12 = t18 * pkin(2);
t4 = t21 * t36 + t38 * t75;
t3 = -t21 * t40 + t36 * t76;
t2 = -t20 * t35 + t39 * t4;
t1 = -t20 * t39 - t35 * t4;
t5 = [(-t42 * mrSges(2,1) + t38 * mrSges(2,2) - m(3) * t70 - m(4) * (t66 + t68) - m(5) * (t51 + t68) - m(6) * t49 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t49) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t82 * t3 + t88 * t76 - t89 * t21 - t100 * t20) * g(2) + (t38 * mrSges(2,1) + t42 * mrSges(2,2) - m(3) * t62 - m(4) * (t57 - t69) - m(5) * (t47 - t69) - m(6) * t45 - t53 * mrSges(6,1) - m(7) * (pkin(5) * t53 + t45) - t85 * mrSges(7,1) + t86 * mrSges(7,2) + t82 * t52 + t88 * t73 + t89 * t19 + t100 * t18) * g(1) (-m(4) * t71 - t90 * (qJ(4) * t74 + t71) + (t99 * t41 + (t94 * pkin(9) - t54 - t84) * t37) * t34) * g(3) + (t95 * (qJ(3) * t19 - t12) + t94 * t12 + t81 * t19 + t80 * t18) * g(2) + (t95 * (qJ(3) * t21 - t14) + t94 * t14 + t81 * t21 + t80 * t20) * g(1) (-g(1) * t20 - g(2) * t18 + g(3) * t74) * (m(4) + t90) t90 * (-g(1) * t21 - g(2) * t19 - g(3) * t77) (t82 * t17 - t91 * (-t67 * t36 + t37 * t75)) * g(3) + (-t91 * t52 - t82 * t53) * g(2) + (t91 * t3 + t82 * t4) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t86 + mrSges(7,2) * t85) - g(3) * ((-t17 * t35 + t39 * t74) * mrSges(7,1) + (-t17 * t39 - t35 * t74) * mrSges(7,2))];
taug  = t5(:);
