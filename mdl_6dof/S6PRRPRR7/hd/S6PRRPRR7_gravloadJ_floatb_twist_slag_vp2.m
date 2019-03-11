% Calculate Gravitation load on the joints for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:16
% EndTime: 2019-03-08 22:30:18
% DurationCPUTime: 0.95s
% Computational Cost: add. (503->110), mult. (1151->150), div. (0->0), fcn. (1339->12), ass. (0->59)
t114 = m(7) * pkin(5);
t113 = -mrSges(4,1) + mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3);
t38 = qJ(5) + qJ(6);
t36 = sin(t38);
t37 = cos(t38);
t40 = sin(qJ(5));
t112 = -t36 * mrSges(7,1) - t37 * mrSges(7,2) - t40 * t114 + mrSges(4,2) - mrSges(5,3);
t43 = cos(qJ(5));
t109 = t40 * mrSges(6,1) + t43 * mrSges(6,2);
t99 = m(5) + m(6) + m(7);
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t105 = t112 * t41 + t113 * t44 - mrSges(3,1);
t104 = -t43 * mrSges(6,1) - t37 * mrSges(7,1) + t40 * mrSges(6,2) + t36 * mrSges(7,2) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t97 = -m(6) * pkin(9) - mrSges(6,3);
t103 = pkin(3) * t99 - t113 - t97;
t102 = -mrSges(6,1) - t114;
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t71 = cos(pkin(11));
t72 = cos(pkin(6));
t55 = t72 * t71;
t70 = sin(pkin(11));
t20 = t70 * t42 - t45 * t55;
t73 = qJ(4) * t41;
t84 = t20 * t44;
t101 = -pkin(3) * t84 - t20 * t73;
t54 = t72 * t70;
t22 = t71 * t42 + t45 * t54;
t83 = t22 * t44;
t100 = -pkin(3) * t83 - t22 * t73;
t95 = -t99 * qJ(4) - t109 + t112;
t92 = t44 * mrSges(6,3) + t109 * t41 - t105;
t35 = pkin(5) * t43 + pkin(4);
t91 = -m(6) * (pkin(4) + pkin(8)) - m(7) * (pkin(8) + t35) + t104;
t21 = t42 * t55 + t70 * t45;
t39 = sin(pkin(6));
t64 = t39 * t71;
t9 = t21 * t41 + t44 * t64;
t88 = (-t20 * t36 + t37 * t9) * mrSges(7,1) + (-t20 * t37 - t36 * t9) * mrSges(7,2);
t23 = -t42 * t54 + t71 * t45;
t63 = t39 * t70;
t11 = t23 * t41 - t44 * t63;
t87 = (t11 * t37 - t22 * t36) * mrSges(7,1) + (-t11 * t36 - t22 * t37) * mrSges(7,2);
t82 = t39 * t42;
t24 = t41 * t82 - t72 * t44;
t81 = t39 * t45;
t86 = (t24 * t37 + t36 * t81) * mrSges(7,1) + (-t24 * t36 + t37 * t81) * mrSges(7,2);
t79 = t40 * t45;
t77 = t43 * t45;
t75 = t44 * t45;
t74 = pkin(2) * t81 + pkin(8) * t82;
t17 = t20 * pkin(2);
t69 = -t17 + t101;
t18 = t22 * pkin(2);
t68 = -t18 + t100;
t66 = t21 * pkin(8) - t17;
t65 = t23 * pkin(8) - t18;
t1 = [(-m(2) - m(3) - m(4) - t99) * g(3) (-m(4) * t74 - t99 * (t39 * pkin(3) * t75 + t73 * t81 + t74) + (t97 * t75 + (-t79 * mrSges(6,1) - t77 * mrSges(6,2)) * t41 + t105 * t45 + (-m(6) * pkin(4) - m(7) * t35 + t104) * t42) * t39) * g(3) + (-m(4) * t66 - m(5) * (t66 + t101) - m(6) * (-pkin(9) * t84 + t69) - m(7) * t69 + t91 * t21 + t92 * t20) * g(2) + (-m(4) * t65 - m(5) * (t65 + t100) - m(6) * (-pkin(9) * t83 + t68) - m(7) * t68 + t91 * t23 + t92 * t22) * g(1) (t95 * (t72 * t41 + t44 * t82) + t103 * t24) * g(3) + (t95 * (t21 * t44 - t41 * t64) + t103 * t9) * g(2) + (t95 * (t23 * t44 + t41 * t63) + t103 * t11) * g(1), t99 * (-g(1) * t11 - g(2) * t9 - g(3) * t24) (-(-t24 * t40 + t39 * t77) * mrSges(6,2) - t86 + t102 * (t24 * t43 + t39 * t79)) * g(3) + (-(-t20 * t43 - t40 * t9) * mrSges(6,2) - t88 + t102 * (-t20 * t40 + t43 * t9)) * g(2) + (-(-t11 * t40 - t22 * t43) * mrSges(6,2) - t87 + t102 * (t11 * t43 - t22 * t40)) * g(1), -g(1) * t87 - g(2) * t88 - g(3) * t86];
taug  = t1(:);
