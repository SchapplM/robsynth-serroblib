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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:18:43
% EndTime: 2018-11-23 15:18:44
% DurationCPUTime: 1.13s
% Computational Cost: add. (1160->114), mult. (1370->145), div. (0->0), fcn. (1339->16), ass. (0->68)
t126 = m(7) * pkin(5) + mrSges(6,1);
t124 = mrSges(4,1) - mrSges(5,2) + mrSges(6,3) - m(7) * (-pkin(10) - pkin(9)) + mrSges(7,3);
t44 = qJ(5) + qJ(6);
t42 = sin(t44);
t43 = cos(t44);
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t123 = t42 * mrSges(7,1) + t49 * mrSges(6,2) + t43 * mrSges(7,2) + t126 * t46 - mrSges(4,2) + mrSges(5,3);
t112 = m(5) + m(6) + m(7);
t103 = pkin(4) + pkin(8);
t105 = -m(6) * t103 - t49 * mrSges(6,1) + t46 * mrSges(6,2) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t49 + t103) - t43 * mrSges(7,1) + t42 * mrSges(7,2);
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t107 = t123 * t47 + t124 * t50 + mrSges(3,1);
t117 = m(6) * pkin(9) + pkin(3) * t112 + t124;
t45 = sin(pkin(11));
t48 = sin(qJ(2));
t85 = pkin(6) + qJ(2);
t71 = cos(t85) / 0.2e1;
t86 = pkin(6) - qJ(2);
t75 = cos(t86);
t55 = t75 / 0.2e1 + t71;
t88 = cos(pkin(11));
t20 = t45 * t48 - t88 * t55;
t90 = qJ(4) * t47;
t97 = t20 * t50;
t116 = -pkin(3) * t97 - t20 * t90;
t23 = t45 * t55 + t88 * t48;
t96 = t23 * t50;
t115 = -pkin(3) * t96 - t23 * t90;
t73 = sin(t85);
t69 = t73 / 0.2e1;
t74 = sin(t86);
t70 = t74 / 0.2e1;
t33 = t69 + t70;
t95 = t33 * t50;
t114 = pkin(3) * t95 + t33 * t90;
t110 = -t112 * qJ(4) - t123;
t60 = t69 - t74 / 0.2e1;
t51 = cos(qJ(2));
t94 = t45 * t51;
t21 = t88 * t60 + t94;
t87 = sin(pkin(6));
t62 = t88 * t87;
t9 = t21 * t47 + t50 * t62;
t102 = (-t20 * t42 + t43 * t9) * mrSges(7,1) + (-t20 * t43 - t42 * t9) * mrSges(7,2);
t76 = t88 * t51;
t24 = -t45 * t60 + t76;
t77 = t45 * t87;
t11 = t24 * t47 - t50 * t77;
t101 = (t11 * t43 - t23 * t42) * mrSges(7,1) + (-t11 * t42 - t23 * t43) * mrSges(7,2);
t34 = t71 - t75 / 0.2e1;
t89 = cos(pkin(6));
t26 = -t34 * t47 - t89 * t50;
t100 = (t26 * t43 + t33 * t42) * mrSges(7,1) + (-t26 * t42 + t33 * t43) * mrSges(7,2);
t17 = t20 * pkin(2);
t84 = -t17 + t116;
t18 = t23 * pkin(2);
t83 = -t18 + t115;
t32 = t33 * pkin(2);
t82 = t32 + t114;
t61 = t70 - t73 / 0.2e1;
t22 = -t88 * t61 + t94;
t80 = t22 * pkin(8) - t17;
t25 = t45 * t61 + t76;
t79 = t25 * pkin(8) - t18;
t78 = -t34 * pkin(8) + t32;
t1 = [(-m(2) - m(3) - m(4) - t112) * g(3) (-m(4) * t78 - m(5) * (t78 + t114) - m(6) * (pkin(9) * t95 + t82) - m(7) * t82 - t105 * t34 - t107 * t33) * g(3) + (-m(4) * t80 - m(5) * (t80 + t116) - m(6) * (-pkin(9) * t97 + t84) - m(7) * t84 + t105 * t22 + t107 * t20) * g(2) + (-m(4) * t79 - m(5) * (t79 + t115) - m(6) * (-pkin(9) * t96 + t83) - m(7) * t83 + t105 * t25 + t107 * t23) * g(1) (t110 * (-t34 * t50 + t89 * t47) + t117 * t26) * g(3) + (t110 * (t21 * t50 - t47 * t62) + t117 * t9) * g(2) + (t110 * (t24 * t50 + t47 * t77) + t117 * t11) * g(1), t112 * (-g(1) * t11 - g(2) * t9 - g(3) * t26) (-(-t26 * t46 + t33 * t49) * mrSges(6,2) - t100 - t126 * (t26 * t49 + t33 * t46)) * g(3) + (-(-t20 * t49 - t46 * t9) * mrSges(6,2) - t102 - t126 * (-t20 * t46 + t49 * t9)) * g(2) + (-(-t11 * t46 - t23 * t49) * mrSges(6,2) - t101 - t126 * (t11 * t49 - t23 * t46)) * g(1), -g(1) * t101 - g(2) * t102 - g(3) * t100];
taug  = t1(:);
