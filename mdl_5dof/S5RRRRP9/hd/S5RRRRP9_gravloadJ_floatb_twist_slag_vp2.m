% Calculate Gravitation load on the joints for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:48
% EndTime: 2019-12-31 22:03:51
% DurationCPUTime: 0.82s
% Computational Cost: add. (344->102), mult. (507->122), div. (0->0), fcn. (487->8), ass. (0->56)
t103 = mrSges(5,1) + mrSges(6,1);
t102 = mrSges(5,2) - mrSges(6,3);
t106 = mrSges(3,2) - mrSges(5,3) - mrSges(6,2);
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t90 = -t40 * mrSges(3,1) + t106 * t37;
t105 = t37 * mrSges(4,3) + mrSges(2,1) - t90;
t104 = m(5) + m(6);
t35 = qJ(3) + qJ(4);
t30 = sin(t35);
t31 = cos(t35);
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t100 = m(4) * pkin(2) + t39 * mrSges(4,1) - t36 * mrSges(4,2) - t102 * t30 + t103 * t31;
t98 = -m(6) * qJ(5) - mrSges(6,3);
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t66 = t41 * t36;
t17 = t38 * t39 - t40 * t66;
t29 = t39 * pkin(3) + pkin(2);
t21 = t40 * t29;
t42 = -pkin(8) - pkin(7);
t72 = t37 * t42;
t97 = t104 * (t21 - t72);
t96 = -m(3) - m(4);
t93 = mrSges(2,2) - mrSges(3,3);
t68 = t41 * t30;
t13 = -t38 * t31 + t40 * t68;
t67 = t41 * t31;
t14 = t38 * t30 + t40 * t67;
t92 = t102 * t14 + t103 * t13;
t69 = t38 * t40;
t11 = t30 * t69 + t67;
t12 = t31 * t69 - t68;
t91 = t102 * t12 + t103 * t11;
t88 = m(6) * pkin(4) + t103;
t87 = -mrSges(5,2) - t98;
t82 = pkin(3) * t36;
t79 = g(3) * t37;
t78 = mrSges(5,2) * t31;
t71 = t38 * t36;
t65 = t41 * t39;
t64 = t41 * pkin(1) + t38 * pkin(6);
t61 = m(4) * pkin(7) + mrSges(4,3);
t60 = t98 * t31 * t37;
t58 = -t11 * pkin(4) + t12 * qJ(5);
t56 = -t13 * pkin(4) + t14 * qJ(5);
t55 = t40 * pkin(2) + t37 * pkin(7);
t50 = pkin(4) * t31 + qJ(5) * t30;
t49 = t17 * pkin(3);
t15 = t36 * t69 + t65;
t46 = t15 * pkin(3);
t33 = t41 * pkin(6);
t18 = t40 * t65 + t71;
t16 = -t39 * t69 + t66;
t1 = [(-t18 * mrSges(4,1) - t17 * mrSges(4,2) + t96 * t64 - t104 * (pkin(3) * t71 + t64) + t93 * t38 - t88 * t14 - t87 * t13) * g(2) + (-t16 * mrSges(4,1) - t15 * mrSges(4,2) - t104 * (pkin(3) * t66 + t38 * t72 + t33) + t96 * t33 + t88 * t12 + t87 * t11 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t55) - t104 * (-pkin(1) - t21) + t105) * t38) * g(1) + ((-m(4) * t55 - t105 - t97) * g(2) + t93 * g(1)) * t41, (-t61 * t37 - t97 + (-m(6) * t50 - t100) * t40 + t90) * g(3) + (g(1) * t41 + g(2) * t38) * ((t104 * t42 + t106 - t61) * t40 + (mrSges(3,1) + m(5) * t29 - m(6) * (-t29 - t50) + t100) * t37), -g(3) * ((m(6) * (-pkin(4) * t30 - t82) - t30 * mrSges(6,1)) * t37 - t60) + (m(5) * t82 + mrSges(4,1) * t36 + mrSges(5,1) * t30 + mrSges(4,2) * t39 + t78) * t79 + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + m(5) * t46 - m(6) * (-t46 + t58) + t91) * g(2) + (-t17 * mrSges(4,1) + t18 * mrSges(4,2) - m(5) * t49 - m(6) * (t49 + t56) + t92) * g(1), ((t30 * t88 + t78) * t37 + t60) * g(3) + (-m(6) * t58 + t91) * g(2) + (-m(6) * t56 + t92) * g(1), (-g(1) * t13 - g(2) * t11 - t30 * t79) * m(6)];
taug = t1(:);
