% Calculate Gravitation load on the joints for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:19:10
% DurationCPUTime: 0.59s
% Computational Cost: add. (278->97), mult. (348->111), div. (0->0), fcn. (293->8), ass. (0->53)
t31 = sin(qJ(5));
t34 = cos(qJ(5));
t93 = mrSges(6,1) * t31 + mrSges(6,2) * t34;
t30 = qJ(2) + qJ(3);
t28 = cos(t30);
t92 = t93 * t28;
t27 = sin(t30);
t91 = (-mrSges(4,1) + mrSges(5,2)) * t28 + (mrSges(4,2) - mrSges(5,3)) * t27;
t20 = t27 * qJ(4);
t36 = cos(qJ(1));
t69 = t28 * t36;
t89 = pkin(3) * t69 + t36 * t20;
t61 = qJ(4) * t28;
t11 = t36 * t61;
t71 = t27 * t36;
t86 = -m(6) * t11 - mrSges(5,2) * t71 - mrSges(5,3) * t69 - t92 * t36;
t33 = sin(qJ(1));
t72 = t27 * t33;
t9 = t33 * t61;
t85 = -m(6) * t9 - mrSges(5,2) * t72 + (-mrSges(5,3) * t28 - t92) * t33;
t84 = g(1) * t36 + g(2) * t33;
t83 = -t28 * mrSges(6,3) - t93 * t27 + t91;
t32 = sin(qJ(2));
t35 = cos(qJ(2));
t46 = t35 * mrSges(3,1) - t32 * mrSges(3,2);
t82 = -m(3) * pkin(1) - mrSges(2,1) - t46 + t91;
t37 = -pkin(7) - pkin(6);
t81 = -m(3) * pkin(6) - m(6) * (pkin(4) - t37) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t80 = pkin(2) * t32;
t77 = g(3) * t28;
t25 = t28 * pkin(3);
t29 = t35 * pkin(2);
t68 = t33 * t31;
t67 = t33 * t34;
t66 = t36 * t31;
t65 = t36 * t34;
t62 = t25 + t20;
t57 = t29 + t62;
t26 = t29 + pkin(1);
t12 = t36 * t26;
t53 = -t33 * t37 + t12;
t52 = -t26 - t20;
t49 = m(6) * (-pkin(3) - pkin(8)) - mrSges(6,3);
t48 = -pkin(3) * t27 - t80;
t43 = mrSges(4,1) * t27 + mrSges(4,2) * t28;
t41 = t49 * t27;
t38 = -m(6) * t80 + t41;
t24 = t28 * pkin(8);
t4 = -t27 * t68 + t65;
t3 = t27 * t67 + t66;
t2 = t27 * t66 + t67;
t1 = t27 * t65 - t68;
t5 = [(-m(4) * t53 - m(5) * (t53 + t89) - m(6) * (pkin(8) * t69 + t12 + t89) - t2 * mrSges(6,1) - t1 * mrSges(6,2) - mrSges(6,3) * t69 + t82 * t36 + t81 * t33) * g(2) + (-t4 * mrSges(6,1) + t3 * mrSges(6,2) + ((m(4) + m(5)) * t37 + t81) * t36 + (m(4) * t26 - m(5) * (t52 - t25) - m(6) * t52 - t49 * t28 - t82) * t33) * g(1), (-m(5) * (t33 * t48 + t9) - t33 * t38 + t85) * g(2) + (-m(5) * (t36 * t48 + t11) - t36 * t38 + t86) * g(1) + (-t46 - m(4) * t29 - m(5) * t57 - m(6) * (t24 + t57) + t83) * g(3) + (m(4) * t80 + mrSges(3,1) * t32 + mrSges(3,2) * t35 + t43) * t84, t84 * t43 + (-m(5) * (-pkin(3) * t72 + t9) - t33 * t41 + t85) * g(2) + (-m(5) * (-pkin(3) * t71 + t11) - t36 * t41 + t86) * g(1) + (-m(5) * t62 - m(6) * (t24 + t62) + t83) * g(3), (-t27 * t84 + t77) * (m(5) + m(6)), -g(1) * (t1 * mrSges(6,1) - t2 * mrSges(6,2)) - g(2) * (t3 * mrSges(6,1) + t4 * mrSges(6,2)) - (-mrSges(6,1) * t34 + mrSges(6,2) * t31) * t77];
taug = t5(:);
