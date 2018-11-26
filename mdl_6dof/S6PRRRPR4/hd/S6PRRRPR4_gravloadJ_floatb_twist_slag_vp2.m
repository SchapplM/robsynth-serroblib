% Calculate Gravitation load on the joints for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:24:11
% EndTime: 2018-11-23 15:24:13
% DurationCPUTime: 1.11s
% Computational Cost: add. (1311->118), mult. (1441->152), div. (0->0), fcn. (1417->18), ass. (0->60)
t44 = qJ(4) + pkin(12);
t40 = cos(t44);
t49 = cos(qJ(4));
t42 = t49 * pkin(4);
t28 = pkin(5) * t40 + t42;
t41 = qJ(6) + t44;
t36 = sin(t41);
t37 = cos(t41);
t39 = sin(t44);
t46 = sin(qJ(4));
t102 = mrSges(4,1) + m(7) * (pkin(3) + t28) + t37 * mrSges(7,1) - t36 * mrSges(7,2) + m(6) * (t42 + pkin(3)) + t40 * mrSges(6,1) - t39 * mrSges(6,2) + m(5) * pkin(3) + t49 * mrSges(5,1) - t46 * mrSges(5,2);
t45 = -qJ(5) - pkin(9);
t101 = mrSges(4,2) + m(7) * (-pkin(10) + t45) - mrSges(7,3) + m(6) * t45 - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t108 = m(6) + m(7);
t114 = -m(4) - m(5);
t103 = t108 - t114;
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t110 = pkin(2) * t103 - t101 * t47 + t102 * t50 + mrSges(3,1);
t92 = pkin(4) * t46;
t27 = pkin(5) * t39 + t92;
t109 = m(6) * (pkin(8) + t92) + t39 * mrSges(6,1) + t40 * mrSges(6,2) + m(7) * (pkin(8) + t27) + t36 * mrSges(7,1) + t37 * mrSges(7,2) + t46 * mrSges(5,1) + t49 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3) - t114 * pkin(8);
t106 = -m(6) * pkin(4) - mrSges(5,1);
t48 = sin(qJ(2));
t85 = pkin(6) + qJ(2);
t74 = cos(t85) / 0.2e1;
t86 = pkin(6) - qJ(2);
t78 = cos(t86);
t56 = t78 / 0.2e1 + t74;
t87 = sin(pkin(11));
t89 = cos(pkin(11));
t13 = t87 * t48 - t89 * t56;
t76 = sin(t85);
t72 = t76 / 0.2e1;
t77 = sin(t86);
t63 = t72 - t77 / 0.2e1;
t51 = cos(qJ(2));
t79 = t87 * t51;
t14 = t89 * t63 + t79;
t88 = sin(pkin(6));
t66 = t89 * t88;
t8 = t14 * t50 - t47 * t66;
t95 = (t13 * t37 - t36 * t8) * mrSges(7,1) + (-t13 * t36 - t37 * t8) * mrSges(7,2);
t80 = t89 * t51;
t17 = -t87 * t63 + t80;
t65 = t88 * t87;
t10 = t17 * t50 + t47 * t65;
t16 = t89 * t48 + t87 * t56;
t94 = (-t10 * t36 + t16 * t37) * mrSges(7,1) + (-t10 * t37 - t16 * t36) * mrSges(7,2);
t26 = t74 - t78 / 0.2e1;
t90 = cos(pkin(6));
t20 = -t26 * t50 + t90 * t47;
t73 = t77 / 0.2e1;
t25 = t72 + t73;
t93 = (-t20 * t36 - t25 * t37) * mrSges(7,1) + (-t20 * t37 + t25 * t36) * mrSges(7,2);
t64 = t73 - t76 / 0.2e1;
t19 = -t26 * t47 - t90 * t50;
t9 = t17 * t47 - t50 * t65;
t7 = t14 * t47 + t50 * t66;
t1 = [(-m(2) - m(3) - t103) * g(3) (t109 * t26 - t110 * t25) * g(3) + (-t109 * (-t89 * t64 + t79) + t110 * t13) * g(2) + (-t109 * (t87 * t64 + t80) + t110 * t16) * g(1) (t101 * t20 + t102 * t19) * g(3) + (t101 * t8 + t102 * t7) * g(2) + (t101 * t10 + t102 * t9) * g(1) (-(-t20 * t49 + t25 * t46) * mrSges(5,2) - (-t20 * t39 - t25 * t40) * mrSges(6,1) - (-t20 * t40 + t25 * t39) * mrSges(6,2) - m(7) * (-t20 * t27 - t25 * t28) - t93 + t106 * (-t20 * t46 - t25 * t49)) * g(3) + (-(-t13 * t46 - t49 * t8) * mrSges(5,2) - (t13 * t40 - t39 * t8) * mrSges(6,1) - (-t13 * t39 - t40 * t8) * mrSges(6,2) - m(7) * (t13 * t28 - t27 * t8) - t95 + t106 * (t13 * t49 - t46 * t8)) * g(2) + (-(-t10 * t49 - t16 * t46) * mrSges(5,2) - (-t10 * t39 + t16 * t40) * mrSges(6,1) - (-t10 * t40 - t16 * t39) * mrSges(6,2) - m(7) * (-t10 * t27 + t16 * t28) - t94 + t106 * (-t10 * t46 + t16 * t49)) * g(1), t108 * (-g(1) * t9 - g(2) * t7 - g(3) * t19) -g(1) * t94 - g(2) * t95 - g(3) * t93];
taug  = t1(:);
