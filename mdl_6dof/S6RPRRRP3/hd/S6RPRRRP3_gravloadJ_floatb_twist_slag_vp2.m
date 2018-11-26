% Calculate Gravitation load on the joints for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:25:21
% EndTime: 2018-11-23 16:25:22
% DurationCPUTime: 0.78s
% Computational Cost: add. (548->113), mult. (530->133), div. (0->0), fcn. (501->10), ass. (0->61)
t108 = mrSges(6,1) + mrSges(7,1);
t107 = mrSges(6,2) - mrSges(7,3);
t109 = m(6) + m(7);
t106 = mrSges(4,2) - mrSges(7,2) - mrSges(6,3);
t39 = qJ(4) + qJ(5);
t35 = sin(t39);
t36 = cos(t39);
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t105 = m(5) * pkin(3) + t43 * mrSges(5,1) - t40 * mrSges(5,2) - t107 * t35 + t108 * t36;
t38 = qJ(1) + pkin(10);
t33 = sin(t38);
t34 = cos(t38);
t44 = cos(qJ(3));
t75 = t40 * t44;
t17 = t33 * t43 - t34 * t75;
t104 = -m(7) * qJ(6) - mrSges(7,3);
t32 = pkin(4) * t43 + pkin(3);
t27 = t44 * t32;
t41 = sin(qJ(3));
t46 = -pkin(9) - pkin(8);
t71 = t41 * t46;
t103 = t109 * (t27 - t71);
t102 = -m(4) - m(5);
t99 = mrSges(3,2) - mrSges(4,3);
t78 = t35 * t44;
t13 = -t33 * t36 + t34 * t78;
t76 = t36 * t44;
t14 = t33 * t35 + t34 * t76;
t98 = t107 * t14 + t108 * t13;
t11 = t33 * t78 + t34 * t36;
t12 = t33 * t76 - t34 * t35;
t97 = t107 * t12 + t108 * t11;
t96 = -t44 * mrSges(4,1) + t106 * t41;
t94 = m(7) * pkin(5) + t108;
t93 = -mrSges(6,2) - t104;
t92 = t41 * mrSges(5,3) + mrSges(3,1) - t96;
t42 = sin(qJ(1));
t87 = pkin(1) * t42;
t86 = pkin(4) * t40;
t83 = g(3) * t41;
t45 = cos(qJ(1));
t37 = t45 * pkin(1);
t82 = mrSges(6,2) * t36;
t81 = t33 * t40;
t79 = t34 * t40;
t70 = t43 * t44;
t68 = t34 * pkin(2) + t33 * pkin(7) + t37;
t67 = m(5) * pkin(8) + mrSges(5,3);
t66 = t104 * t36 * t41;
t65 = t34 * pkin(7) - t87;
t63 = -t11 * pkin(5) + qJ(6) * t12;
t61 = -t13 * pkin(5) + qJ(6) * t14;
t60 = pkin(3) * t44 + pkin(8) * t41;
t55 = pkin(5) * t36 + qJ(6) * t35;
t54 = t17 * pkin(4);
t15 = t33 * t75 + t34 * t43;
t51 = t15 * pkin(4);
t18 = t34 * t70 + t81;
t16 = -t33 * t70 + t79;
t1 = [(-m(3) * t37 - t45 * mrSges(2,1) - t18 * mrSges(5,1) + t42 * mrSges(2,2) - t17 * mrSges(5,2) + t102 * t68 - t109 * (pkin(4) * t81 + t68) + t99 * t33 - t94 * t14 - t93 * t13 + (-m(5) * t60 - t103 - t92) * t34) * g(2) + (m(3) * t87 + t42 * mrSges(2,1) - t16 * mrSges(5,1) + t45 * mrSges(2,2) - t15 * mrSges(5,2) + t102 * t65 - t109 * (pkin(4) * t79 + t33 * t71 + t65) + t99 * t34 + t94 * t12 + t93 * t11 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t60) - t109 * (-pkin(2) - t27) + t92) * t33) * g(1) (-m(3) - t109 + t102) * g(3) (-t67 * t41 - t103 + (-m(7) * t55 - t105) * t44 + t96) * g(3) + (g(1) * t34 + g(2) * t33) * ((t109 * t46 + t106 - t67) * t44 + (mrSges(4,1) + m(6) * t32 - m(7) * (-t32 - t55) + t105) * t41) -g(3) * ((m(7) * (-pkin(5) * t35 - t86) - t35 * mrSges(7,1)) * t41 - t66) + (m(6) * t86 + mrSges(5,1) * t40 + mrSges(6,1) * t35 + mrSges(5,2) * t43 + t82) * t83 + (t15 * mrSges(5,1) - t16 * mrSges(5,2) + m(6) * t51 - m(7) * (-t51 + t63) + t97) * g(2) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t54 - m(7) * (t54 + t61) + t98) * g(1) ((t94 * t35 + t82) * t41 + t66) * g(3) + (-m(7) * t63 + t97) * g(2) + (-m(7) * t61 + t98) * g(1) (-g(1) * t13 - g(2) * t11 - t35 * t83) * m(7)];
taug  = t1(:);
