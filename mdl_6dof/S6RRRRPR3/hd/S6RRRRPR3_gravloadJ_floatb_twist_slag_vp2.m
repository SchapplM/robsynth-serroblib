% Calculate Gravitation load on the joints for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:13:05
% EndTime: 2018-11-23 18:13:06
% DurationCPUTime: 0.76s
% Computational Cost: add. (592->133), mult. (516->136), div. (0->0), fcn. (430->10), ass. (0->75)
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t123 = mrSges(7,1) * t41 + mrSges(7,2) * t44;
t40 = qJ(2) + qJ(3);
t35 = sin(t40);
t36 = cos(t40);
t37 = qJ(4) + t40;
t32 = sin(t37);
t33 = cos(t37);
t58 = mrSges(5,1) * t32 + mrSges(5,2) * t33;
t122 = mrSges(4,1) * t35 + mrSges(4,2) * t36 + t58;
t121 = t123 * t33;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t112 = g(1) * t46 + g(2) * t43;
t120 = (-mrSges(5,1) + mrSges(6,2)) * t33 + (mrSges(5,2) - mrSges(6,3)) * t32;
t24 = t32 * qJ(5);
t91 = t33 * t46;
t118 = pkin(4) * t91 + t46 * t24;
t72 = t36 * mrSges(4,1) - t35 * mrSges(4,2);
t80 = qJ(5) * t33;
t16 = t46 * t80;
t93 = t32 * t46;
t83 = mrSges(6,2) * t93 + mrSges(6,3) * t91;
t114 = -m(7) * t16 - t83;
t14 = t43 * t80;
t94 = t32 * t43;
t84 = t43 * t33 * mrSges(6,3) + mrSges(6,2) * t94;
t113 = -m(7) * t14 - t84;
t110 = -t33 * mrSges(7,3) - t123 * t32 + t120;
t109 = -t72 + t110;
t45 = cos(qJ(2));
t38 = t45 * pkin(2);
t42 = sin(qJ(2));
t62 = t45 * mrSges(3,1) - t42 * mrSges(3,2);
t108 = -mrSges(2,1) - m(3) * pkin(1) - t62 - m(4) * (t38 + pkin(1)) - t72 + t120;
t47 = -pkin(8) - pkin(7);
t39 = -pkin(9) + t47;
t107 = -m(3) * pkin(7) + m(4) * t47 - m(7) * (pkin(5) - t39) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t106 = t121 * t43;
t103 = pkin(2) * t42;
t102 = pkin(3) * t35;
t31 = pkin(3) * t36;
t99 = g(3) * t33;
t29 = t33 * pkin(4);
t89 = t41 * t43;
t88 = t41 * t46;
t87 = t43 * t44;
t86 = t44 * t46;
t85 = t121 * t46;
t82 = t29 + t24;
t81 = t31 + t38;
t76 = t31 + t82;
t12 = pkin(1) + t81;
t5 = t46 * t12;
t73 = -t39 * t43 + t5;
t69 = -t12 - t24;
t68 = t38 + t76;
t67 = -pkin(4) * t94 + t14;
t66 = -pkin(4) * t93 + t16;
t65 = m(7) * (-pkin(4) - pkin(10)) - mrSges(7,3);
t64 = -pkin(4) * t32 - t102;
t56 = t65 * t32;
t50 = t43 * t56 + t106;
t49 = t46 * t56 + t85;
t48 = -m(7) * t102 + t56;
t28 = t33 * pkin(10);
t13 = -t102 - t103;
t7 = t46 * t13;
t6 = t43 * t13;
t4 = -t32 * t89 + t86;
t3 = t32 * t87 + t88;
t2 = t32 * t88 + t87;
t1 = t32 * t86 - t89;
t8 = [(-m(5) * t73 - m(6) * (t73 + t118) - m(7) * (pkin(10) * t91 + t118 + t5) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - mrSges(7,3) * t91 + t108 * t46 + t107 * t43) * g(2) + (-t4 * mrSges(7,1) + t3 * mrSges(7,2) + ((m(5) + m(6)) * t39 + t107) * t46 + (m(5) * t12 - m(6) * (t69 - t29) - m(7) * t69 - t65 * t33 - t108) * t43) * g(1) (-m(6) * (t6 + t67) - t84 - m(7) * (t14 + t6) - t50) * g(2) + (-m(6) * (t66 + t7) - t83 - m(7) * (t16 + t7) - t49) * g(1) + (-t62 - m(4) * t38 - m(5) * t81 - m(6) * t68 - m(7) * (t28 + t68) + t109) * g(3) + t112 * (m(4) * t103 - m(5) * t13 + mrSges(3,1) * t42 + mrSges(3,2) * t45 + t122) (-m(6) * (t43 * t64 + t14) - t48 * t43 - t106 + t113) * g(2) + (-m(6) * (t46 * t64 + t16) - t48 * t46 - t85 + t114) * g(1) + (-m(5) * t31 - m(6) * t76 - m(7) * (t28 + t76) + t109) * g(3) + t112 * (m(5) * t102 + t122) t112 * t58 + (-m(6) * t67 + t113 - t50) * g(2) + (-m(6) * t66 + t114 - t49) * g(1) + (-m(6) * t82 - m(7) * (t28 + t82) + t110) * g(3) (-t112 * t32 + t99) * (m(6) + m(7)) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t3 + mrSges(7,2) * t4) - (-mrSges(7,1) * t44 + mrSges(7,2) * t41) * t99];
taug  = t8(:);
