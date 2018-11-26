% Calculate Gravitation load on the joints for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2018-11-23 15:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:21:48
% EndTime: 2018-11-23 15:21:49
% DurationCPUTime: 0.88s
% Computational Cost: add. (1455->124), mult. (1799->164), div. (0->0), fcn. (1811->14), ass. (0->87)
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t133 = pkin(4) * t68 + qJ(5) * t65;
t79 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t125 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t132 = mrSges(6,1) + mrSges(7,1) - mrSges(4,2) + mrSges(5,3);
t131 = m(6) + m(7);
t130 = mrSges(3,2) - mrSges(4,3);
t106 = cos(pkin(10));
t103 = pkin(6) + qJ(2);
t90 = sin(t103);
t60 = t90 / 0.2e1;
t104 = pkin(6) - qJ(2);
t91 = sin(t104);
t108 = t60 - t91 / 0.2e1;
t63 = sin(pkin(10));
t70 = cos(qJ(2));
t114 = t63 * t70;
t43 = t106 * t108 + t114;
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t105 = sin(pkin(6));
t80 = t106 * t105;
t22 = -t43 * t66 - t69 * t80;
t129 = t133 * t22;
t93 = t106 * t70;
t46 = -t63 * t108 + t93;
t94 = t63 * t105;
t24 = -t46 * t66 + t69 * t94;
t128 = t133 * t24;
t88 = cos(t103) / 0.2e1;
t92 = cos(t104);
t58 = t88 - t92 / 0.2e1;
t64 = cos(pkin(6));
t48 = t58 * t66 + t64 * t69;
t127 = t133 * t48;
t126 = pkin(4) * t131 + t79;
t124 = t125 * t65 - t79 * t68 - mrSges(4,1);
t123 = -m(7) * (pkin(5) + pkin(9)) - t132;
t122 = t69 * mrSges(4,1) + t132 * t66 + mrSges(3,1);
t119 = pkin(3) * t69;
t67 = sin(qJ(2));
t72 = t92 / 0.2e1 + t88;
t42 = -t106 * t72 + t63 * t67;
t117 = t42 * t66;
t45 = t106 * t67 + t63 * t72;
t116 = t45 * t66;
t87 = t91 / 0.2e1;
t57 = t60 + t87;
t115 = t57 * t66;
t113 = t65 * t69;
t109 = t68 * t69;
t71 = t87 - t90 / 0.2e1;
t44 = -t106 * t71 + t114;
t101 = -t42 * pkin(2) + pkin(8) * t44;
t47 = t63 * t71 + t93;
t100 = -t45 * pkin(2) + pkin(8) * t47;
t99 = t57 * pkin(2) - pkin(8) * t58;
t20 = t22 * pkin(3);
t23 = t43 * t69 - t66 * t80;
t98 = pkin(9) * t23 + t20;
t21 = t24 * pkin(3);
t25 = t46 * t69 + t66 * t94;
t97 = pkin(9) * t25 + t21;
t41 = t48 * pkin(3);
t49 = -t58 * t69 + t64 * t66;
t96 = pkin(9) * t49 + t41;
t83 = -pkin(9) * t117 - t42 * t119 + t101;
t82 = -pkin(9) * t116 - t45 * t119 + t100;
t81 = pkin(9) * t115 + t57 * t119 + t99;
t78 = -qJ(5) * t131 + t125;
t10 = -t42 * t109 + t44 * t65;
t9 = -t42 * t113 - t44 * t68;
t75 = t10 * pkin(4) + qJ(5) * t9 + t83;
t11 = -t45 * t113 - t47 * t68;
t12 = -t45 * t109 + t47 * t65;
t74 = t12 * pkin(4) + qJ(5) * t11 + t82;
t27 = t57 * t113 + t58 * t68;
t28 = t57 * t109 - t58 * t65;
t73 = t28 * pkin(4) + qJ(5) * t27 + t81;
t15 = t49 * t68 - t57 * t65;
t14 = t49 * t65 + t57 * t68;
t6 = t25 * t68 + t45 * t65;
t5 = t25 * t65 - t45 * t68;
t4 = t23 * t68 + t42 * t65;
t3 = t23 * t65 - t42 * t68;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t131) * g(3) (-m(4) * t99 - m(5) * t81 - m(6) * t73 - m(7) * (pkin(5) * t115 + t73) - t130 * t58 - t79 * t28 + t125 * t27 - t122 * t57) * g(3) + (-m(4) * t101 - m(5) * t83 - m(6) * t75 - m(7) * (-pkin(5) * t117 + t75) + t130 * t44 + t125 * t9 - t79 * t10 + t122 * t42) * g(2) + (-m(4) * t100 - m(5) * t82 - m(6) * t74 - m(7) * (-pkin(5) * t116 + t74) + t130 * t47 - t79 * t12 + t125 * t11 + t122 * t45) * g(1) (-m(5) * t96 - m(6) * (t96 + t127) - m(7) * (t41 + t127) + t123 * t49 + t124 * t48) * g(3) + (-m(5) * t98 - m(6) * (t98 + t129) - m(7) * (t20 + t129) + t123 * t23 + t124 * t22) * g(2) + (-m(5) * t97 - m(6) * (t97 + t128) - m(7) * (t21 + t128) + t123 * t25 + t124 * t24) * g(1) (t126 * t14 + t78 * t15) * g(3) + (t126 * t3 + t78 * t4) * g(2) + (t126 * t5 + t78 * t6) * g(1), t131 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) (-g(1) * t6 - g(2) * t4 - g(3) * t15) * m(7)];
taug  = t1(:);
