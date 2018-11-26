% Calculate Gravitation load on the joints for
% S6PRRRPP2
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
% Datum: 2018-11-23 15:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:21:01
% EndTime: 2018-11-23 15:21:02
% DurationCPUTime: 0.89s
% Computational Cost: add. (1441->117), mult. (1783->154), div. (0->0), fcn. (1793->14), ass. (0->78)
t64 = sin(qJ(4));
t67 = cos(qJ(4));
t132 = pkin(4) * t67 + qJ(5) * t64;
t123 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t83 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t131 = mrSges(4,2) + mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t130 = pkin(3) * t68 + pkin(9) * t65;
t129 = m(6) + m(7);
t128 = mrSges(3,2) - mrSges(4,3);
t106 = cos(pkin(10));
t62 = sin(pkin(10));
t69 = cos(qJ(2));
t114 = t62 * t69;
t103 = pkin(6) + qJ(2);
t90 = sin(t103);
t60 = t90 / 0.2e1;
t104 = pkin(6) - qJ(2);
t91 = sin(t104);
t78 = t60 - t91 / 0.2e1;
t43 = t106 * t78 + t114;
t105 = sin(pkin(6));
t79 = t106 * t105;
t22 = -t43 * t65 - t68 * t79;
t127 = t132 * t22;
t93 = t106 * t69;
t46 = -t62 * t78 + t93;
t94 = t62 * t105;
t24 = -t46 * t65 + t68 * t94;
t126 = t132 * t24;
t88 = cos(t103) / 0.2e1;
t92 = cos(t104);
t58 = t88 - t92 / 0.2e1;
t63 = cos(pkin(6));
t48 = t58 * t65 + t63 * t68;
t125 = t132 * t48;
t124 = t129 * pkin(4) + t83;
t122 = t123 * t64 - t83 * t67 - mrSges(4,1);
t121 = -m(7) * (pkin(9) - qJ(6)) + t131;
t119 = -t68 * mrSges(4,1) - mrSges(3,1) + (m(7) * qJ(6) + t131) * t65;
t113 = t64 * t68;
t109 = t67 * t68;
t66 = sin(qJ(2));
t71 = t92 / 0.2e1 + t88;
t42 = -t106 * t71 + t62 * t66;
t87 = t91 / 0.2e1;
t70 = t87 - t90 / 0.2e1;
t44 = -t106 * t70 + t114;
t100 = -t42 * pkin(2) + pkin(8) * t44;
t45 = t106 * t66 + t62 * t71;
t47 = t62 * t70 + t93;
t99 = -t45 * pkin(2) + pkin(8) * t47;
t57 = t60 + t87;
t98 = t57 * pkin(2) - pkin(8) * t58;
t20 = t22 * pkin(3);
t23 = t43 * t68 - t65 * t79;
t97 = pkin(9) * t23 + t20;
t21 = t24 * pkin(3);
t25 = t46 * t68 + t65 * t94;
t96 = pkin(9) * t25 + t21;
t41 = t48 * pkin(3);
t49 = -t58 * t68 + t63 * t65;
t95 = pkin(9) * t49 + t41;
t82 = -t130 * t42 + t100;
t81 = -t130 * t45 + t99;
t80 = t130 * t57 + t98;
t77 = -t129 * qJ(5) + t123;
t28 = t109 * t57 - t58 * t64;
t27 = t113 * t57 + t58 * t67;
t14 = t49 * t64 + t57 * t67;
t12 = -t109 * t45 + t47 * t64;
t11 = -t113 * t45 - t47 * t67;
t10 = -t109 * t42 + t44 * t64;
t9 = -t113 * t42 - t44 * t67;
t5 = t25 * t64 - t45 * t67;
t3 = t23 * t64 - t42 * t67;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t129) * g(3) (-m(4) * t98 - m(5) * t80 - t129 * (t28 * pkin(4) + qJ(5) * t27 + t80) - t128 * t58 - t83 * t28 + t123 * t27 + t119 * t57) * g(3) + (-m(4) * t100 - m(5) * t82 - t129 * (t10 * pkin(4) + qJ(5) * t9 + t82) + t128 * t44 + t123 * t9 - t83 * t10 - t119 * t42) * g(2) + (-m(4) * t99 - m(5) * t81 - t129 * (t12 * pkin(4) + qJ(5) * t11 + t81) + t128 * t47 - t83 * t12 + t123 * t11 - t119 * t45) * g(1) (-m(5) * t95 - m(6) * (t95 + t125) - m(7) * (t41 + t125) + t121 * t49 + t122 * t48) * g(3) + (-m(5) * t97 - m(6) * (t97 + t127) - m(7) * (t20 + t127) + t121 * t23 + t122 * t22) * g(2) + (-m(5) * t96 - m(6) * (t96 + t126) - m(7) * (t21 + t126) + t121 * t25 + t122 * t24) * g(1) (t77 * (t49 * t67 - t57 * t64) + t124 * t14) * g(3) + (t77 * (t23 * t67 + t42 * t64) + t124 * t3) * g(2) + (t77 * (t25 * t67 + t45 * t64) + t124 * t5) * g(1), t129 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) (-g(1) * t24 - g(2) * t22 - g(3) * t48) * m(7)];
taug  = t1(:);
