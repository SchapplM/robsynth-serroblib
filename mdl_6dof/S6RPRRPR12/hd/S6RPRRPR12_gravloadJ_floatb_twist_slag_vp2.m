% Calculate Gravitation load on the joints for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:09
% EndTime: 2019-03-09 05:48:13
% DurationCPUTime: 1.36s
% Computational Cost: add. (1098->137), mult. (2969->186), div. (0->0), fcn. (3757->14), ass. (0->78)
t105 = cos(pkin(7));
t117 = cos(qJ(1));
t101 = sin(pkin(12));
t115 = sin(qJ(1));
t104 = cos(pkin(12));
t106 = cos(pkin(6));
t91 = t106 * t104;
t78 = t115 * t101 - t117 * t91;
t102 = sin(pkin(7));
t103 = sin(pkin(6));
t87 = t103 * t102;
t146 = t78 * t105 + t117 * t87;
t135 = mrSges(5,2) - mrSges(6,3);
t56 = sin(qJ(6));
t59 = cos(qJ(6));
t145 = t56 * mrSges(7,1) + t59 * mrSges(7,2) - t135;
t116 = cos(qJ(3));
t89 = t106 * t101;
t44 = t104 * t115 + t117 * t89;
t58 = sin(qJ(3));
t27 = -t116 * t44 + t146 * t58;
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t88 = t105 * t103;
t67 = -t78 * t102 + t117 * t88;
t9 = t27 * t57 - t60 * t67;
t10 = t27 * t60 + t57 * t67;
t136 = m(6) + m(7);
t141 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t73 = t101 * t117 + t115 * t91;
t62 = t73 * t102 + t115 * t88;
t125 = m(7) * pkin(11) + t141;
t137 = pkin(4) * t136 + t125;
t24 = t146 * t116 + t44 * t58;
t107 = qJ(5) * t57;
t114 = t24 * t60;
t134 = -pkin(4) * t114 - t24 * t107;
t131 = t73 * t105 - t115 * t87;
t45 = t104 * t117 - t115 * t89;
t28 = t131 * t116 + t45 * t58;
t113 = t28 * t60;
t133 = -pkin(4) * t113 - t28 * t107;
t130 = t106 * t102 + t104 * t88;
t86 = t103 * t101;
t36 = -t130 * t116 + t58 * t86;
t112 = t36 * t60;
t132 = -pkin(4) * t112 - t36 * t107;
t124 = -t136 * qJ(5) - t145;
t123 = -m(7) * (pkin(5) + pkin(10)) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t122 = t141 * t60 + t145 * t57 + mrSges(4,1);
t121 = -t59 * mrSges(7,1) + t56 * mrSges(7,2) + t123;
t119 = t24 * pkin(10);
t118 = t28 * pkin(10);
t94 = t103 * t115;
t108 = t117 * pkin(1) + qJ(2) * t94;
t18 = t24 * pkin(3);
t100 = -pkin(10) * t27 - t18;
t20 = t28 * pkin(3);
t29 = t45 * t116 - t131 * t58;
t99 = t29 * pkin(10) - t20;
t35 = t36 * pkin(3);
t37 = t116 * t86 + t130 * t58;
t98 = t37 * pkin(10) - t35;
t95 = t117 * t103;
t97 = -pkin(1) * t115 + qJ(2) * t95;
t72 = -t104 * t87 + t105 * t106;
t69 = -t44 * pkin(2) + t67 * pkin(9) + t97;
t66 = t27 * pkin(3) + t69;
t65 = t45 * pkin(2) + t62 * pkin(9) + t108;
t64 = t29 * pkin(3) + t65;
t63 = t10 * pkin(4) + t9 * qJ(5) + t66;
t11 = t29 * t57 - t60 * t62;
t12 = t29 * t60 + t57 * t62;
t61 = t12 * pkin(4) + t11 * qJ(5) + t64;
t22 = t37 * t57 - t60 * t72;
t2 = t11 * t56 + t28 * t59;
t1 = t11 * t59 - t28 * t56;
t3 = [(-mrSges(2,1) * t117 + mrSges(2,2) * t115 - m(3) * t108 - t45 * mrSges(3,1) + mrSges(3,2) * t73 - mrSges(3,3) * t94 - m(4) * t65 - t29 * mrSges(4,1) - mrSges(4,3) * t62 - m(5) * (t64 + t118) - m(6) * (t61 + t118) - m(7) * t61 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t135 * t11 + t123 * t28 - t125 * t12) * g(2) + (mrSges(2,1) * t115 + mrSges(2,2) * t117 - m(3) * t97 + t44 * mrSges(3,1) - mrSges(3,2) * t78 - mrSges(3,3) * t95 - m(4) * t69 - t27 * mrSges(4,1) - mrSges(4,3) * t67 - m(5) * (t66 - t119) - m(6) * (t63 - t119) - m(7) * t63 - t145 * t9 - t121 * t24 - t125 * t10) * g(1) (-g(1) * t94 + g(2) * t95 - g(3) * t106) * (m(3) + m(4) + m(5) + t136) (-m(5) * t98 - m(6) * (t98 + t132) - m(7) * (-pkin(11) * t112 + t132 - t35) + t121 * t37 + t122 * t36) * g(3) + (-m(5) * t100 - m(6) * (t100 + t134) - m(7) * (-pkin(11) * t114 + t134 - t18) - t121 * t27 + t122 * t24) * g(2) + (-m(5) * t99 - m(6) * (t99 + t133) - m(7) * (-pkin(11) * t113 + t133 - t20) + t121 * t29 + t122 * t28) * g(1) (t124 * (t37 * t60 + t57 * t72) + t137 * t22) * g(3) + (-t124 * t10 - t137 * t9) * g(2) + (t137 * t11 + t124 * t12) * g(1), t136 * (-g(1) * t11 + g(2) * t9 - g(3) * t22) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t24 * t56 - t59 * t9) * mrSges(7,1) + (-t24 * t59 + t56 * t9) * mrSges(7,2)) - g(3) * ((t22 * t59 - t36 * t56) * mrSges(7,1) + (-t22 * t56 - t36 * t59) * mrSges(7,2))];
taug  = t3(:);
