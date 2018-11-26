% Calculate Gravitation load on the joints for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:22:06
% EndTime: 2018-11-23 16:22:07
% DurationCPUTime: 1.42s
% Computational Cost: add. (3665->134), mult. (3888->176), div. (0->0), fcn. (3874->24), ass. (0->87)
t60 = sin(pkin(13));
t61 = cos(pkin(13));
t137 = -m(7) * (pkin(5) * t61 + pkin(4)) - m(6) * pkin(4) - mrSges(6,1) * t61 + mrSges(6,2) * t60 - mrSges(5,1);
t59 = pkin(13) + qJ(6);
t56 = sin(t59);
t57 = cos(t59);
t166 = -mrSges(7,1) * t57 + mrSges(7,2) * t56 + t137;
t156 = -m(5) - m(6);
t139 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(11) - qJ(5)) - mrSges(7,3);
t136 = cos(qJ(1));
t124 = pkin(7) + qJ(3);
t108 = sin(t124) / 0.2e1;
t125 = pkin(7) - qJ(3);
t116 = sin(t125);
t147 = t108 - t116 / 0.2e1;
t135 = sin(qJ(1));
t122 = pkin(6) + pkin(12);
t104 = sin(t122) / 0.2e1;
t123 = pkin(6) - pkin(12);
t114 = sin(t123);
t46 = t104 - t114 / 0.2e1;
t62 = cos(pkin(12));
t40 = t135 * t62 + t136 * t46;
t67 = cos(qJ(3));
t126 = sin(pkin(12));
t105 = cos(t123) / 0.2e1;
t115 = cos(t122);
t89 = t105 + t115 / 0.2e1;
t82 = t135 * t126 - t136 * t89;
t157 = -t147 * t82 + t40 * t67;
t118 = cos(t125);
t110 = t118 / 0.2e1;
t117 = cos(t124);
t103 = t110 - t117 / 0.2e1;
t128 = sin(pkin(6));
t95 = t103 * t128;
t19 = t136 * t95 - t157;
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t129 = cos(pkin(7));
t106 = t129 * t128;
t127 = sin(pkin(7));
t74 = t136 * t106 - t82 * t127;
t165 = t19 * t66 + t74 * t64;
t164 = t19 * t64 - t74 * t66;
t79 = t136 * t126 + t135 * t89;
t68 = t135 * t106 + t79 * t127;
t155 = m(6) + m(7);
t146 = m(5) + t155;
t162 = pkin(3) * t146 - t139 * t64 - t166 * t66 + mrSges(4,1);
t138 = -t60 * mrSges(6,1) - t61 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3) - m(7) * (pkin(5) * t60 + pkin(10));
t160 = -t56 * mrSges(7,1) - t57 * mrSges(7,2) + t138;
t47 = t105 - t115 / 0.2e1;
t88 = t104 + t114 / 0.2e1;
t159 = t147 * t88 + t47 * t67;
t41 = -t135 * t46 + t136 * t62;
t158 = -t147 * t79 + t41 * t67;
t65 = sin(qJ(3));
t91 = t108 + t116 / 0.2e1;
t86 = t91 * t128;
t109 = t117 / 0.2e1;
t94 = t110 + t109;
t15 = t136 * t86 + t40 * t65 + t82 * t94;
t140 = t156 * pkin(10) + t160;
t111 = t128 * t135;
t131 = t136 * pkin(1) + qJ(2) * t111;
t130 = cos(pkin(6));
t112 = t136 * t128;
t113 = -t135 * pkin(1) + qJ(2) * t112;
t93 = t109 - t118 / 0.2e1;
t87 = t93 * t128;
t78 = -t88 * t127 + t130 * t129;
t76 = -t40 * pkin(2) + t74 * pkin(9) + t113;
t73 = t19 * pkin(3) + t76;
t72 = t41 * pkin(2) + t68 * pkin(9) + t131;
t21 = t135 * t95 + t158;
t71 = t21 * pkin(3) + t72;
t25 = t130 * t103 + t159;
t24 = -t130 * t91 + t47 * t65 - t88 * t94;
t20 = -t135 * t86 + t41 * t65 + t79 * t94;
t10 = t25 * t66 + t78 * t64;
t9 = t25 * t64 - t78 * t66;
t8 = t21 * t66 + t68 * t64;
t7 = t21 * t64 - t68 * t66;
t2 = t20 * t56 + t57 * t8;
t1 = t20 * t57 - t56 * t8;
t3 = [(-m(3) * t131 - m(4) * t72 - m(7) * t71 - t136 * mrSges(2,1) - t41 * mrSges(3,1) - t21 * mrSges(4,1) - t2 * mrSges(7,1) + t135 * mrSges(2,2) + t79 * mrSges(3,2) - t1 * mrSges(7,2) - mrSges(3,3) * t111 - t68 * mrSges(4,3) + t156 * (t20 * pkin(10) + t71) + t137 * t8 + t138 * t20 + t139 * t7) * g(2) + (-m(3) * t113 - m(4) * t76 - m(7) * t73 + t135 * mrSges(2,1) + t40 * mrSges(3,1) - t19 * mrSges(4,1) + t136 * mrSges(2,2) - t82 * mrSges(3,2) - mrSges(3,3) * t112 - t74 * mrSges(4,3) + t156 * (-pkin(10) * t15 + t73) + t166 * t165 - t160 * t15 + t139 * t164) * g(1) (-g(1) * t111 + g(2) * t112 - g(3) * t130) * (m(3) + m(4) + t146) (t140 * (-t130 * t93 + t159) + t162 * t24) * g(3) + (t140 * (t136 * t87 + t157) + t162 * t15) * g(2) + (t140 * (-t135 * t87 + t158) + t162 * t20) * g(1) (t139 * t10 - t166 * t9) * g(3) + (-t139 * t165 + t164 * t166) * g(2) + (t139 * t8 - t166 * t7) * g(1), t155 * (-g(1) * t7 + g(2) * t164 - g(3) * t9) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t15 * t57 + t165 * t56) * mrSges(7,1) + (-t15 * t56 + t165 * t57) * mrSges(7,2)) - g(3) * ((-t10 * t56 + t24 * t57) * mrSges(7,1) + (-t10 * t57 - t24 * t56) * mrSges(7,2))];
taug  = t3(:);
