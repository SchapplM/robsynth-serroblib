% Calculate Gravitation load on the joints for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:52
% EndTime: 2019-03-08 22:42:55
% DurationCPUTime: 1.16s
% Computational Cost: add. (631->108), mult. (1326->153), div. (0->0), fcn. (1561->12), ass. (0->61)
t111 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t127 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t55 = qJ(4) + pkin(11);
t53 = sin(t55);
t54 = cos(t55);
t132 = t111 * t54 + t127 * t53;
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t131 = -m(5) * pkin(3) - t61 * mrSges(5,1) + t58 * mrSges(5,2) - mrSges(4,1);
t130 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(5,3);
t118 = m(6) + m(7);
t52 = pkin(4) * t61 + pkin(3);
t128 = t118 * t52;
t115 = -mrSges(6,3) - mrSges(7,2);
t126 = -m(4) - t118;
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t125 = t130 * t59 + t131 * t62 - mrSges(3,1);
t124 = -t58 * mrSges(5,1) - t61 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t60 = sin(qJ(2));
t63 = cos(qJ(2));
t88 = sin(pkin(10));
t90 = cos(pkin(6));
t74 = t90 * t88;
t89 = cos(pkin(10));
t38 = -t60 * t74 + t63 * t89;
t56 = sin(pkin(6));
t81 = t56 * t88;
t18 = t38 * t62 + t59 * t81;
t37 = t60 * t89 + t63 * t74;
t123 = -t18 * t58 + t37 * t61;
t75 = t90 * t89;
t36 = t60 * t75 + t63 * t88;
t82 = t56 * t89;
t16 = t36 * t62 - t59 * t82;
t35 = t60 * t88 - t63 * t75;
t122 = -t16 * t58 + t35 * t61;
t108 = pkin(4) * t58;
t121 = -m(5) * pkin(8) - t118 * t108 - t111 * t53 + t127 * t54 + t124;
t57 = -qJ(5) - pkin(9);
t120 = -t125 + (-t118 * t57 - t115) * t59 + (t128 + t132) * t62;
t119 = -m(4) - m(5);
t113 = -t131 + t132;
t112 = t115 + t130;
t99 = t56 * t60;
t98 = t56 * t63;
t94 = t62 * t63;
t91 = pkin(2) * t98 + pkin(8) * t99;
t87 = t59 * t98;
t86 = t53 * t98;
t40 = t59 * t90 + t62 * t99;
t73 = -t40 * t58 - t61 * t98;
t39 = t59 * t99 - t62 * t90;
t34 = t37 * pkin(2);
t33 = t35 * pkin(2);
t17 = t38 * t59 - t62 * t81;
t15 = t36 * t59 + t62 * t82;
t11 = t40 * t53 + t54 * t98;
t3 = t18 * t53 - t37 * t54;
t1 = t16 * t53 - t35 * t54;
t2 = [(-m(2) - m(3) - t118 + t119) * g(3) (t119 * t91 + t115 * t87 - t118 * (t99 * t108 - t57 * t87 + t91) - t127 * (-t54 * t99 + t62 * t86) + (-t94 * t128 - t111 * (t53 * t60 + t54 * t94) + t125 * t63 + t124 * t60) * t56) * g(3) + (m(5) * t33 + t126 * (pkin(8) * t36 - t33) + t121 * t36 + t120 * t35) * g(2) + (m(5) * t34 + t126 * (pkin(8) * t38 - t34) + t121 * t38 + t120 * t37) * g(1) (-t118 * (-t39 * t52 - t40 * t57) + t112 * t40 + t113 * t39) * g(3) + (-t118 * (-t15 * t52 - t16 * t57) + t112 * t16 + t113 * t15) * g(2) + (-t118 * (-t17 * t52 - t18 * t57) + t112 * t18 + t113 * t17) * g(1) (-t73 * mrSges(5,1) - (-t40 * t61 + t58 * t98) * mrSges(5,2) - t127 * (t40 * t54 - t86) + t111 * t11) * g(3) + (-t122 * mrSges(5,1) - (-t16 * t61 - t35 * t58) * mrSges(5,2) - t127 * (t16 * t54 + t35 * t53) + t111 * t1) * g(2) + (-t123 * mrSges(5,1) - (-t18 * t61 - t37 * t58) * mrSges(5,2) - t127 * (t18 * t54 + t37 * t53) + t111 * t3) * g(1) + (-g(1) * t123 - g(2) * t122 - g(3) * t73) * t118 * pkin(4), t118 * (-g(1) * t17 - g(2) * t15 - g(3) * t39) (-g(1) * t3 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
