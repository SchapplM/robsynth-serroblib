% Calculate Gravitation load on the joints for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:06:53
% EndTime: 2018-11-23 17:06:54
% DurationCPUTime: 1.03s
% Computational Cost: add. (1340->128), mult. (1408->165), div. (0->0), fcn. (1363->16), ass. (0->66)
t117 = mrSges(5,2) - mrSges(6,3);
t119 = m(6) + m(7);
t76 = -qJ(5) * t119 + t117;
t62 = sin(qJ(6));
t65 = cos(qJ(6));
t80 = -t62 * mrSges(7,1) - t65 * mrSges(7,2);
t114 = t76 + t80;
t123 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t115 = m(7) * pkin(10) + t123;
t120 = pkin(4) * t119 + t115;
t58 = sin(pkin(11));
t60 = cos(pkin(11));
t118 = -m(4) * pkin(2) - t60 * mrSges(4,1) + t58 * mrSges(4,2) - mrSges(3,1);
t68 = -m(4) * qJ(3) - m(7) * pkin(5) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t112 = t65 * mrSges(7,1) - t62 * mrSges(7,2) - t68;
t57 = pkin(11) + qJ(4);
t54 = sin(t57);
t55 = cos(t57);
t113 = -t118 + t123 * t55 + (-t80 - t117) * t54;
t110 = cos(qJ(1));
t63 = sin(qJ(2));
t64 = sin(qJ(1));
t95 = pkin(6) + qJ(2);
t83 = cos(t95) / 0.2e1;
t96 = pkin(6) - qJ(2);
t87 = cos(t96);
t67 = t87 / 0.2e1 + t83;
t28 = -t110 * t67 + t63 * t64;
t109 = t28 * t55;
t31 = t110 * t63 + t64 * t67;
t108 = t31 * t55;
t82 = sin(t95) / 0.2e1;
t86 = sin(t96);
t40 = t82 + t86 / 0.2e1;
t107 = t40 * t55;
t59 = sin(pkin(6));
t103 = t59 * t64;
t53 = pkin(3) * t60 + pkin(2);
t61 = -pkin(9) - qJ(3);
t41 = t82 - t86 / 0.2e1;
t66 = cos(qJ(2));
t74 = t110 * t41 + t64 * t66;
t102 = -t28 * t53 - t61 * t74;
t73 = t110 * t66 - t64 * t41;
t101 = -t31 * t53 - t61 * t73;
t42 = t83 - t87 / 0.2e1;
t100 = t40 * t53 + t42 * t61;
t99 = t110 * pkin(1) + pkin(8) * t103;
t98 = qJ(5) * t54;
t97 = cos(pkin(6));
t93 = t59 * t110;
t92 = -pkin(1) * t64 + pkin(8) * t93;
t8 = -t54 * t93 + t55 * t74;
t90 = -pkin(4) * t109 - t28 * t98 + t102;
t89 = -pkin(4) * t108 - t31 * t98 + t101;
t88 = pkin(4) * t107 + t40 * t98 + t100;
t85 = t58 * t93;
t84 = t58 * pkin(3) * t103 - t31 * t61 + t53 * t73 + t99;
t75 = pkin(3) * t85 + t28 * t61 - t53 * t74 + t92;
t7 = t54 * t74 + t55 * t93;
t22 = -t42 * t54 - t55 * t97;
t12 = t103 * t54 + t55 * t73;
t11 = -t103 * t55 + t54 * t73;
t2 = t11 * t62 + t31 * t65;
t1 = t11 * t65 - t31 * t62;
t3 = [(-t110 * mrSges(2,1) - m(5) * t84 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t118 * t73 + (mrSges(2,2) + (-mrSges(4,1) * t58 - mrSges(4,2) * t60 - mrSges(3,3)) * t59) * t64 + t76 * t11 - t115 * t12 + t68 * t31 + (-m(3) - m(4)) * t99 - t119 * (t12 * pkin(4) + t84)) * g(2) + (t64 * mrSges(2,1) + t110 * mrSges(2,2) - m(3) * t92 + t74 * mrSges(3,1) - mrSges(3,3) * t93 - m(4) * (-pkin(2) * t74 + t92) - (-t60 * t74 + t85) * mrSges(4,1) - (t58 * t74 + t60 * t93) * mrSges(4,2) - m(5) * t75 + t115 * t8 - t114 * t7 + t112 * t28 + t119 * (pkin(4) * t8 - t75)) * g(1) (-m(5) * t100 - m(6) * t88 - m(7) * (pkin(10) * t107 + t88) + t112 * t42 - t113 * t40) * g(3) + (-m(5) * t102 - m(6) * t90 - m(7) * (-pkin(10) * t109 + t90) - t112 * t74 + t113 * t28) * g(2) + (-m(5) * t101 - m(6) * t89 - m(7) * (-pkin(10) * t108 + t89) - t112 * t73 + t113 * t31) * g(1) (-g(1) * t31 - g(2) * t28 + g(3) * t40) * (m(4) + m(5) + t119) (t114 * (-t42 * t55 + t54 * t97) + t120 * t22) * g(3) + (t114 * t8 + t120 * t7) * g(2) + (t120 * t11 + t114 * t12) * g(1), t119 * (-g(1) * t11 - g(2) * t7 - g(3) * t22) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t28 * t62 + t65 * t7) * mrSges(7,1) + (-t28 * t65 - t62 * t7) * mrSges(7,2)) - g(3) * ((t22 * t65 + t40 * t62) * mrSges(7,1) + (-t22 * t62 + t40 * t65) * mrSges(7,2))];
taug  = t3(:);
