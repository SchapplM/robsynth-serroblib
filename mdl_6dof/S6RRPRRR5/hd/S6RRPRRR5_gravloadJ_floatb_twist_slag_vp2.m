% Calculate Gravitation load on the joints for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:04
% EndTime: 2019-03-09 13:38:07
% DurationCPUTime: 1.30s
% Computational Cost: add. (971->135), mult. (2259->184), div. (0->0), fcn. (2818->14), ass. (0->70)
t62 = cos(qJ(5));
t51 = pkin(5) * t62 + pkin(4);
t55 = qJ(5) + qJ(6);
t53 = sin(t55);
t54 = cos(t55);
t58 = sin(qJ(5));
t122 = m(6) * pkin(4) + m(7) * t51 + t62 * mrSges(6,1) + t54 * mrSges(7,1) - t58 * mrSges(6,2) - t53 * mrSges(7,2) + mrSges(5,1);
t74 = m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t132 = -m(5) - m(6);
t129 = t53 * mrSges(7,1) + t62 * mrSges(6,2) + t54 * mrSges(7,2);
t59 = sin(qJ(4));
t63 = cos(qJ(4));
t126 = -t122 * t63 + t74 * t59 - mrSges(4,1);
t115 = pkin(5) * t58;
t128 = mrSges(4,2) - mrSges(5,3);
t124 = t58 * mrSges(6,1) - t128 - t132 * pkin(9) + m(7) * (pkin(9) + t115) + t129;
t113 = cos(qJ(2));
t60 = sin(qJ(2));
t98 = sin(pkin(12));
t99 = cos(pkin(12));
t72 = t113 * t98 + t60 * t99;
t61 = sin(qJ(1));
t108 = t61 * t60;
t57 = cos(pkin(6));
t64 = cos(qJ(1));
t90 = t64 * t113;
t35 = -t57 * t90 + t108;
t97 = m(7) - t132;
t39 = -t113 * t99 + t60 * t98;
t123 = -m(7) * pkin(5) - mrSges(6,1);
t120 = -t97 * pkin(9) + t128;
t100 = t72 * t57;
t26 = -t100 * t61 - t64 * t39;
t21 = -t100 * t64 + t61 * t39;
t56 = sin(pkin(6));
t111 = t56 * t64;
t12 = -t59 * t111 - t21 * t63;
t69 = t57 * t39;
t22 = -t61 * t72 - t64 * t69;
t117 = (-t12 * t53 - t22 * t54) * mrSges(7,1) + (-t12 * t54 + t22 * t53) * mrSges(7,2);
t112 = t56 * t61;
t16 = t112 * t59 + t26 * t63;
t25 = t61 * t69 - t64 * t72;
t5 = -t16 * t53 - t25 * t54;
t6 = t16 * t54 - t25 * t53;
t116 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t33 = t72 * t56;
t28 = t33 * t63 + t57 * t59;
t32 = t39 * t56;
t114 = (-t28 * t53 + t32 * t54) * mrSges(7,1) + (-t28 * t54 - t32 * t53) * mrSges(7,2);
t104 = t64 * t60;
t96 = t113 * pkin(2);
t52 = t96 + pkin(1);
t45 = t64 * t52;
t102 = t26 * pkin(3) + t45;
t95 = m(4) + t97;
t94 = -m(3) * pkin(1) - mrSges(2,1);
t91 = t61 * t113;
t11 = -t63 * t111 + t21 * t59;
t49 = t56 * t96;
t85 = t35 * pkin(2);
t7 = -t16 * t58 - t25 * t62;
t37 = -t57 * t91 - t104;
t73 = t37 * pkin(2);
t68 = mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t56 + t95 * (pkin(2) * t57 * t60 + (-pkin(8) - qJ(3)) * t56);
t38 = -t108 * t57 + t90;
t36 = -t104 * t57 - t91;
t15 = -t112 * t63 + t26 * t59;
t8 = t16 * t62 - t25 * t58;
t1 = [(-t38 * mrSges(3,1) - t37 * mrSges(3,2) - m(4) * t45 - t26 * mrSges(4,1) - m(5) * t102 - t16 * mrSges(5,1) - m(6) * (pkin(4) * t16 + t102) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t51 + t102) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t94 * t64 + t68 * t61 + (m(7) * t115 - t120) * t25 + t74 * t15) * g(2) + (-t36 * mrSges(3,1) - t35 * mrSges(3,2) + (t52 * t95 - t94) * t61 + t68 * t64 + t122 * t12 + (t123 * t58 + t120 - t129) * t22 + t74 * t11 - (t97 * pkin(3) + mrSges(4,1)) * t21) * g(1) (-(mrSges(3,1) * t113 - mrSges(3,2) * t60) * t56 - m(4) * t49 - t97 * (-t32 * pkin(3) + t49) - t124 * t33 - t126 * t32) * g(3) + (m(4) * t85 + mrSges(3,1) * t35 - mrSges(3,2) * t36 - t97 * (t22 * pkin(3) - t85) + t126 * t22 + t124 * t21) * g(2) + (-m(4) * t73 - mrSges(3,1) * t37 + mrSges(3,2) * t38 - t97 * (t25 * pkin(3) + t73) + t126 * t25 - t124 * t26) * g(1) (-g(3) * t57 + (-g(1) * t61 + g(2) * t64) * t56) * t95 (t74 * t28 - t122 * (-t33 * t59 + t57 * t63)) * g(3) + (-t11 * t122 + t12 * t74) * g(2) + (t122 * t15 + t16 * t74) * g(1) (-(-t28 * t62 - t32 * t58) * mrSges(6,2) - t114 + t123 * (-t28 * t58 + t32 * t62)) * g(3) + (-(-t12 * t62 + t22 * t58) * mrSges(6,2) - t117 + t123 * (-t12 * t58 - t22 * t62)) * g(2) + (mrSges(6,2) * t8 + t123 * t7 - t116) * g(1), -g(1) * t116 - g(2) * t117 - g(3) * t114];
taug  = t1(:);
