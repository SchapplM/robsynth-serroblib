% Calculate Gravitation load on the joints for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:17
% EndTime: 2019-03-09 11:17:21
% DurationCPUTime: 1.22s
% Computational Cost: add. (577->128), mult. (1116->170), div. (0->0), fcn. (1260->12), ass. (0->63)
t76 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t113 = m(6) + m(7);
t114 = m(4) + m(5);
t118 = -t114 - t113;
t120 = t118 * qJ(3);
t53 = sin(qJ(6));
t57 = cos(qJ(6));
t112 = m(7) * pkin(5) + mrSges(7,1) * t57 - mrSges(7,2) * t53 + mrSges(6,1);
t50 = qJ(4) + pkin(11);
t47 = sin(t50);
t48 = cos(t50);
t54 = sin(qJ(4));
t58 = cos(qJ(4));
t74 = t54 * mrSges(5,1) + t58 * mrSges(5,2);
t96 = mrSges(3,2) - mrSges(4,3);
t119 = -t112 * t47 - t76 * t48 - t74 + t96;
t69 = -m(5) * pkin(9) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t107 = t53 * mrSges(7,1) + t57 * mrSges(7,2) - t69;
t55 = sin(qJ(2));
t56 = sin(qJ(1));
t59 = cos(qJ(2));
t60 = cos(qJ(1));
t92 = cos(pkin(6));
t79 = t60 * t92;
t30 = t55 * t56 - t59 * t79;
t51 = sin(pkin(6));
t97 = t51 * t60;
t116 = t30 * t58 + t54 * t97;
t80 = t56 * t92;
t32 = t60 * t55 + t59 * t80;
t99 = t51 * t56;
t9 = t32 * t58 - t54 * t99;
t98 = t51 * t59;
t64 = -t92 * t54 - t58 * t98;
t110 = t96 + t120;
t106 = t119 + t120;
t105 = pkin(4) * t54;
t104 = t30 * t54;
t102 = t32 * t54;
t100 = t51 * t55;
t94 = pkin(2) * t98 + qJ(3) * t100;
t93 = t60 * pkin(1) + pkin(8) * t99;
t33 = -t55 * t80 + t59 * t60;
t86 = t33 * pkin(2) + t93;
t81 = -t56 * pkin(1) + pkin(8) * t97;
t77 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(3,3);
t31 = t55 * t79 + t56 * t59;
t75 = t31 * pkin(2) - t81;
t46 = pkin(4) * t58 + pkin(3);
t52 = -qJ(5) - pkin(9);
t70 = pkin(4) * t102 - t33 * t52 + t46 * t99 + t86;
t7 = -t30 * t47 + t48 * t97;
t5 = t30 * t48 + t47 * t97;
t40 = t58 * t97;
t28 = t32 * pkin(2);
t26 = t30 * pkin(2);
t15 = -t47 * t98 + t48 * t92;
t10 = t58 * t99 + t102;
t4 = t32 * t47 + t48 * t99;
t3 = -t32 * t48 + t47 * t99;
t2 = t33 * t53 + t4 * t57;
t1 = t33 * t57 - t4 * t53;
t6 = [(-t60 * mrSges(2,1) - m(3) * t93 - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t70 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t70) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t110 * t32 + t76 * t3 + (t51 * t77 + mrSges(2,2)) * t56 + t69 * t33 - t114 * t86) * g(2) + (t56 * mrSges(2,1) - m(3) * t81 - t40 * mrSges(5,1) + t76 * t5 - t112 * t7 + (-t110 + t74) * t30 + (mrSges(2,2) + (mrSges(5,2) * t54 + t77) * t51) * t60 + t107 * t31 + t114 * t75 + t113 * (pkin(4) * t104 - t31 * t52 - t46 * t97 + t75)) * g(1) (-t114 * t94 - t113 * (t100 * t105 + t94) + ((t113 * t52 - t107) * t59 + t119 * t55) * t51) * g(3) + (-t113 * (t31 * t105 + t30 * t52 - t26) + t114 * t26 + t106 * t31 + t107 * t30) * g(2) + (-t113 * (t33 * t105 + t32 * t52 - t28) + t114 * t28 + t106 * t33 + t107 * t32) * g(1) -(-g(1) * t32 - g(2) * t30 + g(3) * t98) * t118 (-t64 * mrSges(5,1) - (t54 * t98 - t58 * t92) * mrSges(5,2) - t112 * (-t47 * t92 - t48 * t98) + t76 * t15) * g(3) + (-t116 * mrSges(5,1) - (t40 - t104) * mrSges(5,2) - t76 * t7 - t112 * t5) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) + t112 * t3 + t4 * t76) * g(1) + (-g(1) * t9 - g(2) * t116 - g(3) * t64) * t113 * pkin(4), t113 * (-g(1) * t33 - g(2) * t31 - g(3) * t100) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t31 * t57 + t53 * t7) * mrSges(7,1) + (-t31 * t53 + t57 * t7) * mrSges(7,2)) - g(3) * ((t100 * t57 - t15 * t53) * mrSges(7,1) + (-t100 * t53 - t15 * t57) * mrSges(7,2))];
taug  = t6(:);
