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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:26
% EndTime: 2019-03-09 11:04:29
% DurationCPUTime: 1.20s
% Computational Cost: add. (692->125), mult. (1192->169), div. (0->0), fcn. (1363->12), ass. (0->63)
t53 = sin(qJ(6));
t56 = cos(qJ(6));
t117 = t53 * mrSges(7,1) + t56 * mrSges(7,2);
t110 = mrSges(5,2) - mrSges(6,3);
t113 = m(6) + m(7);
t64 = -qJ(5) * t113 + t110;
t103 = t64 - t117;
t58 = -m(4) * qJ(3) - m(7) * pkin(5) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t121 = t56 * mrSges(7,1) - t53 * mrSges(7,2) - t58;
t120 = m(5) + t113;
t119 = -mrSges(5,1) + mrSges(6,2);
t49 = sin(pkin(11));
t51 = cos(pkin(11));
t112 = -m(4) * pkin(2) - t51 * mrSges(4,1) + t49 * mrSges(4,2) - mrSges(3,1);
t48 = pkin(11) + qJ(4);
t45 = sin(t48);
t46 = cos(t48);
t118 = -t110 * t45 - t119 * t46 - t112;
t107 = -m(7) * pkin(10) - mrSges(7,3);
t104 = -t107 - t119;
t114 = pkin(4) * t113 + t104;
t101 = t46 * mrSges(7,3) + t117 * t45 + t118;
t97 = cos(qJ(1));
t54 = sin(qJ(2));
t55 = sin(qJ(1));
t57 = cos(qJ(2));
t81 = cos(pkin(6));
t70 = t81 * t97;
t28 = t54 * t55 - t57 * t70;
t96 = t28 * t46;
t75 = t55 * t81;
t30 = t54 * t97 + t57 * t75;
t95 = t30 * t46;
t50 = sin(pkin(6));
t90 = t50 * t54;
t89 = t50 * t55;
t88 = t50 * t57;
t87 = t53 * t57;
t86 = t56 * t57;
t29 = t54 * t70 + t55 * t57;
t44 = pkin(3) * t51 + pkin(2);
t52 = -pkin(9) - qJ(3);
t85 = -t28 * t44 - t29 * t52;
t31 = -t54 * t75 + t57 * t97;
t84 = -t30 * t44 - t31 * t52;
t83 = t97 * pkin(1) + pkin(8) * t89;
t82 = qJ(5) * t45;
t78 = t50 * t97;
t77 = -pkin(1) * t55 + pkin(8) * t78;
t8 = t29 * t46 - t45 * t78;
t74 = -pkin(4) * t96 - t28 * t82 + t85;
t73 = -pkin(4) * t95 - t30 * t82 + t84;
t72 = t49 * t78;
t71 = t49 * pkin(3) * t89 - t30 * t52 + t31 * t44 + t83;
t63 = pkin(3) * t72 + t28 * t52 - t29 * t44 + t77;
t7 = t29 * t45 + t46 * t78;
t33 = t44 * t88;
t22 = t45 * t90 - t46 * t81;
t12 = t31 * t46 + t45 * t89;
t11 = t31 * t45 - t46 * t89;
t2 = t11 * t53 + t30 * t56;
t1 = t11 * t56 - t30 * t53;
t3 = [(-t97 * mrSges(2,1) - m(5) * t71 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t112 * t31 + (mrSges(2,2) + (-mrSges(4,1) * t49 - mrSges(4,2) * t51 - mrSges(3,3)) * t50) * t55 + t64 * t11 - t104 * t12 + t58 * t30 + (-m(3) - m(4)) * t83 - t113 * (t12 * pkin(4) + t71)) * g(2) + (t55 * mrSges(2,1) + t97 * mrSges(2,2) - m(3) * t77 + t29 * mrSges(3,1) - mrSges(3,3) * t78 - m(4) * (-pkin(2) * t29 + t77) - (-t29 * t51 + t72) * mrSges(4,1) - (t29 * t49 + t51 * t78) * mrSges(4,2) - m(5) * t63 + t104 * t8 - t103 * t7 + t121 * t28 + t113 * (pkin(4) * t8 - t63)) * g(1) (-m(5) * t85 - m(6) * t74 - m(7) * (-pkin(10) * t96 + t74) - t121 * t29 + t101 * t28) * g(2) + (-m(5) * t84 - m(6) * t73 - m(7) * (-pkin(10) * t95 + t73) - t121 * t31 + t101 * t30) * g(1) + (-m(5) * t33 - t113 * (t33 + (pkin(4) * t46 + t82) * t88) + ((-t87 * mrSges(7,1) - t86 * mrSges(7,2)) * t45 + (t107 * t46 - t118) * t57 + (t120 * t52 - t121) * t54) * t50) * g(3) (-g(1) * t30 - g(2) * t28 + g(3) * t88) * (m(4) + t120) (t103 * (t45 * t81 + t46 * t90) + t114 * t22) * g(3) + (t103 * t8 + t114 * t7) * g(2) + (t103 * t12 + t11 * t114) * g(1), t113 * (-g(1) * t11 - g(2) * t7 - g(3) * t22) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t28 * t53 + t56 * t7) * mrSges(7,1) + (-t28 * t56 - t53 * t7) * mrSges(7,2)) - g(3) * ((t22 * t56 + t50 * t87) * mrSges(7,1) + (-t22 * t53 + t50 * t86) * mrSges(7,2))];
taug  = t3(:);
