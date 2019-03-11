% Calculate Gravitation load on the joints for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:35
% EndTime: 2019-03-08 20:22:37
% DurationCPUTime: 0.88s
% Computational Cost: add. (661->109), mult. (1372->160), div. (0->0), fcn. (1671->14), ass. (0->62)
t130 = mrSges(6,2) - mrSges(7,3);
t58 = sin(qJ(6));
t61 = cos(qJ(6));
t129 = -t61 * mrSges(7,1) + t58 * mrSges(7,2) - mrSges(6,1);
t122 = -m(6) - m(7);
t123 = -m(4) - m(5);
t128 = t122 + t123;
t53 = sin(pkin(12));
t60 = sin(qJ(2));
t63 = cos(qJ(2));
t95 = cos(pkin(12));
t39 = -t63 * t53 - t60 * t95;
t57 = cos(pkin(6));
t102 = t57 * t63;
t54 = sin(pkin(11));
t56 = cos(pkin(11));
t126 = t56 * t102 - t54 * t60;
t55 = sin(pkin(6));
t62 = cos(qJ(4));
t107 = t55 * t62;
t72 = -t60 * t53 + t63 * t95;
t96 = t39 * t57;
t24 = t54 * t96 + t56 * t72;
t59 = sin(qJ(4));
t125 = t54 * t107 - t24 * t59;
t37 = t39 * t55;
t124 = t37 * t59 + t57 * t62;
t52 = qJ(4) + qJ(5);
t50 = sin(t52);
t51 = cos(t52);
t30 = t37 * t50 + t51 * t57;
t31 = -t37 * t51 + t50 * t57;
t121 = t129 * t30 + t130 * t31;
t111 = t54 * t55;
t13 = t111 * t51 - t24 * t50;
t14 = t111 * t50 + t24 * t51;
t120 = t129 * t13 + t130 * t14;
t109 = t55 * t56;
t19 = -t54 * t72 + t56 * t96;
t11 = -t109 * t51 + t19 * t50;
t12 = -t109 * t50 - t19 * t51;
t119 = t129 * t11 + t130 * t12;
t118 = -m(5) * pkin(3) - t62 * mrSges(5,1) + t59 * mrSges(5,2) - mrSges(4,1) + (-m(7) * pkin(5) + t129) * t51 + (-m(7) * pkin(10) + t130) * t50;
t117 = m(5) * pkin(8) + t58 * mrSges(7,1) + t61 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t108 = t55 * t59;
t104 = t57 * t60;
t46 = t55 * t63 * pkin(2);
t88 = t11 * pkin(5) + pkin(10) * t12;
t87 = t13 * pkin(5) + pkin(10) * t14;
t86 = t30 * pkin(5) + pkin(10) * t31;
t83 = t126 * pkin(2);
t82 = t125 * pkin(4);
t81 = t124 * pkin(4);
t75 = -t107 * t56 + t19 * t59;
t69 = t75 * pkin(4);
t67 = t57 * t72;
t64 = -pkin(9) - pkin(8);
t49 = pkin(4) * t62 + pkin(3);
t36 = t72 * t55;
t23 = t39 * t56 - t54 * t67;
t20 = t54 * t39 + t56 * t67;
t1 = [(-m(2) - m(3) + t128) * g(3) (-(mrSges(3,1) * t63 - mrSges(3,2) * t60) * t55 + t122 * (t36 * t49 + t37 * t64 + t46) + t123 * t46 + t117 * t37 + t118 * t36) * g(3) + (-t126 * mrSges(3,1) - (-t104 * t56 - t54 * t63) * mrSges(3,2) + t123 * t83 + t122 * (t19 * t64 + t20 * t49 + t83) + t118 * t20 + t117 * t19) * g(2) + (-(t104 * t54 - t56 * t63) * mrSges(3,2) + t122 * (t23 * t49 - t24 * t64) + t118 * t23 - t117 * t24 + (t128 * pkin(2) - mrSges(3,1)) * (-t54 * t102 - t56 * t60)) * g(1) -(-g(3) * t57 + (-g(1) * t54 + g(2) * t56) * t55) * t128 (-t124 * mrSges(5,1) - (t37 * t62 - t57 * t59) * mrSges(5,2) - m(6) * t81 - m(7) * (t81 + t86) + t121) * g(3) + (-t75 * mrSges(5,1) - (t108 * t56 + t19 * t62) * mrSges(5,2) - m(6) * t69 - m(7) * (t69 + t88) + t119) * g(2) + (-t125 * mrSges(5,1) - (-t108 * t54 - t24 * t62) * mrSges(5,2) - m(6) * t82 - m(7) * (t82 + t87) + t120) * g(1) (-m(7) * t86 + t121) * g(3) + (-m(7) * t88 + t119) * g(2) + (-m(7) * t87 + t120) * g(1), -g(1) * ((-t14 * t58 - t23 * t61) * mrSges(7,1) + (-t14 * t61 + t23 * t58) * mrSges(7,2)) - g(2) * ((-t12 * t58 - t20 * t61) * mrSges(7,1) + (-t12 * t61 + t20 * t58) * mrSges(7,2)) - g(3) * ((-t31 * t58 - t36 * t61) * mrSges(7,1) + (-t31 * t61 + t36 * t58) * mrSges(7,2))];
taug  = t1(:);
