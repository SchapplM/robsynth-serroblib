% Calculate Gravitation load on the joints for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:42
% EndTime: 2019-03-09 12:27:45
% DurationCPUTime: 1.26s
% Computational Cost: add. (760->128), mult. (1320->175), div. (0->0), fcn. (1524->12), ass. (0->60)
t108 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t59 = cos(qJ(5));
t46 = pkin(5) * t59 + pkin(4);
t119 = -m(6) * pkin(4) - m(7) * t46 - mrSges(5,1);
t64 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t118 = mrSges(6,1) + mrSges(7,1);
t111 = mrSges(6,2) + mrSges(7,2);
t104 = m(7) * pkin(5);
t56 = sin(qJ(5));
t63 = t56 * t104 - t108;
t51 = sin(pkin(11));
t53 = cos(pkin(11));
t113 = -m(4) * pkin(2) - t53 * mrSges(4,1) + t51 * mrSges(4,2) - mrSges(3,1);
t50 = pkin(11) + qJ(4);
t47 = sin(t50);
t48 = cos(t50);
t117 = t119 * t48 + t47 * t64 + t113;
t57 = sin(qJ(2));
t58 = sin(qJ(1));
t60 = cos(qJ(2));
t101 = cos(qJ(1));
t82 = cos(pkin(6));
t71 = t82 * t101;
t32 = t57 * t71 + t58 * t60;
t52 = sin(pkin(6));
t79 = t52 * t101;
t14 = t32 * t48 - t47 * t79;
t31 = t57 * t58 - t60 * t71;
t115 = t14 * t56 - t31 * t59;
t114 = -t14 * t59 - t31 * t56;
t110 = m(5) + m(6) + m(7);
t109 = -t111 * t56 + t118 * t59 - t119;
t107 = -t104 - t118;
t98 = t32 * t56;
t75 = t58 * t82;
t34 = t101 * t60 - t57 * t75;
t97 = t34 * t56;
t94 = t48 * t56;
t93 = t48 * t59;
t92 = t52 * t57;
t91 = t52 * t58;
t90 = t52 * t60;
t89 = t56 * t60;
t88 = t59 * t60;
t83 = t101 * pkin(1) + pkin(8) * t91;
t77 = -pkin(1) * t58 + pkin(8) * t79;
t74 = t51 * t79;
t33 = t101 * t57 + t60 * t75;
t45 = pkin(3) * t53 + pkin(2);
t55 = -pkin(9) - qJ(3);
t73 = t51 * pkin(3) * t91 - t33 * t55 + t34 * t45 + t83;
t18 = t34 * t48 + t47 * t91;
t5 = -t18 * t56 + t33 * t59;
t68 = pkin(3) * t74 + t31 * t55 - t32 * t45 + t77;
t13 = t32 * t47 + t48 * t79;
t26 = t47 * t82 + t48 * t92;
t25 = t47 * t92 - t48 * t82;
t17 = t34 * t47 - t48 * t91;
t6 = t18 * t59 + t33 * t56;
t1 = [(-t101 * mrSges(2,1) - m(5) * t73 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t73) - m(7) * (t18 * t46 + t73) - t118 * t6 + t113 * t34 + (mrSges(2,2) + (-mrSges(4,1) * t51 - mrSges(4,2) * t53 - mrSges(3,3)) * t52) * t58 - t111 * t5 - t63 * t33 + t64 * t17 + (-m(3) - m(4)) * t83) * g(2) + (t58 * mrSges(2,1) + t101 * mrSges(2,2) - m(3) * t77 + t32 * mrSges(3,1) - mrSges(3,3) * t79 - m(4) * (-pkin(2) * t32 + t77) - (-t32 * t53 + t74) * mrSges(4,1) - (t32 * t51 + t53 * t79) * mrSges(4,2) - m(5) * t68 + t14 * mrSges(5,1) - m(6) * (-pkin(4) * t14 + t68) - m(7) * (-t14 * t46 + t68) - t118 * t114 - t111 * t115 + t63 * t31 - t64 * t13) * g(1) (-t98 * t104 - t110 * (-t31 * t45 - t32 * t55) - t118 * (-t31 * t93 + t98) - t111 * (t31 * t94 + t32 * t59) + t108 * t32 - t117 * t31) * g(2) + (-t97 * t104 - t111 * (t33 * t94 + t34 * t59) - t110 * (-t33 * t45 - t34 * t55) - t118 * (-t33 * t93 + t97) + t108 * t34 - t117 * t33) * g(1) + (-t110 * t45 * t90 + (t117 * t60 + (t111 * t89 - t118 * t88) * t48 + (t110 * t55 - t111 * t59 - t118 * t56 - t63) * t57) * t52) * g(3) (-g(1) * t33 - g(2) * t31 + g(3) * t90) * (m(4) + t110) (t109 * t25 + t64 * t26) * g(3) + (t109 * t13 + t64 * t14) * g(2) + (t109 * t17 + t64 * t18) * g(1) (-t111 * (-t26 * t59 + t52 * t89) + t107 * (-t26 * t56 - t52 * t88)) * g(3) + (-t107 * t115 - t111 * t114) * g(2) + (t107 * t5 + t111 * t6) * g(1) (-g(1) * t17 - g(2) * t13 - g(3) * t25) * m(7)];
taug  = t1(:);
