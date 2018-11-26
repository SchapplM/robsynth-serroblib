% Calculate Gravitation load on the joints for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:05:57
% EndTime: 2018-11-23 17:05:58
% DurationCPUTime: 1.24s
% Computational Cost: add. (1418->114), mult. (1460->148), div. (0->0), fcn. (1424->18), ass. (0->55)
t49 = pkin(12) + qJ(6);
t44 = sin(t49);
t46 = cos(t49);
t51 = sin(pkin(12));
t54 = cos(pkin(12));
t63 = -mrSges(5,1) - m(6) * pkin(4) - mrSges(6,1) * t54 + mrSges(6,2) * t51 - m(7) * (pkin(5) * t54 + pkin(4));
t124 = -mrSges(7,1) * t46 + mrSges(7,2) * t44 + t63;
t67 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t52 = sin(pkin(11));
t55 = cos(pkin(11));
t115 = -m(4) * pkin(2) - t55 * mrSges(4,1) + t52 * mrSges(4,2) - mrSges(3,1);
t50 = pkin(11) + qJ(4);
t45 = sin(t50);
t47 = cos(t50);
t106 = -t124 * t47 - t67 * t45 - t115;
t107 = -m(4) * qJ(3) - mrSges(6,2) * t54 - t51 * (m(7) * pkin(5) + mrSges(6,1)) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t122 = t44 * mrSges(7,1) + t46 * mrSges(7,2) - t107;
t116 = m(6) + m(7);
t111 = m(5) + t116;
t118 = m(4) + t111;
t101 = cos(qJ(1));
t53 = sin(pkin(6));
t59 = sin(qJ(1));
t97 = t53 * t59;
t93 = cos(pkin(6));
t92 = pkin(6) - qJ(2);
t91 = pkin(6) + qJ(2);
t88 = t53 * t101;
t86 = -pkin(1) * t59 + pkin(8) * t88;
t78 = sin(t91) / 0.2e1;
t82 = sin(t92);
t30 = t78 - t82 / 0.2e1;
t60 = cos(qJ(2));
t73 = t101 * t30 + t59 * t60;
t4 = -t45 * t88 + t47 * t73;
t83 = cos(t92);
t81 = t52 * t88;
t79 = cos(t91) / 0.2e1;
t72 = t101 * t60 - t59 * t30;
t3 = t45 * t73 + t47 * t88;
t64 = t83 / 0.2e1 + t79;
t58 = sin(qJ(2));
t57 = -pkin(9) - qJ(3);
t43 = pkin(3) * t55 + pkin(2);
t31 = t79 - t83 / 0.2e1;
t29 = t78 + t82 / 0.2e1;
t22 = t101 * t58 + t59 * t64;
t19 = -t101 * t64 + t58 * t59;
t14 = -t31 * t47 + t45 * t93;
t13 = -t31 * t45 - t47 * t93;
t8 = t45 * t97 + t47 * t72;
t7 = t45 * t72 - t47 * t97;
t2 = t22 * t44 + t46 * t8;
t1 = t22 * t46 - t44 * t8;
t5 = [(-t101 * mrSges(2,1) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t115 * t72 + (mrSges(2,2) + (-mrSges(4,1) * t52 - mrSges(4,2) * t55 - mrSges(3,3)) * t53) * t59 + t63 * t8 + t67 * t7 + t107 * t22 - t111 * (t52 * pkin(3) * t97 - t22 * t57 + t43 * t72) + (-m(3) - t118) * (t101 * pkin(1) + pkin(8) * t97)) * g(2) + (t59 * mrSges(2,1) + t101 * mrSges(2,2) - m(3) * t86 + t73 * mrSges(3,1) - mrSges(3,3) * t88 - m(4) * (-pkin(2) * t73 + t86) - (-t55 * t73 + t81) * mrSges(4,1) - (t52 * t73 + t55 * t88) * mrSges(4,2) - t67 * t3 - t124 * t4 + t122 * t19 + t111 * (-pkin(3) * t81 - t19 * t57 + t43 * t73 - t86)) * g(1) (-t111 * (t29 * t43 + t31 * t57) + t122 * t31 - t106 * t29) * g(3) + (-t111 * (-t19 * t43 - t57 * t73) - t122 * t73 + t106 * t19) * g(2) + (-t111 * (-t22 * t43 - t57 * t72) - t122 * t72 + t106 * t22) * g(1) (-g(1) * t22 - g(2) * t19 + g(3) * t29) * t118 (-t124 * t13 + t67 * t14) * g(3) + (-t124 * t3 + t67 * t4) * g(2) + (-t124 * t7 + t67 * t8) * g(1), t116 * (-g(1) * t7 - g(2) * t3 - g(3) * t13) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t19 * t46 - t4 * t44) * mrSges(7,1) + (-t19 * t44 - t4 * t46) * mrSges(7,2)) - g(3) * ((-t14 * t44 - t29 * t46) * mrSges(7,1) + (-t14 * t46 + t29 * t44) * mrSges(7,2))];
taug  = t5(:);
