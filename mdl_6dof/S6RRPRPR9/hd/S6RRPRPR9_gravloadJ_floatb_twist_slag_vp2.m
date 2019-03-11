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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:37
% EndTime: 2019-03-09 10:55:41
% DurationCPUTime: 1.27s
% Computational Cost: add. (755->108), mult. (1239->147), div. (0->0), fcn. (1424->14), ass. (0->49)
t40 = pkin(12) + qJ(6);
t35 = sin(t40);
t37 = cos(t40);
t42 = sin(pkin(12));
t45 = cos(pkin(12));
t54 = -mrSges(5,1) - m(6) * pkin(4) - t45 * mrSges(6,1) + t42 * mrSges(6,2) - m(7) * (pkin(5) * t45 + pkin(4));
t113 = -t37 * mrSges(7,1) + t35 * mrSges(7,2) + t54;
t57 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t93 = -(m(7) * pkin(5) + mrSges(6,1)) * t42 - m(4) * qJ(3) - mrSges(6,2) * t45 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t112 = t35 * mrSges(7,1) + t37 * mrSges(7,2) - t93;
t43 = sin(pkin(11));
t46 = cos(pkin(11));
t101 = -m(4) * pkin(2) - t46 * mrSges(4,1) + t43 * mrSges(4,2) - mrSges(3,1);
t41 = pkin(11) + qJ(4);
t36 = sin(t41);
t38 = cos(t41);
t92 = -t113 * t38 - t57 * t36 - t101;
t102 = m(6) + m(7);
t96 = m(5) + t102;
t105 = m(4) + t96;
t86 = cos(qJ(1));
t44 = sin(pkin(6));
t49 = sin(qJ(2));
t82 = t44 * t49;
t50 = sin(qJ(1));
t81 = t44 * t50;
t51 = cos(qJ(2));
t80 = t44 * t51;
t77 = cos(pkin(6));
t74 = t44 * t86;
t72 = -pkin(1) * t50 + pkin(8) * t74;
t66 = t77 * t86;
t20 = t49 * t66 + t50 * t51;
t4 = t20 * t38 - t36 * t74;
t69 = t50 * t77;
t68 = t43 * t74;
t3 = t20 * t36 + t38 * t74;
t48 = -pkin(9) - qJ(3);
t34 = pkin(3) * t46 + pkin(2);
t22 = -t49 * t69 + t51 * t86;
t21 = t49 * t86 + t51 * t69;
t19 = t49 * t50 - t51 * t66;
t14 = t36 * t77 + t38 * t82;
t13 = t36 * t82 - t38 * t77;
t8 = t22 * t38 + t36 * t81;
t7 = t22 * t36 - t38 * t81;
t2 = t21 * t35 + t37 * t8;
t1 = t21 * t37 - t35 * t8;
t5 = [(-t86 * mrSges(2,1) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t101 * t22 + (mrSges(2,2) + (-mrSges(4,1) * t43 - mrSges(4,2) * t46 - mrSges(3,3)) * t44) * t50 + t54 * t8 + t57 * t7 + t93 * t21 - t96 * (t43 * pkin(3) * t81 - t21 * t48 + t22 * t34) + (-m(3) - t105) * (t86 * pkin(1) + pkin(8) * t81)) * g(2) + (t50 * mrSges(2,1) + t86 * mrSges(2,2) - m(3) * t72 + t20 * mrSges(3,1) - mrSges(3,3) * t74 - m(4) * (-pkin(2) * t20 + t72) - (-t20 * t46 + t68) * mrSges(4,1) - (t20 * t43 + t46 * t74) * mrSges(4,2) - t57 * t3 - t113 * t4 + t112 * t19 + t96 * (-pkin(3) * t68 - t19 * t48 + t20 * t34 - t72)) * g(1) (-t96 * (-t19 * t34 - t20 * t48) - t112 * t20 + t92 * t19) * g(2) + (-t96 * (-t21 * t34 - t22 * t48) - t112 * t22 + t92 * t21) * g(1) + (-t96 * t34 * t80 + (-t92 * t51 + (t96 * t48 - t112) * t49) * t44) * g(3) (-g(1) * t21 - g(2) * t19 + g(3) * t80) * t105 (-t113 * t13 + t57 * t14) * g(3) + (-t113 * t3 + t57 * t4) * g(2) + (-t113 * t7 + t57 * t8) * g(1), t102 * (-g(1) * t7 - g(2) * t3 - g(3) * t13) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t19 * t37 - t35 * t4) * mrSges(7,1) + (-t19 * t35 - t37 * t4) * mrSges(7,2)) - g(3) * ((-t14 * t35 - t37 * t80) * mrSges(7,1) + (-t14 * t37 + t35 * t80) * mrSges(7,2))];
taug  = t5(:);
