% Calculate Gravitation load on the joints for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:03
% EndTime: 2019-12-31 20:29:05
% DurationCPUTime: 0.71s
% Computational Cost: add. (226->97), mult. (521->124), div. (0->0), fcn. (531->8), ass. (0->50)
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t78 = m(6) * pkin(4) + t28 * mrSges(6,1) - t25 * mrSges(6,2) + mrSges(5,1);
t27 = sin(qJ(1));
t29 = cos(qJ(2));
t54 = qJ(3) * t29;
t13 = t27 * t54;
t26 = sin(qJ(2));
t64 = -pkin(2) - pkin(3);
t75 = t26 * t64;
t77 = t27 * t75 + t13;
t30 = cos(qJ(1));
t14 = t30 * t54;
t76 = t30 * t75 + t14;
t74 = (mrSges(3,1) + mrSges(4,1)) * t29 + (-mrSges(3,2) + mrSges(4,3)) * t26;
t73 = -mrSges(2,1) - t74;
t50 = m(6) * pkin(8) + mrSges(6,3);
t72 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t71 = -mrSges(5,2) + t50;
t23 = t30 * pkin(6);
t19 = t26 * qJ(3);
t44 = -pkin(1) - t19;
t70 = t27 * (t64 * t29 + t44) - t30 * pkin(7) + t23;
t60 = sin(qJ(4));
t61 = cos(qJ(4));
t8 = t26 * t61 - t29 * t60;
t69 = g(1) * t30 + g(2) * t27;
t3 = t8 * t27;
t7 = t26 * t60 + t29 * t61;
t4 = t7 * t27;
t67 = -t4 * mrSges(5,2) + t78 * t3;
t66 = -t8 * mrSges(5,2) - t78 * t7;
t46 = t30 * t60;
t47 = t30 * t61;
t5 = -t26 * t46 - t29 * t47;
t6 = -t26 * t47 + t29 * t46;
t65 = t5 * mrSges(5,2) - t78 * t6;
t22 = t29 * pkin(2);
t57 = t29 * t30;
t56 = t22 + t19;
t55 = t30 * pkin(1) + t27 * pkin(6);
t52 = t29 * pkin(3) + t56;
t43 = pkin(2) * t57 + t30 * t19 + t55;
t39 = t4 * t25 - t30 * t28;
t38 = -t30 * t25 - t4 * t28;
t34 = pkin(3) * t57 - t27 * pkin(7) + t43;
t32 = t29 * mrSges(4,3) + (-m(4) * pkin(2) - mrSges(4,1)) * t26;
t2 = -t27 * t25 - t5 * t28;
t1 = t5 * t25 - t27 * t28;
t9 = [(-m(3) * t55 - m(4) * t43 - m(5) * t34 + t5 * mrSges(5,1) - m(6) * (-t5 * pkin(4) + t34) - t2 * mrSges(6,1) - t1 * mrSges(6,2) - t71 * t6 + t73 * t30 + t72 * t27) * g(2) + (-t70 * m(5) + t4 * mrSges(5,1) - t38 * mrSges(6,1) - t39 * mrSges(6,2) - (-t4 * pkin(4) + t70) * m(6) - t71 * t3 + (-m(3) - m(4)) * t23 + (m(3) * pkin(1) - m(4) * (t44 - t22) - t73) * t27 + t72 * t30) * g(1), t69 * (mrSges(3,1) * t26 + mrSges(3,2) * t29) + (-m(4) * t13 - t32 * t27 - m(5) * t77 - m(6) * (-t4 * pkin(8) + t77) + t4 * mrSges(6,3) + t67) * g(2) + (-m(4) * t14 - t32 * t30 - m(5) * t76 - m(6) * (t5 * pkin(8) + t76) - t5 * mrSges(6,3) + t65) * g(1) + (-m(4) * t56 - m(5) * t52 - m(6) * (-t8 * pkin(8) + t52) + t8 * mrSges(6,3) + t66 - t74) * g(3), (t29 * g(3) - t69 * t26) * (m(4) + m(5) + m(6)), (-t50 * t8 - t66) * g(3) + (-t50 * t4 - t67) * g(2) + (t50 * t5 - t65) * g(1), -g(1) * (t1 * mrSges(6,1) - t2 * mrSges(6,2)) - g(2) * (-t39 * mrSges(6,1) + t38 * mrSges(6,2)) - g(3) * (-t25 * mrSges(6,1) - t28 * mrSges(6,2)) * t8];
taug = t9(:);
