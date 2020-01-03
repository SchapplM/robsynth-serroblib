% Calculate Gravitation load on the joints for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:59:55
% EndTime: 2019-12-31 19:59:57
% DurationCPUTime: 0.60s
% Computational Cost: add. (264->73), mult. (364->88), div. (0->0), fcn. (333->8), ass. (0->44)
t70 = mrSges(5,1) + mrSges(6,1);
t69 = -mrSges(5,2) + mrSges(6,3);
t63 = m(5) + m(6);
t68 = mrSges(5,3) + mrSges(6,2);
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t67 = t69 * t20 + t70 * t23;
t22 = sin(qJ(1));
t66 = g(2) * t22;
t65 = m(4) + t63;
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t36 = t24 * mrSges(3,1) - t21 * mrSges(3,2);
t64 = -m(3) * pkin(1) - mrSges(2,1) - t36;
t18 = qJ(2) + pkin(8);
t15 = sin(t18);
t16 = cos(t18);
t62 = -t16 * pkin(3) - t15 * pkin(7);
t31 = pkin(4) * t23 + qJ(5) * t20;
t54 = pkin(2) * t21;
t61 = t63 * t54 - t68 * t16 + (-m(6) * (-pkin(3) - t31) + m(5) * pkin(3) + t67) * t15;
t59 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t34 = t16 * mrSges(4,1) - t15 * mrSges(4,2);
t58 = t68 * t15 + t34;
t57 = pkin(7) * t63;
t56 = m(6) * pkin(4) + t70;
t55 = m(6) * qJ(5) + t69;
t51 = g(3) * t15;
t17 = t24 * pkin(2);
t48 = t22 * t20;
t47 = t22 * t23;
t25 = cos(qJ(1));
t46 = t25 * t15;
t45 = t25 * t16;
t43 = t25 * t20;
t42 = t25 * t23;
t14 = t17 + pkin(1);
t19 = -qJ(3) - pkin(6);
t38 = t25 * t14 - t22 * t19;
t4 = t16 * t42 + t48;
t3 = t16 * t43 - t47;
t2 = t16 * t47 - t43;
t1 = t16 * t48 + t42;
t5 = [(-m(4) * t38 - t68 * t46 - t63 * (pkin(3) * t45 + pkin(7) * t46 + t38) - t56 * t4 - t55 * t3 + (-t34 + t64) * t25 + t59 * t22) * g(2) + (t56 * t2 + t55 * t1 + (m(4) * t14 - t63 * (-t14 + t62) + t58 - t64) * t22 + (t65 * t19 + t59) * t25) * g(1), (t61 * t25 - t45 * t57) * g(1) + (-m(4) * t17 - t36 - t63 * (t17 - t62) + (-m(6) * t31 - t67) * t16 - t58) * g(3) + (m(4) * t54 + mrSges(3,1) * t21 + mrSges(4,1) * t15 + mrSges(3,2) * t24 + mrSges(4,2) * t16) * (g(1) * t25 + t66) + (-t16 * t57 + t61) * t66, (-g(1) * t22 + g(2) * t25) * t65, (t56 * t20 - t55 * t23) * t51 + (t56 * t1 - t55 * t2) * g(2) + (t56 * t3 - t55 * t4) * g(1), (-g(1) * t3 - g(2) * t1 - t20 * t51) * m(6)];
taug = t5(:);
