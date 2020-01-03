% Calculate Gravitation load on the joints for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:47
% EndTime: 2019-12-31 20:35:49
% DurationCPUTime: 0.76s
% Computational Cost: add. (448->97), mult. (807->137), div. (0->0), fcn. (920->12), ass. (0->44)
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t77 = m(6) * pkin(4) + t45 * mrSges(6,1) - t42 * mrSges(6,2) + mrSges(5,1);
t58 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t52 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t75 = t42 * mrSges(6,1) + t45 * mrSges(6,2) - t52;
t37 = pkin(10) + qJ(4);
t34 = sin(t37);
t35 = cos(t37);
t38 = sin(pkin(10));
t40 = cos(pkin(10));
t78 = -m(4) * pkin(2) - t40 * mrSges(4,1) + t38 * mrSges(4,2) - mrSges(3,1);
t80 = t58 * t34 - t77 * t35 + t78;
t79 = m(5) + m(6);
t70 = cos(qJ(1));
t39 = sin(pkin(5));
t43 = sin(qJ(2));
t69 = t39 * t43;
t44 = sin(qJ(1));
t68 = t39 * t44;
t46 = cos(qJ(2));
t67 = t39 * t46;
t65 = t70 * pkin(1) + pkin(7) * t68;
t64 = cos(pkin(5));
t62 = t39 * t70;
t61 = -pkin(1) * t44 + pkin(7) * t62;
t55 = t64 * t70;
t20 = t43 * t55 + t44 * t46;
t4 = t20 * t35 - t34 * t62;
t3 = -t20 * t34 - t35 * t62;
t59 = t44 * t64;
t57 = t38 * t62;
t21 = t70 * t43 + t46 * t59;
t22 = -t43 * t59 + t70 * t46;
t33 = pkin(3) * t40 + pkin(2);
t41 = -pkin(8) - qJ(3);
t56 = t38 * pkin(3) * t68 - t21 * t41 + t22 * t33 + t65;
t19 = t43 * t44 - t46 * t55;
t14 = t64 * t34 + t35 * t69;
t8 = t22 * t35 + t34 * t68;
t7 = t22 * t34 - t35 * t68;
t2 = t21 * t42 + t45 * t8;
t1 = t21 * t45 - t42 * t8;
t5 = [(-t70 * mrSges(2,1) - m(5) * t56 - t8 * mrSges(5,1) - m(6) * (pkin(4) * t8 + t56) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t58 * t7 + t78 * t22 + (mrSges(2,2) + (-mrSges(4,1) * t38 - mrSges(4,2) * t40 - mrSges(3,3)) * t39) * t44 + t52 * t21 + (-m(3) - m(4)) * t65) * g(2) + (t44 * mrSges(2,1) + t70 * mrSges(2,2) - m(3) * t61 + t20 * mrSges(3,1) - mrSges(3,3) * t62 - m(4) * (-pkin(2) * t20 + t61) - (-t20 * t40 + t57) * mrSges(4,1) - (t20 * t38 + t40 * t62) * mrSges(4,2) + t58 * t3 + t77 * t4 + t75 * t19 + t79 * (-pkin(3) * t57 - t19 * t41 + t20 * t33 - t61)) * g(1), (-t79 * (-t19 * t33 - t20 * t41) - t75 * t20 - t80 * t19) * g(2) + (-t79 * (-t21 * t33 - t22 * t41) - t75 * t22 - t80 * t21) * g(1) + (-t79 * t33 * t67 + (t80 * t46 + (t79 * t41 - t75) * t43) * t39) * g(3), (m(4) + t79) * (-g(1) * t21 - g(2) * t19 + g(3) * t67), (t58 * t14 - t77 * (-t34 * t69 + t64 * t35)) * g(3) + (-t77 * t3 + t58 * t4) * g(2) + (t58 * t8 + t77 * t7) * g(1), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((t19 * t45 - t4 * t42) * mrSges(6,1) + (-t19 * t42 - t4 * t45) * mrSges(6,2)) - g(3) * ((-t14 * t42 - t45 * t67) * mrSges(6,1) + (-t14 * t45 + t42 * t67) * mrSges(6,2))];
taug = t5(:);
