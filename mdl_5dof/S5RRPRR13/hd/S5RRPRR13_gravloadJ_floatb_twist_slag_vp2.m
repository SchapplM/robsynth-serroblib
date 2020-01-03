% Calculate Gravitation load on the joints for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:09
% EndTime: 2019-12-31 20:32:11
% DurationCPUTime: 0.57s
% Computational Cost: add. (298->91), mult. (389->105), div. (0->0), fcn. (353->10), ass. (0->52)
t77 = mrSges(5,3) + mrSges(6,3);
t28 = cos(pkin(9));
t17 = t28 * pkin(3) + pkin(2);
t26 = pkin(9) + qJ(4);
t19 = cos(t26);
t13 = pkin(4) * t19 + t17;
t20 = qJ(5) + t26;
t15 = sin(t20);
t16 = cos(t20);
t18 = sin(t26);
t76 = -m(5) * t17 - m(6) * t13 - t19 * mrSges(5,1) - t16 * mrSges(6,1) + t18 * mrSges(5,2) + t15 * mrSges(6,2);
t29 = -pkin(7) - qJ(3);
t25 = -pkin(8) + t29;
t75 = m(5) * t29 + m(6) * t25 - t77;
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t74 = g(1) * t33 + g(2) * t31;
t32 = cos(qJ(2));
t27 = sin(pkin(9));
t40 = m(4) * pkin(2) + t28 * mrSges(4,1) - t27 * mrSges(4,2);
t73 = t40 * t32;
t72 = m(6) * pkin(4) + mrSges(5,1);
t71 = m(4) + m(5) + m(6);
t70 = -m(3) - t71;
t30 = sin(qJ(2));
t46 = t32 * mrSges(3,1) - t30 * mrSges(3,2);
t69 = t77 * t30 + mrSges(2,1) + t46;
t59 = t27 * pkin(3);
t63 = pkin(4) * t18;
t67 = -m(5) * t59 - m(6) * (t59 + t63) + mrSges(2,2) - mrSges(3,3) - t27 * mrSges(4,1) - t28 * mrSges(4,2);
t54 = t33 * t16;
t56 = t31 * t32;
t5 = t15 * t56 + t54;
t55 = t33 * t15;
t6 = -t16 * t56 + t55;
t65 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = t31 * t16 - t32 * t55;
t8 = t31 * t15 + t32 * t54;
t64 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t60 = g(3) * t30;
t53 = t33 * t18;
t52 = t33 * t19;
t47 = m(4) * qJ(3) + mrSges(4,3);
t43 = -mrSges(6,1) * t15 - mrSges(6,2) * t16;
t42 = t32 * t13 - t30 * t25;
t41 = t32 * t17 - t30 * t29;
t11 = t31 * t19 - t32 * t53;
t9 = t18 * t56 + t52;
t37 = t47 * t30 + t73;
t12 = t31 * t18 + t32 * t52;
t10 = -t19 * t56 + t53;
t1 = [(-t12 * mrSges(5,1) - t8 * mrSges(6,1) - t11 * mrSges(5,2) - t7 * mrSges(6,2) + t70 * (t33 * pkin(1) + t31 * pkin(6)) + t67 * t31 + (-m(5) * t41 - m(6) * t42 - t37 - t69) * t33) * g(2) + (-t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (m(3) * pkin(1) - m(4) * (-t30 * qJ(3) - pkin(1)) + t30 * mrSges(4,3) + t73 - m(5) * (-pkin(1) - t41) - m(6) * (-pkin(1) - t42) + t69) * t31 + (t70 * pkin(6) + t67) * t33) * g(1), (-t37 - t46) * g(3) + (t76 * g(3) + t74 * (mrSges(3,2) - t47 + t75)) * t32 + (t75 * g(3) + t74 * (mrSges(3,1) + t40 - t76)) * t30, (t32 * g(3) - t74 * t30) * t71, (m(6) * t63 + mrSges(5,1) * t18 + mrSges(5,2) * t19 - t43) * t60 + (-t10 * mrSges(5,2) + t72 * t9 - t65) * g(2) + (t12 * mrSges(5,2) - t72 * t11 - t64) * g(1), -g(1) * t64 - g(2) * t65 - t43 * t60];
taug = t1(:);
