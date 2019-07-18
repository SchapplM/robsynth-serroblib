% Calculate Gravitation load on the joints for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:44
% EndTime: 2019-07-18 17:20:47
% DurationCPUTime: 0.40s
% Computational Cost: add. (182->68), mult. (264->82), div. (0->0), fcn. (217->8), ass. (0->41)
t65 = m(4) * pkin(1) + mrSges(3,1) + mrSges(4,1);
t64 = mrSges(3,2) + mrSges(4,2);
t15 = qJ(2) + qJ(4);
t14 = cos(t15);
t12 = t14 * mrSges(5,1);
t13 = sin(t15);
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t63 = t13 * mrSges(5,2) + t64 * t18 - t65 * t21 - t12;
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t62 = g(1) * t22 + g(2) * t19;
t23 = pkin(2) + pkin(1);
t43 = t21 * t23;
t51 = t13 * pkin(4);
t61 = -m(5) * t43 - m(6) * (t43 + t51) + t63;
t59 = m(5) + m(6);
t57 = -m(4) * qJ(3) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t17 = sin(qJ(5));
t49 = mrSges(6,2) * t17;
t38 = t13 * t49;
t55 = (mrSges(6,3) * t14 + t38) * t19;
t42 = t22 * mrSges(6,3);
t54 = t14 * t42 + t22 * t38;
t20 = cos(qJ(5));
t50 = mrSges(6,1) * t20;
t11 = t13 * mrSges(6,3);
t46 = t18 * t23;
t45 = t19 * t17;
t44 = t19 * t20;
t41 = t22 * t17;
t40 = t22 * t20;
t36 = -m(6) * pkin(4) + mrSges(5,2);
t30 = t11 + (-t49 + t50) * t14;
t25 = m(6) * (pkin(4) * t14 - t46) - t13 * t50;
t16 = -pkin(3) - qJ(3);
t4 = t14 * t40 + t45;
t3 = -t14 * t41 + t44;
t2 = -t14 * t44 + t41;
t1 = t14 * t45 + t40;
t5 = [(-t4 * mrSges(6,1) - t3 * mrSges(6,2) - t13 * t42 - t59 * (-t19 * t16 + t22 * t43) + (-m(6) * t51 - mrSges(2,1) + t63) * t22 + t57 * t19) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (t59 * t16 + t57) * t22 + (mrSges(2,1) + t11 - t61) * t19) * g(1), -g(1) * (t25 * t22 + t54) - g(2) * (t25 * t19 + t55) + (-t30 + t61) * g(3) + t62 * (m(5) * t46 + mrSges(5,1) * t13 + mrSges(5,2) * t14 + t65 * t18 + t64 * t21), (-g(1) * t19 + g(2) * t22) * (m(4) + t59), -g(1) * t54 - g(2) * t55 + (t36 * t13 - t12 - t30) * g(3) + t62 * (t36 * t14 + (mrSges(5,1) + t50) * t13), -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t17 - mrSges(6,2) * t20) * t13];
taug  = t5(:);
