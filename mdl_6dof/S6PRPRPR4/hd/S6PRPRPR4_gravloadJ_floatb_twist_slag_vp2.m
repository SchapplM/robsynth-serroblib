% Calculate Gravitation load on the joints for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:05
% EndTime: 2019-03-08 19:39:06
% DurationCPUTime: 0.73s
% Computational Cost: add. (505->75), mult. (817->105), div. (0->0), fcn. (920->14), ass. (0->41)
t27 = pkin(12) + qJ(6);
t23 = sin(t27);
t25 = cos(t27);
t29 = sin(pkin(12));
t32 = cos(pkin(12));
t68 = mrSges(5,1) + m(7) * (pkin(5) * t32 + pkin(4)) + t25 * mrSges(7,1) - t23 * mrSges(7,2) + m(6) * pkin(4) + t32 * mrSges(6,1) - t29 * mrSges(6,2);
t67 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3);
t28 = pkin(11) + qJ(4);
t24 = sin(t28);
t26 = cos(t28);
t33 = cos(pkin(11));
t66 = m(4) * pkin(2) + t33 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2) + mrSges(3,1) + t68 * t26 - t67 * t24;
t65 = -m(4) * qJ(3) - t23 * mrSges(7,1) - t32 * mrSges(6,2) - t25 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t29;
t75 = m(6) + m(7);
t69 = m(5) + t75;
t31 = sin(pkin(6));
t36 = sin(qJ(2));
t60 = t31 * t36;
t37 = cos(qJ(2));
t59 = t31 * t37;
t58 = cos(pkin(6));
t57 = cos(pkin(10));
t56 = sin(pkin(10));
t55 = m(4) + t69;
t51 = t31 * t57;
t50 = t31 * t56;
t46 = t58 * t57;
t45 = t58 * t56;
t35 = -pkin(8) - qJ(3);
t22 = pkin(3) * t33 + pkin(2);
t14 = -t36 * t45 + t57 * t37;
t13 = t57 * t36 + t37 * t45;
t12 = t36 * t46 + t56 * t37;
t11 = t56 * t36 - t37 * t46;
t8 = t58 * t24 + t26 * t60;
t7 = t24 * t60 - t58 * t26;
t4 = t14 * t26 + t24 * t50;
t3 = t14 * t24 - t26 * t50;
t2 = t12 * t26 - t24 * t51;
t1 = t12 * t24 + t26 * t51;
t5 = [(-m(2) - m(3) - t55) * g(3) (-t69 * (-t11 * t22 - t12 * t35) + t65 * t12 + t66 * t11) * g(2) + (-t69 * (-t13 * t22 - t14 * t35) + t65 * t14 + t66 * t13) * g(1) + (-t69 * t22 * t59 + (-t66 * t37 + (t69 * t35 + t65) * t36) * t31) * g(3) (-g(1) * t13 - g(2) * t11 + g(3) * t59) * t55 (t67 * t8 + t68 * t7) * g(3) + (t68 * t1 + t67 * t2) * g(2) + (t68 * t3 + t67 * t4) * g(1), t75 * (-g(1) * t3 - g(2) * t1 - g(3) * t7) -g(1) * ((t13 * t25 - t23 * t4) * mrSges(7,1) + (-t13 * t23 - t25 * t4) * mrSges(7,2)) - g(2) * ((t11 * t25 - t2 * t23) * mrSges(7,1) + (-t11 * t23 - t2 * t25) * mrSges(7,2)) - g(3) * ((-t8 * t23 - t25 * t59) * mrSges(7,1) + (t23 * t59 - t8 * t25) * mrSges(7,2))];
taug  = t5(:);
