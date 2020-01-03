% Calculate Gravitation load on the joints for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:36
% EndTime: 2019-12-31 19:31:38
% DurationCPUTime: 0.53s
% Computational Cost: add. (248->77), mult. (298->83), div. (0->0), fcn. (256->10), ass. (0->43)
t15 = sin(pkin(9));
t16 = cos(pkin(9));
t63 = -m(5) * pkin(3) - t16 * mrSges(5,1) + t15 * mrSges(5,2) - mrSges(4,1);
t51 = m(5) + m(6);
t62 = m(4) + t51;
t61 = mrSges(4,2) - mrSges(6,3);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t57 = g(1) * t22 + g(2) * t20;
t14 = qJ(2) + pkin(8);
t9 = sin(t14);
t60 = t9 * t57;
t11 = cos(t14);
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t59 = -t21 * mrSges(3,1) + t19 * mrSges(3,2) + t63 * t11 + t61 * t9;
t58 = m(5) * qJ(4) + mrSges(5,3);
t56 = -m(3) * pkin(1) - mrSges(2,1) + t59;
t17 = -qJ(3) - pkin(6);
t53 = -m(3) * pkin(6) + m(5) * t17 - t15 * mrSges(5,1) - t16 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t49 = pkin(4) * t15;
t13 = pkin(9) + qJ(5);
t8 = sin(t13);
t46 = t20 * t8;
t12 = t21 * pkin(2);
t45 = t22 * t8;
t44 = t9 * mrSges(5,3);
t18 = -pkin(7) - qJ(4);
t42 = t9 * t18;
t10 = cos(t13);
t41 = t20 * t10;
t40 = t22 * t10;
t39 = t9 * qJ(4);
t6 = t16 * pkin(4) + pkin(3);
t34 = t11 * t6 - t42;
t30 = m(6) * t6 + t10 * mrSges(6,1) - t8 * mrSges(6,2);
t7 = t12 + pkin(1);
t5 = t22 * t7;
t4 = t11 * t40 + t46;
t3 = -t11 * t45 + t41;
t2 = -t11 * t41 + t45;
t1 = t11 * t46 + t40;
t23 = [(-m(5) * t5 - t4 * mrSges(6,1) - t3 * mrSges(6,2) + (-m(4) - m(6)) * (-t20 * t17 + t5) + (-m(6) * t49 + t53) * t20 + (-m(6) * t34 - t58 * t9 + t56) * t22) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (m(4) * t17 - m(6) * (-t17 + t49) + t53) * t22 + (m(4) * t7 - m(5) * (-t7 - t39) + t44 - m(6) * (-t34 - t7) - t56) * t20) * g(1), (t30 - t63) * t60 + (-m(4) * t12 - m(5) * (t12 + t39) - t44 - m(6) * (t12 - t42) + t59 - t30 * t11) * g(3) + (mrSges(3,2) * t21 + (m(6) * t18 - t58 + t61) * t11 + (t62 * pkin(2) + mrSges(3,1)) * t19) * t57, (-g(1) * t20 + g(2) * t22) * t62, (t11 * g(3) - t60) * t51, -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t8 - mrSges(6,2) * t10) * t9];
taug = t23(:);
