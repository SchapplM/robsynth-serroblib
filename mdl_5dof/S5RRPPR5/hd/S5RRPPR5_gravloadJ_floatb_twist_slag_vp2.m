% Calculate Gravitation load on the joints for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:28:55
% DurationCPUTime: 0.42s
% Computational Cost: add. (228->75), mult. (292->83), div. (0->0), fcn. (254->8), ass. (0->39)
t19 = qJ(2) + pkin(8);
t16 = sin(t19);
t17 = cos(t19);
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t62 = -t25 * mrSges(3,1) + t22 * mrSges(3,2) - (mrSges(4,1) + mrSges(5,1)) * t17 - (-mrSges(4,2) + mrSges(5,3)) * t16;
t51 = pkin(2) * t22;
t52 = -pkin(3) - pkin(4);
t53 = m(5) + m(6);
t61 = -m(5) * (-pkin(3) * t16 - t51) + t16 * mrSges(5,1) - m(6) * (t16 * t52 - t51) + (-qJ(4) * t53 - mrSges(5,3)) * t17;
t13 = t16 * qJ(4);
t26 = cos(qJ(1));
t47 = t26 * t17;
t60 = pkin(3) * t47 + t26 * t13;
t23 = sin(qJ(1));
t57 = g(1) * t26 + g(2) * t23;
t55 = -m(3) * pkin(1) - mrSges(2,1) + t62;
t20 = -qJ(3) - pkin(6);
t54 = mrSges(2,2) - m(3) * pkin(6) - mrSges(3,3) - mrSges(4,3) - mrSges(5,2) - m(6) * (-pkin(7) - t20) + mrSges(6,3);
t14 = t17 * pkin(3);
t18 = t25 * pkin(2);
t24 = cos(qJ(5));
t48 = t16 * t24;
t44 = t14 + t13 + t18;
t15 = t18 + pkin(1);
t12 = t26 * t15;
t42 = -t23 * t20 + t12;
t41 = -t15 - t13;
t21 = sin(qJ(5));
t6 = -t17 * t21 + t48;
t1 = t6 * t23;
t31 = t16 * t21 + t17 * t24;
t2 = t31 * t23;
t39 = t1 * mrSges(6,1) - t2 * mrSges(6,2);
t3 = t21 * t47 - t26 * t48;
t4 = t31 * t26;
t38 = -t3 * mrSges(6,1) - t4 * mrSges(6,2);
t37 = -mrSges(6,1) * t31 - t6 * mrSges(6,2);
t5 = [(-m(4) * t42 - m(5) * (t42 + t60) - m(6) * (pkin(4) * t47 + t12 + t60) - t4 * mrSges(6,1) + t3 * mrSges(6,2) + t55 * t26 + t54 * t23) * g(2) + (t2 * mrSges(6,1) + t1 * mrSges(6,2) + ((m(4) + m(5)) * t20 + t54) * t26 + (m(4) * t15 - m(5) * (t41 - t14) - m(6) * (t17 * t52 + t41) - t55) * t23) * g(1), (t61 * t23 + t39) * g(2) + (t61 * t26 + t38) * g(1) + (-m(4) * t18 - m(5) * t44 - m(6) * (t17 * pkin(4) + t44) + t37 + t62) * g(3) + (m(4) * t51 + mrSges(3,1) * t22 + mrSges(4,1) * t16 + mrSges(3,2) * t25 + mrSges(4,2) * t17) * t57, (-g(1) * t23 + g(2) * t26) * (m(4) + t53), (g(3) * t17 - t16 * t57) * t53, -g(1) * t38 - g(2) * t39 - g(3) * t37];
taug = t5(:);
