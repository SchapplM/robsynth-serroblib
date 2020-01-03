% Calculate Gravitation load on the joints for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:33
% DurationCPUTime: 0.46s
% Computational Cost: add. (152->67), mult. (331->85), div. (0->0), fcn. (304->6), ass. (0->41)
t59 = mrSges(5,1) + mrSges(6,1);
t58 = -mrSges(5,2) + mrSges(6,3);
t52 = -mrSges(5,3) - mrSges(6,2);
t53 = m(5) + m(6);
t57 = t53 * pkin(7);
t56 = mrSges(2,2) - mrSges(3,3);
t15 = sin(qJ(4));
t18 = cos(qJ(4));
t55 = t58 * t15 + t59 * t18;
t54 = -m(3) - m(4);
t25 = pkin(4) * t18 + qJ(5) * t15;
t51 = m(5) * pkin(3) - m(6) * (-pkin(3) - t25) + t55;
t50 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t19 = cos(qJ(3));
t16 = sin(qJ(3));
t28 = t16 * mrSges(4,1) + t19 * mrSges(4,2);
t49 = -t52 * t19 - t28;
t48 = m(6) * pkin(4) + t59;
t47 = m(6) * qJ(5) + t58;
t46 = -pkin(1) - pkin(6);
t17 = sin(qJ(1));
t44 = g(1) * t17;
t20 = cos(qJ(1));
t43 = g(2) * t20;
t42 = g(3) * t19;
t41 = t16 * t17;
t40 = t17 * t15;
t39 = t17 * t18;
t38 = t17 * t19;
t35 = t19 * t20;
t34 = t20 * t15;
t33 = t20 * t18;
t32 = t20 * pkin(1) + t17 * qJ(2);
t31 = t20 * pkin(6) + t32;
t29 = mrSges(4,1) * t19 - mrSges(4,2) * t16;
t11 = t20 * qJ(2);
t4 = t16 * t33 - t40;
t3 = t16 * t34 + t39;
t2 = t16 * t39 + t34;
t1 = t16 * t40 - t33;
t5 = [(-m(3) * t32 - m(4) * t31 - t53 * (pkin(3) * t41 - pkin(7) * t38 + t31) - t48 * t2 - t47 * t1 + t50 * t20 + (t49 + t56) * t17) * g(2) + (-t53 * (t20 * t16 * pkin(3) - pkin(7) * t35 + t46 * t17 + t11) - t48 * t4 - t52 * t35 - t47 * t3 + t54 * t11 + (-t28 + t56) * t20 + (m(3) * pkin(1) - m(4) * t46 - t50) * t17) * g(1), (-t44 + t43) * (t53 - t54), -t29 * t44 + (-t53 * (pkin(3) * t38 + pkin(7) * t41) + ((-m(6) * t25 - t55) * t19 + t52 * t16) * t17) * g(1) + (t29 + t51 * t19 + (-t52 + t57) * t16) * t43 + (t51 * t16 - t19 * t57 - t49) * g(3), (t48 * t15 - t47 * t18) * t42 + (-t48 * t3 + t47 * t4) * g(2) + (t48 * t1 - t47 * t2) * g(1), (-g(1) * t1 + g(2) * t3 - t15 * t42) * m(6)];
taug = t5(:);
