% Calculate Gravitation load on the joints for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:48
% EndTime: 2019-03-09 02:02:49
% DurationCPUTime: 0.58s
% Computational Cost: add. (306->74), mult. (354->93), div. (0->0), fcn. (318->8), ass. (0->43)
t61 = mrSges(6,1) + mrSges(7,1);
t60 = -mrSges(6,2) + mrSges(7,3);
t56 = -mrSges(6,3) - mrSges(7,2);
t19 = sin(qJ(5));
t22 = cos(qJ(5));
t59 = t60 * t19 + t61 * t22;
t58 = -m(4) - m(5);
t57 = m(6) + m(7);
t29 = pkin(5) * t22 + qJ(6) * t19;
t55 = m(6) * pkin(4) - m(7) * (-pkin(4) - t29) + t59;
t54 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t20 = sin(qJ(4));
t23 = cos(qJ(4));
t53 = -t20 * mrSges(5,1) + (-mrSges(5,2) - t56) * t23;
t52 = m(7) * pkin(5) + t61;
t51 = m(7) * qJ(6) + t60;
t50 = mrSges(3,2) - mrSges(4,3) + t53;
t49 = -pkin(2) - pkin(7);
t18 = qJ(1) + pkin(9);
t14 = sin(t18);
t47 = g(1) * t14;
t15 = cos(t18);
t46 = g(2) * t15;
t45 = g(3) * t23;
t21 = sin(qJ(1));
t44 = t21 * pkin(1);
t16 = t23 * pkin(8);
t24 = cos(qJ(1));
t17 = t24 * pkin(1);
t43 = t14 * t20;
t42 = t14 * t23;
t41 = t20 * t19;
t40 = t20 * t22;
t37 = t15 * pkin(2) + t14 * qJ(3) + t17;
t36 = t57 - t58;
t35 = t15 * qJ(3) - t44;
t34 = t15 * pkin(7) + t37;
t33 = mrSges(5,1) * t23 - mrSges(5,2) * t20;
t4 = -t14 * t19 + t15 * t40;
t3 = t14 * t22 + t15 * t41;
t2 = t14 * t40 + t15 * t19;
t1 = t14 * t41 - t15 * t22;
t5 = [(-m(3) * t17 - m(4) * t37 - m(5) * t34 - t24 * mrSges(2,1) + t21 * mrSges(2,2) - t57 * (pkin(4) * t43 - pkin(8) * t42 + t34) - t52 * t2 - t51 * t1 + t54 * t15 + t50 * t14) * g(2) + (m(3) * t44 + t21 * mrSges(2,1) + t24 * mrSges(2,2) - t57 * (t49 * t14 + t35) - t52 * t4 + t58 * t35 - t51 * t3 + (m(4) * pkin(2) - m(5) * t49 - t54) * t14 + (-t57 * (t20 * pkin(4) - t16) + t50) * t15) * g(1) (-m(3) - t36) * g(3) (-t47 + t46) * t36, -t33 * t47 + (-t57 * (pkin(4) * t42 + pkin(8) * t43) + ((-m(7) * t29 - t59) * t23 + t56 * t20) * t14) * g(1) + (t33 + t55 * t23 + (t57 * pkin(8) - t56) * t20) * t46 + (-t57 * t16 + t55 * t20 - t53) * g(3) (t52 * t19 - t51 * t22) * t45 + (-t52 * t3 + t51 * t4) * g(2) + (t52 * t1 - t51 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - t19 * t45) * m(7)];
taug  = t5(:);
