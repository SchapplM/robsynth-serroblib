% Calculate Gravitation load on the joints for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:46
% EndTime: 2019-12-05 15:44:47
% DurationCPUTime: 0.29s
% Computational Cost: add. (193->50), mult. (193->64), div. (0->0), fcn. (153->10), ass. (0->31)
t20 = qJ(2) + pkin(9);
t18 = qJ(4) + t20;
t14 = sin(t18);
t15 = cos(t18);
t25 = cos(qJ(5));
t58 = -t25 * mrSges(6,1) - mrSges(5,1);
t27 = mrSges(5,2) * t15 + (m(6) * pkin(4) - t58) * t14;
t59 = -m(5) - m(6);
t23 = sin(qJ(5));
t41 = t23 * mrSges(6,2);
t57 = -t14 * t41 + t15 * (-m(6) * pkin(7) - mrSges(6,3));
t16 = sin(t20);
t17 = cos(t20);
t24 = sin(qJ(2));
t26 = cos(qJ(2));
t51 = pkin(2) * t24;
t56 = m(4) * t51 + mrSges(3,1) * t24 + mrSges(4,1) * t16 + mrSges(3,2) * t26 + mrSges(4,2) * t17 + t59 * (-pkin(3) * t16 - t51) + t27;
t53 = (t41 + t58) * t15 + (-mrSges(6,3) + mrSges(5,2)) * t14;
t19 = t26 * pkin(2);
t21 = sin(pkin(8));
t45 = t21 * t23;
t44 = t21 * t25;
t22 = cos(pkin(8));
t43 = t22 * t23;
t42 = t22 * t25;
t39 = t15 * pkin(4) + t14 * pkin(7);
t38 = pkin(3) * t17 + t19;
t37 = m(4) - t59;
t34 = t57 * t21;
t33 = t57 * t22;
t1 = [(-m(2) - m(3) - t37) * g(3), (-mrSges(3,1) * t26 + t24 * mrSges(3,2) - m(4) * t19 - t17 * mrSges(4,1) + t16 * mrSges(4,2) - m(5) * t38 - m(6) * (t38 + t39) + t53) * g(3) + (t56 * t21 + t34) * g(2) + (t56 * t22 + t33) * g(1), (-g(1) * t21 + g(2) * t22) * t37, (-m(6) * t39 + t53) * g(3) + (t27 * t21 + t34) * g(2) + (t27 * t22 + t33) * g(1), -g(1) * ((-t15 * t43 + t44) * mrSges(6,1) + (-t15 * t42 - t45) * mrSges(6,2)) - g(2) * ((-t15 * t45 - t42) * mrSges(6,1) + (-t15 * t44 + t43) * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t23 - mrSges(6,2) * t25) * t14];
taug = t1(:);
