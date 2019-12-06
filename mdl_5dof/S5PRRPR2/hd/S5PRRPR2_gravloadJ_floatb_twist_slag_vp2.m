% Calculate Gravitation load on the joints for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:18
% EndTime: 2019-12-05 16:17:19
% DurationCPUTime: 0.24s
% Computational Cost: add. (250->52), mult. (170->62), div. (0->0), fcn. (146->8), ass. (0->29)
t55 = -mrSges(5,3) + mrSges(4,2);
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t54 = -t28 * mrSges(5,1) - mrSges(4,1) + (mrSges(5,2) - mrSges(6,3)) * t27;
t52 = pkin(4) * t28 + pkin(7) * t27;
t51 = m(5) + m(6);
t26 = pkin(8) + qJ(2);
t25 = qJ(3) + t26;
t21 = sin(t25);
t22 = cos(t25);
t30 = cos(qJ(5));
t29 = sin(qJ(5));
t42 = t28 * t29;
t5 = t21 * t42 + t22 * t30;
t41 = t28 * t30;
t6 = -t21 * t41 + t22 * t29;
t50 = -t6 * mrSges(6,1) - t5 * mrSges(6,2) + t55 * t22 + (-m(6) * (-pkin(3) - t52) - t54) * t21;
t7 = t21 * t30 - t22 * t42;
t8 = t21 * t29 + t22 * t41;
t49 = -t8 * mrSges(6,1) - t7 * mrSges(6,2) + t55 * t21 + t54 * t22;
t23 = sin(t26);
t48 = pkin(2) * t23;
t24 = cos(t26);
t20 = pkin(2) * t24;
t40 = t22 * pkin(3) + t21 * qJ(4);
t15 = t22 * qJ(4);
t38 = -pkin(3) * t21 + t15;
t36 = t52 * t22 + t40;
t1 = [(-m(2) - m(3) - m(4) - t51) * g(3), (-mrSges(3,1) * t24 + mrSges(3,2) * t23 - m(4) * t20 - m(5) * (t20 + t40) - m(6) * (t20 + t36) + t49) * g(2) + (mrSges(3,1) * t23 + mrSges(3,2) * t24 + m(4) * t48 - m(5) * (t38 - t48) - m(6) * (t15 - t48) + t50) * g(1), (-m(5) * t40 - m(6) * t36 + t49) * g(2) + (-m(5) * t38 - m(6) * t15 + t50) * g(1), t51 * (-g(1) * t21 + g(2) * t22), -g(1) * (mrSges(6,1) * t7 - mrSges(6,2) * t8) - g(2) * (-mrSges(6,1) * t5 + mrSges(6,2) * t6) - g(3) * (-mrSges(6,1) * t29 - mrSges(6,2) * t30) * t27];
taug = t1(:);
