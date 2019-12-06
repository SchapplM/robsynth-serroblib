% Calculate Gravitation load on the joints for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:46
% EndTime: 2019-12-05 15:21:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (166->50), mult. (158->62), div. (0->0), fcn. (136->8), ass. (0->22)
t31 = m(5) + m(6);
t25 = m(4) + t31;
t15 = sin(pkin(9));
t17 = cos(pkin(9));
t33 = -mrSges(5,2) * t17 + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t15;
t16 = sin(pkin(8));
t18 = cos(pkin(8));
t32 = -mrSges(3,1) + (-m(5) * pkin(3) - t17 * mrSges(5,1) + t15 * mrSges(5,2) - mrSges(4,1)) * t18 + (mrSges(4,2) - mrSges(6,3)) * t16;
t14 = pkin(7) + qJ(2);
t10 = sin(t14);
t28 = t10 * t18;
t12 = cos(t14);
t27 = t12 * t18;
t23 = -t16 * (-pkin(6) - qJ(4)) + t18 * (pkin(4) * t17 + pkin(3));
t13 = pkin(9) + qJ(5);
t11 = cos(t13);
t9 = sin(t13);
t4 = t10 * t9 + t11 * t27;
t3 = t10 * t11 - t9 * t27;
t2 = -t11 * t28 + t12 * t9;
t1 = t11 * t12 + t9 * t28;
t5 = [(-m(2) - m(3) - t25) * g(3), (-t4 * mrSges(6,1) - t3 * mrSges(6,2) - t25 * (t12 * pkin(2) + t10 * qJ(3)) + t33 * t10 + (-(m(5) * qJ(4) + mrSges(5,3)) * t16 - m(6) * t23 + t32) * t12) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (m(4) * pkin(2) - m(5) * (-qJ(4) * t16 - pkin(2)) + t16 * mrSges(5,3) - m(6) * (-pkin(2) - t23) - t32) * t10 + (-qJ(3) * t25 + t33) * t12) * g(1), (-g(1) * t10 + g(2) * t12) * t25, (g(3) * t18 + t16 * (-g(1) * t12 - g(2) * t10)) * t31, -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t9 - mrSges(6,2) * t11) * t16];
taug = t5(:);
