% Calculate Gravitation load on the joints for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:11
% EndTime: 2019-12-05 15:14:12
% DurationCPUTime: 0.26s
% Computational Cost: add. (168->43), mult. (195->57), div. (0->0), fcn. (167->8), ass. (0->27)
t14 = sin(qJ(4));
t15 = cos(qJ(4));
t11 = qJ(4) + qJ(5);
t8 = sin(t11);
t9 = cos(t11);
t45 = m(5) * pkin(3) + t15 * mrSges(5,1) - t14 * mrSges(5,2) + mrSges(4,1) + m(6) * (pkin(4) * t15 + pkin(3)) + t9 * mrSges(6,1) - t8 * mrSges(6,2);
t44 = mrSges(4,2) + m(6) * (-pkin(7) - pkin(6)) - mrSges(6,3) - m(5) * pkin(6) - mrSges(5,3);
t43 = m(6) * pkin(4) + mrSges(5,1);
t10 = pkin(9) + qJ(3);
t6 = sin(t10);
t40 = g(3) * t6;
t13 = cos(pkin(8));
t32 = t13 * t9;
t33 = t13 * t8;
t12 = sin(pkin(8));
t34 = t12 * t9;
t35 = t12 * t8;
t7 = cos(t10);
t39 = (-t7 * t35 - t32) * mrSges(6,1) + (-t7 * t34 + t33) * mrSges(6,2);
t38 = (-t7 * t33 + t34) * mrSges(6,1) + (-t7 * t32 - t35) * mrSges(6,2);
t31 = t12 * t14;
t30 = t12 * t15;
t29 = t13 * t14;
t28 = t13 * t15;
t27 = m(3) + m(4) + m(5) + m(6);
t23 = -mrSges(6,1) * t8 - mrSges(6,2) * t9;
t1 = [(-m(2) - t27) * g(3), (-g(1) * t12 + g(2) * t13) * t27, (t44 * t6 - t45 * t7) * g(3) + (g(1) * t13 + g(2) * t12) * (t44 * t7 + t45 * t6), (mrSges(5,2) * t15 + t43 * t14 - t23) * t40 + (-(-t7 * t30 + t29) * mrSges(5,2) - t39 - t43 * (-t7 * t31 - t28)) * g(2) + (-(-t7 * t28 - t31) * mrSges(5,2) - t38 - t43 * (-t7 * t29 + t30)) * g(1), -g(1) * t38 - g(2) * t39 - t23 * t40];
taug = t1(:);
