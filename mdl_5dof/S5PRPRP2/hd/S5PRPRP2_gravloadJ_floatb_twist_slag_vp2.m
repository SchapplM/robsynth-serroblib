% Calculate Gravitation load on the joints for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:32
% DurationCPUTime: 0.28s
% Computational Cost: add. (170->47), mult. (188->58), div. (0->0), fcn. (167->6), ass. (0->25)
t29 = m(6) * pkin(4);
t34 = -mrSges(5,1) - mrSges(6,1);
t33 = mrSges(3,2) - mrSges(4,3);
t32 = mrSges(5,2) + mrSges(6,2);
t20 = m(4) + m(5) + m(6);
t31 = t29 - t34;
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t30 = t13 * mrSges(4,1) + mrSges(3,1) + (-mrSges(4,2) + mrSges(5,3) + mrSges(6,3)) * t12;
t15 = sin(qJ(4));
t11 = pkin(7) + qJ(2);
t9 = sin(t11);
t26 = t9 * t15;
t10 = cos(t11);
t25 = t10 * t15;
t22 = t13 * t15;
t16 = cos(qJ(4));
t21 = t13 * t16;
t19 = pkin(3) * t13 + pkin(6) * t12;
t18 = -t12 * (-qJ(5) - pkin(6)) + t13 * (pkin(4) * t16 + pkin(3));
t1 = t10 * t16 + t9 * t22;
t3 = -t10 * t22 + t16 * t9;
t4 = t10 * t21 + t26;
t2 = -t9 * t21 + t25;
t5 = [(-m(2) - m(3) - t20) * g(3), (-t26 * t29 + t33 * t9 + t34 * t4 - t32 * t3 - t20 * (t10 * pkin(2) + t9 * qJ(3)) + (-m(5) * t19 - m(6) * t18 - t30) * t10) * g(2) + (-t25 * t29 + t34 * t2 - t32 * t1 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t19) - m(6) * (-pkin(2) - t18) + t30) * t9 + (-t20 * qJ(3) + t33) * t10) * g(1), (-g(1) * t9 + g(2) * t10) * t20, (t31 * t15 + t32 * t16) * g(3) * t12 + (t31 * t1 - t32 * t2) * g(2) + (-t31 * t3 + t32 * t4) * g(1), (g(3) * t13 + (-g(1) * t10 - g(2) * t9) * t12) * m(6)];
taug = t5(:);
