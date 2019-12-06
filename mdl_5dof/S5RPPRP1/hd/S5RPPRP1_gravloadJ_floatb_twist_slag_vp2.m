% Calculate Gravitation load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:36
% EndTime: 2019-12-05 17:35:38
% DurationCPUTime: 0.32s
% Computational Cost: add. (181->51), mult. (203->64), div. (0->0), fcn. (179->8), ass. (0->27)
t32 = m(6) * pkin(4);
t37 = mrSges(5,1) + mrSges(6,1);
t36 = mrSges(3,2) - mrSges(4,3);
t35 = mrSges(5,2) + mrSges(6,2);
t24 = -m(4) - m(5) - m(6);
t34 = t32 + t37;
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t17 = cos(qJ(4));
t33 = mrSges(3,1) + m(4) * pkin(2) + mrSges(4,1) * t13 - m(5) * (-pkin(3) * t13 - pkin(2)) - m(6) * (-t13 * (pkin(4) * t17 + pkin(3)) - pkin(2)) + (-mrSges(4,2) + m(5) * pkin(6) + mrSges(5,3) - m(6) * (-qJ(5) - pkin(6)) + mrSges(6,3)) * t12;
t16 = sin(qJ(1));
t31 = pkin(1) * t16;
t18 = cos(qJ(1));
t30 = pkin(1) * t18;
t15 = sin(qJ(4));
t11 = qJ(1) + pkin(7);
t9 = sin(t11);
t28 = t15 * t9;
t10 = cos(t11);
t27 = t10 * t15;
t26 = t13 * t15;
t25 = t13 * t17;
t1 = t10 * t17 + t26 * t9;
t3 = t10 * t26 - t9 * t17;
t4 = -t10 * t25 - t28;
t2 = t25 * t9 - t27;
t5 = [(-t27 * t32 + m(3) * t31 + t16 * mrSges(2,1) + mrSges(2,2) * t18 + t24 * (t10 * qJ(3) - t31) + t37 * t2 + t36 * t10 - t35 * t1 + t33 * t9) * g(3) + (t28 * t32 + m(3) * t30 + mrSges(2,1) * t18 - t16 * mrSges(2,2) - t36 * t9 - t37 * t4 - t35 * t3 + t24 * (-t9 * qJ(3) - t30) + t33 * t10) * g(2), (-m(3) + t24) * g(1), (g(2) * t10 + g(3) * t9) * t24, (t15 * t34 + t17 * t35) * g(1) * t12 + (t3 * t34 - t35 * t4) * g(3) + (-t1 * t34 - t2 * t35) * g(2), (g(1) * t13 + (g(2) * t9 - g(3) * t10) * t12) * m(6)];
taug = t5(:);
