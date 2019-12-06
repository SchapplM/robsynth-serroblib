% Calculate Gravitation load on the joints for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:46
% EndTime: 2019-12-05 18:44:47
% DurationCPUTime: 0.36s
% Computational Cost: add. (291->61), mult. (271->57), div. (0->0), fcn. (206->8), ass. (0->32)
t21 = sin(qJ(2));
t20 = qJ(2) + qJ(3);
t14 = sin(t20);
t15 = cos(t20);
t17 = qJ(4) + t20;
t11 = sin(t17);
t12 = cos(t17);
t44 = mrSges(5,2) + mrSges(6,2);
t45 = mrSges(5,1) + mrSges(6,1);
t35 = t44 * t11 - t45 * t12;
t30 = -t15 * mrSges(4,1) + t14 * mrSges(4,2) + t35;
t58 = -t21 * mrSges(3,2) - t30;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t54 = g(1) * t24 + g(2) * t22;
t23 = cos(qJ(2));
t18 = t23 * pkin(2);
t10 = pkin(3) * t15;
t43 = t10 + t18;
t8 = pkin(4) * t12;
t42 = t8 + t43;
t53 = mrSges(2,1) + m(4) * (t18 + pkin(1)) + m(6) * (pkin(1) + t42) + m(5) * (pkin(1) + t43) + m(3) * pkin(1) + t23 * mrSges(3,1) + t58;
t25 = -pkin(7) - pkin(6);
t19 = -pkin(8) + t25;
t52 = mrSges(2,2) + m(6) * (-qJ(5) + t19) - mrSges(6,3) + m(5) * t19 - mrSges(5,3) + m(4) * t25 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t51 = pkin(3) * t14;
t48 = t21 * pkin(2);
t40 = m(4) * pkin(2) + mrSges(3,1);
t36 = t44 * t12;
t3 = -pkin(4) * t11 - t51;
t28 = mrSges(4,2) * t15 + t45 * t11 + t36;
t1 = [(t52 * t22 - t53 * t24) * g(2) + (t53 * t22 + t52 * t24) * g(1), (-m(5) * t43 - m(6) * t42 - t40 * t23 - t58) * g(3) + t54 * (-m(5) * (-t48 - t51) - m(6) * (t3 - t48) + mrSges(4,1) * t14 + mrSges(3,2) * t23 + t40 * t21 + t28), (-m(5) * t10 - m(6) * (t8 + t10) + t30) * g(3) + t54 * (-m(6) * t3 + (m(5) * pkin(3) + mrSges(4,1)) * t14 + t28), (-m(6) * t8 + t35) * g(3) + t54 * (t36 + (m(6) * pkin(4) + t45) * t11), (-g(1) * t22 + g(2) * t24) * m(6)];
taug = t1(:);
