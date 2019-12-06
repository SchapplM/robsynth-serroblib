% Calculate Gravitation load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:49
% EndTime: 2019-12-05 18:31:50
% DurationCPUTime: 0.21s
% Computational Cost: add. (271->46), mult. (182->46), div. (0->0), fcn. (135->10), ass. (0->32)
t22 = qJ(4) + qJ(5);
t20 = cos(t22);
t13 = t20 * mrSges(6,1);
t26 = cos(qJ(4));
t24 = sin(qJ(4));
t42 = mrSges(5,2) * t24;
t49 = -t42 + m(5) * pkin(3) + m(6) * (pkin(4) * t26 + pkin(3)) + t26 * mrSges(5,1) + mrSges(4,1) + t13;
t48 = -m(5) * pkin(7) - mrSges(5,3) - mrSges(6,3) + mrSges(4,2) + m(6) * (-pkin(8) - pkin(7));
t36 = m(6) * pkin(4) + mrSges(5,1);
t45 = -mrSges(5,2) * t26 - t24 * t36;
t23 = qJ(1) + qJ(2);
t17 = pkin(9) + t23;
t14 = sin(t17);
t18 = sin(t22);
t38 = t14 * t18;
t39 = mrSges(6,2) * t20;
t44 = mrSges(6,1) * t38 + t14 * t39;
t15 = cos(t17);
t43 = g(3) * t15;
t40 = mrSges(6,2) * t18;
t37 = m(4) + m(5) + m(6);
t35 = t13 - t40;
t34 = -mrSges(6,1) * t18 - t39;
t33 = pkin(2) * t37 + mrSges(3,1);
t32 = mrSges(2,1) + (m(3) + t37) * pkin(1);
t19 = sin(t23);
t21 = cos(t23);
t30 = -t19 * mrSges(3,2) + t33 * t21 + (-t40 + t49) * t15 - t48 * t14;
t29 = mrSges(3,2) * t21 - mrSges(6,2) * t38 + t49 * t14 + t48 * t15 + t33 * t19;
t27 = cos(qJ(1));
t25 = sin(qJ(1));
t1 = [(mrSges(2,2) * t27 + t25 * t32 + t29) * g(3) + (-mrSges(2,2) * t25 + t27 * t32 + t30) * g(2), g(2) * t30 + g(3) * t29, -t37 * g(1), (t45 * t14 - t44) * g(2) + (-t26 * t36 - t35 + t42) * g(1) + (-t34 - t45) * t43, -g(1) * t35 - g(2) * t44 - t34 * t43];
taug = t1(:);
