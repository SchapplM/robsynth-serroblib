% Calculate Gravitation load on the joints for
% S5RRPRR3
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:18
% EndTime: 2019-12-05 18:30:18
% DurationCPUTime: 0.20s
% Computational Cost: add. (284->35), mult. (156->36), div. (0->0), fcn. (112->10), ass. (0->26)
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t37 = mrSges(6,1) * t19 - mrSges(6,2) * t17;
t36 = m(6) * pkin(4) + mrSges(5,1) + t37;
t35 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t33 = m(5) + m(6);
t16 = qJ(1) + qJ(2);
t30 = m(4) + t33;
t13 = pkin(9) + t16;
t29 = t33 * pkin(3) + mrSges(4,1);
t27 = t30 * pkin(2) + mrSges(3,1);
t25 = mrSges(2,1) + (m(3) + t30) * pkin(1);
t12 = qJ(4) + t13;
t7 = sin(t12);
t8 = cos(t12);
t24 = t35 * t8 + t36 * t7;
t23 = -t35 * t7 + t36 * t8;
t10 = sin(t13);
t11 = cos(t13);
t14 = sin(t16);
t15 = cos(t16);
t22 = -t14 * mrSges(3,2) - t10 * mrSges(4,2) + t29 * t11 + t27 * t15 + t23;
t21 = mrSges(3,2) * t15 + t11 * mrSges(4,2) + t29 * t10 + t27 * t14 + t24;
t20 = cos(qJ(1));
t18 = sin(qJ(1));
t1 = [(mrSges(2,2) * t20 + t25 * t18 + t21) * g(3) + (-t18 * mrSges(2,2) + t25 * t20 + t22) * g(2), t22 * g(2) + t21 * g(3), -t30 * g(1), t23 * g(2) + t24 * g(3), -g(1) * t37 + (-g(2) * t7 + g(3) * t8) * (mrSges(6,1) * t17 + mrSges(6,2) * t19)];
taug = t1(:);
