% Calculate Gravitation load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:37
% EndTime: 2019-12-05 17:41:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (185->43), mult. (145->44), div. (0->0), fcn. (105->10), ass. (0->27)
t28 = -m(4) - m(5) - m(6);
t42 = -m(3) + t28;
t45 = -t42 * pkin(1) + mrSges(2,1);
t15 = pkin(9) + qJ(4);
t12 = qJ(5) + t15;
t5 = sin(t12);
t6 = cos(t12);
t26 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t8 = sin(t15);
t44 = -t8 * mrSges(5,2) + t26;
t43 = mrSges(6,1) * t5 + mrSges(6,2) * t6;
t10 = cos(t15);
t27 = m(6) * pkin(4) + mrSges(5,1);
t41 = -mrSges(5,2) * t10 - t27 * t8;
t18 = cos(pkin(9));
t7 = t18 * pkin(3) + pkin(2);
t40 = mrSges(3,1) + m(4) * pkin(2) + t18 * mrSges(4,1) - sin(pkin(9)) * mrSges(4,2) + m(5) * t7 + t10 * mrSges(5,1) + m(6) * (pkin(4) * t10 + t7) + t44;
t19 = -pkin(6) - qJ(3);
t38 = m(4) * qJ(3) - m(5) * t19 - m(6) * (-pkin(7) + t19) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t16 = qJ(1) + pkin(8);
t9 = sin(t16);
t37 = t43 * t9;
t11 = cos(t16);
t32 = g(3) * t11;
t21 = cos(qJ(1));
t20 = sin(qJ(1));
t1 = [(mrSges(2,2) * t21 - t38 * t11 + t45 * t20 + t40 * t9) * g(3) + (-t20 * mrSges(2,2) + t40 * t11 + t45 * t21 + t38 * t9) * g(2), t42 * g(1), (g(2) * t11 + g(3) * t9) * t28, (t41 * t9 - t37) * g(2) + (-t27 * t10 - t44) * g(1) + (t43 - t41) * t32, -g(1) * t26 - g(2) * t37 + t32 * t43];
taug = t1(:);
