% Calculate Gravitation load on the joints for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:17
% DurationCPUTime: 0.22s
% Computational Cost: add. (158->42), mult. (171->43), div. (0->0), fcn. (130->8), ass. (0->24)
t15 = sin(qJ(2));
t14 = qJ(2) + qJ(3);
t10 = cos(t14);
t11 = qJ(4) + t14;
t6 = sin(t11);
t7 = cos(t11);
t30 = t7 * mrSges(5,1) - t6 * mrSges(5,2);
t9 = sin(t14);
t26 = -t10 * mrSges(4,1) + t9 * mrSges(4,2) - t30;
t43 = -t15 * mrSges(3,2) - t26;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t41 = g(1) * t18 + g(2) * t16;
t17 = cos(qJ(2));
t12 = t17 * pkin(2);
t5 = pkin(3) * t10;
t34 = t5 + t12;
t40 = m(3) * pkin(1) + t17 * mrSges(3,1) + mrSges(2,1) + m(4) * (t12 + pkin(1)) + m(5) * (pkin(1) + t34) + t43;
t19 = -pkin(6) - pkin(5);
t39 = mrSges(2,2) + m(5) * (-pkin(7) + t19) - mrSges(5,3) + m(4) * t19 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t31 = m(4) * pkin(2) + mrSges(3,1);
t27 = mrSges(5,1) * t6 + mrSges(5,2) * t7;
t22 = mrSges(4,2) * t10 + t27;
t1 = [(t39 * t16 - t40 * t18) * g(2) + (t40 * t16 + t39 * t18) * g(1), (-m(5) * t34 - t31 * t17 - t43) * g(3) + t41 * (-m(5) * (-pkin(2) * t15 - pkin(3) * t9) + mrSges(4,1) * t9 + mrSges(3,2) * t17 + t31 * t15 + t22), (-m(5) * t5 + t26) * g(3) + t41 * ((m(5) * pkin(3) + mrSges(4,1)) * t9 + t22), -g(3) * t30 + t41 * t27];
taug = t1(:);
