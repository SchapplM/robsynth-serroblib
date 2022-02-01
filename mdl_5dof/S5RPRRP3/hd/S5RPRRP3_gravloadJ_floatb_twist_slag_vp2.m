% Calculate Gravitation load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:35
% EndTime: 2022-01-23 09:29:36
% DurationCPUTime: 0.30s
% Computational Cost: add. (203->47), mult. (182->46), div. (0->0), fcn. (134->8), ass. (0->26)
t16 = sin(qJ(3));
t15 = qJ(3) + qJ(4);
t10 = cos(t15);
t29 = mrSges(5,2) + mrSges(6,2);
t30 = mrSges(5,1) + mrSges(6,1);
t9 = sin(t15);
t26 = -t30 * t10 + t29 * t9;
t43 = -t16 * mrSges(4,2) - t26;
t38 = m(3) + m(4) + m(5) + m(6);
t40 = pkin(1) * t38 + mrSges(2,1);
t14 = qJ(1) + pkin(8);
t7 = sin(t14);
t8 = cos(t14);
t39 = g(1) * t8 + g(2) * t7;
t18 = cos(qJ(3));
t11 = t18 * pkin(3);
t5 = pkin(4) * t10;
t32 = t5 + t11;
t37 = m(4) * pkin(2) + t18 * mrSges(4,1) + mrSges(3,1) + m(6) * (pkin(2) + t32) + m(5) * (t11 + pkin(2)) + t43;
t20 = -pkin(7) - pkin(6);
t36 = -m(4) * pkin(6) + m(5) * t20 + m(6) * (-qJ(5) + t20) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t28 = m(5) * pkin(3) + mrSges(4,1);
t27 = t29 * t10;
t19 = cos(qJ(1));
t17 = sin(qJ(1));
t1 = [(mrSges(2,2) * t17 - t40 * t19 + t36 * t7 - t37 * t8) * g(2) + (mrSges(2,2) * t19 + t40 * t17 + t36 * t8 + t37 * t7) * g(1), -t38 * g(3), (-m(6) * t32 - t28 * t18 - t43) * g(3) + t39 * (-m(6) * (-pkin(3) * t16 - pkin(4) * t9) + mrSges(4,2) * t18 + t28 * t16 + t30 * t9 + t27), (-m(6) * t5 + t26) * g(3) + t39 * (t27 + (m(6) * pkin(4) + t30) * t9), (-g(1) * t7 + g(2) * t8) * m(6)];
taug = t1(:);
