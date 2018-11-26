% Calculate Gravitation load on the joints for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:30
% EndTime: 2018-11-23 15:38:31
% DurationCPUTime: 0.35s
% Computational Cost: add. (197->65), mult. (373->78), div. (0->0), fcn. (399->8), ass. (0->34)
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t39 = sin(qJ(1));
t40 = cos(qJ(1));
t5 = -t33 * t39 - t34 * t40;
t6 = t33 * t40 - t34 * t39;
t47 = -g(1) * t6 + g(2) * t5;
t51 = -m(6) - m(7);
t32 = m(5) - t51;
t50 = m(4) + t32;
t49 = mrSges(2,1) + mrSges(3,1);
t48 = mrSges(2,2) - mrSges(3,3);
t46 = mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t17 = cos(qJ(5));
t15 = sin(qJ(5));
t25 = t15 * mrSges(6,1) + t17 * mrSges(6,2);
t43 = t17 * mrSges(7,3) + mrSges(4,2) - mrSges(5,3) - t25 + (-m(5) - m(6)) * qJ(4);
t14 = sin(qJ(6));
t38 = t14 * t15;
t16 = cos(qJ(6));
t37 = t15 * t16;
t35 = t40 * pkin(1) + t39 * qJ(2);
t31 = t40 * pkin(2) + t35;
t30 = m(7) * pkin(8) + mrSges(7,3);
t29 = -t17 * pkin(8) + qJ(4);
t28 = -t5 * pkin(3) + t31;
t27 = -pkin(1) * t39 + t40 * qJ(2);
t24 = mrSges(7,1) * t14 + mrSges(7,2) * t16;
t22 = m(7) * pkin(5) + t16 * mrSges(7,1) - t14 * mrSges(7,2);
t20 = -pkin(2) * t39 + t27;
t19 = t22 * t15;
t2 = -t5 * t14 + t37 * t6;
t1 = -t5 * t16 - t38 * t6;
t3 = [(-m(3) * t35 - m(4) * t31 - m(5) * t28 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t49 * t40 + t48 * t39 + t51 * (-t5 * pkin(7) + t28) + t46 * t5 + (-m(7) * (t15 * pkin(5) + t29) + t43) * t6) * g(2) + (-m(3) * t27 - m(4) * t20 + t48 * t40 + t49 * t39 - t32 * (t6 * pkin(3) + t20) + (t51 * pkin(7) - t24 - t46) * t6 + (-m(7) * t29 - t19 + t43) * t5) * g(1) (-t39 * g(1) + t40 * g(2)) * (m(3) + t50) t50 * g(3), t47 * t32 (t17 * t30 - t19 - t25) * g(3) + ((mrSges(6,1) + t22) * t17 + (-mrSges(6,2) + t30) * t15) * t47, -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * ((-t6 * t16 + t38 * t5) * mrSges(7,1) + (t6 * t14 + t37 * t5) * mrSges(7,2)) - g(3) * t24 * t17];
taug  = t3(:);
