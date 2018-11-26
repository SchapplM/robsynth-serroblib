% Calculate Gravitation load on the joints for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2018-11-23 15:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:55
% EndTime: 2018-11-23 15:38:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (178->66), mult. (338->81), div. (0->0), fcn. (345->8), ass. (0->32)
t47 = -m(6) - m(7);
t46 = -mrSges(5,2) + mrSges(6,3);
t16 = sin(qJ(5));
t19 = cos(qJ(5));
t15 = sin(qJ(6));
t18 = cos(qJ(6));
t22 = m(7) * pkin(5) + t18 * mrSges(7,1) - t15 * mrSges(7,2);
t28 = t19 * mrSges(6,1) - t16 * mrSges(6,2);
t31 = m(7) * pkin(8) + mrSges(7,3);
t45 = -t31 * t16 - t22 * t19 - t28;
t44 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t43 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t39 = t20 * pkin(1) + t17 * qJ(2);
t38 = t15 * t19;
t37 = t18 * t19;
t36 = -pkin(1) - qJ(3);
t35 = sin(pkin(9));
t34 = -m(5) + t47;
t33 = t20 * qJ(3) + t39;
t32 = -m(4) + t34;
t30 = t17 * pkin(3) + t33;
t26 = -mrSges(7,1) * t15 - mrSges(7,2) * t18;
t10 = t20 * qJ(2);
t25 = t20 * pkin(3) + t36 * t17 + t10;
t14 = cos(pkin(9));
t6 = t17 * t14 + t20 * t35;
t5 = t20 * t14 - t17 * t35;
t2 = -t5 * t15 + t6 * t37;
t1 = -t5 * t18 - t6 * t38;
t3 = [(-m(3) * t39 - m(4) * t33 - m(5) * t30 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t46 * t5 + t47 * (t6 * pkin(4) - t5 * pkin(7) + t30) + (-mrSges(5,1) - t28 - m(7) * (t19 * pkin(5) + t16 * pkin(8)) - t16 * mrSges(7,3)) * t6 + t44 * t20 + t43 * t17) * g(2) + (-m(5) * t25 + t47 * (t5 * pkin(4) + t25) + (-m(3) - m(4)) * t10 + (t47 * pkin(7) + t26 - t46) * t6 + (-mrSges(5,1) + t45) * t5 + t43 * t20 + (m(3) * pkin(1) - m(4) * t36 - t44) * t17) * g(1) (-g(1) * t17 + g(2) * t20) * (m(3) - t32) (g(1) * t20 + g(2) * t17) * t32, t34 * g(3), t45 * g(3) + ((mrSges(6,2) - t31) * t19 + (mrSges(6,1) + t22) * t16) * (g(1) * t6 - g(2) * t5) -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * ((-t6 * t18 + t5 * t38) * mrSges(7,1) + (t6 * t15 + t37 * t5) * mrSges(7,2)) - g(3) * t26 * t16];
taug  = t3(:);
