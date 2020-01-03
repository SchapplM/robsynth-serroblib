% Calculate Gravitation load on the joints for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:20
% DurationCPUTime: 0.30s
% Computational Cost: add. (153->58), mult. (298->74), div. (0->0), fcn. (319->8), ass. (0->29)
t45 = -m(5) - m(6);
t44 = mrSges(2,1) + mrSges(3,1);
t43 = mrSges(2,2) - mrSges(3,3);
t42 = mrSges(4,2) - mrSges(5,3);
t40 = m(4) - t45;
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t14 = sin(qJ(5));
t16 = cos(qJ(5));
t21 = m(6) * pkin(4) + mrSges(6,1) * t16 - mrSges(6,2) * t14;
t25 = t17 * mrSges(5,1) - t15 * mrSges(5,2);
t27 = m(6) * pkin(7) + mrSges(6,3);
t38 = t15 * t27 + t17 * t21 + t25;
t35 = cos(qJ(1));
t34 = sin(qJ(1));
t33 = t14 * t17;
t32 = t16 * t17;
t31 = t35 * pkin(1) + t34 * qJ(2);
t30 = cos(pkin(8));
t29 = sin(pkin(8));
t28 = t35 * pkin(2) + t31;
t26 = -pkin(1) * t34 + t35 * qJ(2);
t23 = mrSges(6,1) * t14 + mrSges(6,2) * t16;
t20 = -pkin(2) * t34 + t26;
t6 = t29 * t35 - t30 * t34;
t5 = -t29 * t34 - t30 * t35;
t2 = t14 * t6 - t32 * t5;
t1 = t16 * t6 + t33 * t5;
t3 = [(-m(3) * t31 - m(4) * t28 - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t42 * t6 - t44 * t35 + t43 * t34 + t45 * (-t5 * pkin(3) + t6 * pkin(6) + t28) + (mrSges(4,1) + t25 - m(6) * (-pkin(4) * t17 - pkin(7) * t15) + t15 * mrSges(6,3)) * t5) * g(2) + (-m(3) * t26 - m(4) * t20 + t43 * t35 + t44 * t34 + t45 * (t6 * pkin(3) + t20) + (-mrSges(4,1) - t38) * t6 + (pkin(6) * t45 - t23 + t42) * t5) * g(1), (-g(1) * t34 + g(2) * t35) * (m(3) + t40), t40 * g(3), t38 * g(3) + ((mrSges(5,2) - t27) * t17 + (mrSges(5,1) + t21) * t15) * (-g(1) * t5 - g(2) * t6), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((-t16 * t5 + t33 * t6) * mrSges(6,1) + (t14 * t5 + t32 * t6) * mrSges(6,2)) - g(3) * t23 * t15];
taug = t3(:);
