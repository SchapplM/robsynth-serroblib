% Calculate Gravitation load on the joints for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:22
% DurationCPUTime: 0.22s
% Computational Cost: add. (116->40), mult. (186->46), div. (0->0), fcn. (182->8), ass. (0->20)
t36 = m(5) + m(6);
t35 = mrSges(2,1) + mrSges(3,1);
t34 = mrSges(2,2) - mrSges(3,3);
t33 = m(4) + t36;
t15 = cos(pkin(8));
t13 = pkin(8) + qJ(5);
t7 = sin(t13);
t8 = cos(t13);
t21 = -mrSges(6,1) * t8 + mrSges(6,2) * t7;
t32 = mrSges(4,1) + m(5) * pkin(3) + t15 * mrSges(5,1) - sin(pkin(8)) * mrSges(5,2) + m(6) * (pkin(4) * t15 + pkin(3)) - t21;
t30 = -m(5) * qJ(4) + m(6) * (-pkin(6) - qJ(4)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t27 = sin(qJ(1));
t28 = cos(qJ(1));
t29 = t28 * pkin(1) + t27 * qJ(2);
t26 = cos(pkin(7));
t25 = sin(pkin(7));
t22 = -t27 * pkin(1) + t28 * qJ(2);
t2 = t28 * t25 - t27 * t26;
t1 = -t27 * t25 - t28 * t26;
t3 = [(-m(3) * t29 - t35 * t28 + t34 * t27 - t33 * (t28 * pkin(2) + t29) + t30 * t2 + t32 * t1) * g(2) + (-m(3) * t22 + t34 * t28 + t35 * t27 - t32 * t2 - t33 * (-t27 * pkin(2) + t22) + t30 * t1) * g(1), (-t27 * g(1) + t28 * g(2)) * (m(3) + t33), t33 * g(3), t36 * (-g(1) * t2 + g(2) * t1), -g(3) * t21 + (-g(1) * t1 - g(2) * t2) * (mrSges(6,1) * t7 + mrSges(6,2) * t8)];
taug = t3(:);
