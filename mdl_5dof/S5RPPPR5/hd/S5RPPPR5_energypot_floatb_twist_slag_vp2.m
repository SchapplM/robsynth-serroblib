% Calculate potential energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:18
% EndTime: 2019-12-31 17:46:18
% DurationCPUTime: 0.29s
% Computational Cost: add. (124->55), mult. (153->44), div. (0->0), fcn. (141->8), ass. (0->21)
t35 = -m(1) - m(2);
t34 = -mrSges(2,1) - mrSges(3,1);
t33 = mrSges(2,2) - mrSges(3,3);
t32 = -m(4) - m(5) - m(6);
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t15 = pkin(8) + qJ(5);
t7 = sin(t15);
t8 = cos(t15);
t31 = m(5) * pkin(3) + t18 * mrSges(5,1) - t17 * mrSges(5,2) + mrSges(4,1) + m(6) * (pkin(4) * t18 + pkin(3)) + t8 * mrSges(6,1) - t7 * mrSges(6,2);
t30 = m(5) * qJ(4) - m(6) * (-pkin(6) - qJ(4)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t29 = sin(qJ(1));
t28 = cos(pkin(7));
t27 = sin(pkin(7));
t16 = pkin(5) + r_base(3);
t20 = cos(qJ(1));
t26 = t20 * pkin(1) + t29 * qJ(2) + r_base(1);
t24 = t29 * pkin(1) - t20 * qJ(2) + r_base(2);
t2 = t20 * t27 - t29 * t28;
t1 = -t20 * t28 - t29 * t27;
t3 = (-m(1) * r_base(3) + t7 * mrSges(6,1) + mrSges(5,2) * t18 + t8 * mrSges(6,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t17 + t32 * (-qJ(3) + t16) + (-m(2) - m(3)) * t16) * g(3) + (-m(3) * t24 - mrSges(1,2) + t35 * r_base(2) + t34 * t29 - t33 * t20 + t32 * (t29 * pkin(2) + t24) + t31 * t2 + t30 * t1) * g(2) + (-m(3) * t26 - mrSges(1,1) + t35 * r_base(1) + t33 * t29 + t34 * t20 + t32 * (t20 * pkin(2) + t26) - t30 * t2 + t31 * t1) * g(1);
U = t3;
