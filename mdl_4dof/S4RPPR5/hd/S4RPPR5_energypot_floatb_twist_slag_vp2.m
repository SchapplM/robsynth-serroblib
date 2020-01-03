% Calculate potential energy for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:41
% EndTime: 2019-12-31 16:39:41
% DurationCPUTime: 0.20s
% Computational Cost: add. (80->44), mult. (102->36), div. (0->0), fcn. (86->6), ass. (0->18)
t29 = -m(1) - m(2);
t28 = -m(4) - m(5);
t27 = -mrSges(2,1) - mrSges(3,1);
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t26 = m(5) * pkin(3) + t14 * mrSges(5,1) - t13 * mrSges(5,2) + mrSges(4,1);
t25 = mrSges(2,2) - mrSges(3,3);
t24 = m(5) * pkin(5) - mrSges(4,2) + mrSges(5,3);
t23 = sin(qJ(1));
t22 = cos(pkin(6));
t21 = sin(pkin(6));
t12 = pkin(4) + r_base(3);
t15 = cos(qJ(1));
t20 = t15 * pkin(1) + t23 * qJ(2) + r_base(1);
t18 = t23 * pkin(1) - qJ(2) * t15 + r_base(2);
t2 = t15 * t21 - t23 * t22;
t1 = -t15 * t22 - t23 * t21;
t3 = (-m(1) * r_base(3) + mrSges(5,1) * t13 + mrSges(5,2) * t14 - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t28 * (-qJ(3) + t12) + (-m(2) - m(3)) * t12) * g(3) + (-m(3) * t18 - mrSges(1,2) + t29 * r_base(2) + t27 * t23 + t26 * t2 + t28 * (t23 * pkin(2) + t18) - t25 * t15 + t24 * t1) * g(2) + (-m(3) * t20 - mrSges(1,1) + t29 * r_base(1) + t25 * t23 + t28 * (t15 * pkin(2) + t20) - t24 * t2 + t27 * t15 + t26 * t1) * g(1);
U = t3;
