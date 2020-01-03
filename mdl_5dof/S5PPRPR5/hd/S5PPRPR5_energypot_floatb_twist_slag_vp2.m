% Calculate potential energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:22
% EndTime: 2019-12-31 17:33:22
% DurationCPUTime: 0.30s
% Computational Cost: add. (113->53), mult. (144->41), div. (0->0), fcn. (132->6), ass. (0->20)
t34 = m(5) + m(6);
t33 = -m(1) - m(2);
t31 = -mrSges(2,1) - mrSges(3,1);
t30 = mrSges(2,2) - mrSges(3,3);
t29 = m(6) * pkin(6) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t16 = sin(qJ(5));
t17 = cos(qJ(5));
t28 = t16 * mrSges(6,1) + t17 * mrSges(6,2) + t34 * qJ(4) - mrSges(4,2) + mrSges(5,3);
t27 = cos(qJ(3));
t26 = sin(qJ(3));
t25 = sin(pkin(7));
t14 = qJ(1) + r_base(3);
t15 = cos(pkin(7));
t24 = t15 * pkin(1) + t25 * qJ(2) + r_base(1);
t23 = t15 * pkin(2) + t24;
t21 = t25 * pkin(1) - qJ(2) * t15 + r_base(2);
t20 = t25 * pkin(2) + t21;
t4 = t15 * t26 - t25 * t27;
t3 = -t15 * t27 - t25 * t26;
t1 = (-m(1) * r_base(3) + m(6) * pkin(4) + t17 * mrSges(6,1) - t16 * mrSges(6,2) + mrSges(5,1) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(4) - t34) * (-pkin(5) + t14) + (-m(2) - m(3)) * t14) * g(3) + (-m(3) * t21 - m(4) * t20 - mrSges(1,2) + t33 * r_base(2) + t31 * t25 - t34 * (-t4 * pkin(3) + t20) - t30 * t15 + t29 * t4 + t28 * t3) * g(2) + (-m(3) * t24 - m(4) * t23 - mrSges(1,1) + t33 * r_base(1) + t30 * t25 - t34 * (-t3 * pkin(3) + t23) + t31 * t15 - t28 * t4 + t29 * t3) * g(1);
U = t1;
