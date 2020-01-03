% Calculate potential energy for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:39
% EndTime: 2019-12-31 18:17:40
% DurationCPUTime: 0.23s
% Computational Cost: add. (143->55), mult. (88->41), div. (0->0), fcn. (52->8), ass. (0->21)
t31 = -m(1) - m(2);
t30 = m(5) + m(6);
t15 = sin(qJ(5));
t17 = cos(qJ(5));
t29 = -t15 * mrSges(6,1) - t17 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t28 = -m(6) * pkin(7) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t27 = pkin(5) + r_base(3);
t14 = qJ(1) + pkin(8);
t16 = sin(qJ(1));
t26 = t16 * pkin(1) + r_base(2);
t18 = cos(qJ(1));
t25 = t18 * pkin(1) + r_base(1);
t9 = sin(t14);
t24 = pkin(2) * t9 + t26;
t10 = cos(t14);
t23 = pkin(2) * t10 + t25;
t22 = qJ(2) + t27;
t11 = qJ(3) + t14;
t7 = cos(t11);
t6 = sin(t11);
t1 = (-m(1) * r_base(3) - m(2) * t27 - m(3) * t22 - m(6) * pkin(4) - t17 * mrSges(6,1) + t15 * mrSges(6,2) - mrSges(5,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(4) - t30) * (pkin(6) + t22)) * g(3) + (-m(3) * t26 - m(4) * t24 - t16 * mrSges(2,1) - t9 * mrSges(3,1) - mrSges(2,2) * t18 - t10 * mrSges(3,2) - mrSges(1,2) + t31 * r_base(2) - t30 * (t6 * pkin(3) + t24) + (t30 * qJ(4) - t29) * t7 + t28 * t6) * g(2) + (-m(3) * t25 - m(4) * t23 - mrSges(2,1) * t18 - t10 * mrSges(3,1) + t16 * mrSges(2,2) + t9 * mrSges(3,2) - mrSges(1,1) + t31 * r_base(1) - t30 * (t7 * pkin(3) + t6 * qJ(4) + t23) + t28 * t7 + t29 * t6) * g(1);
U = t1;
