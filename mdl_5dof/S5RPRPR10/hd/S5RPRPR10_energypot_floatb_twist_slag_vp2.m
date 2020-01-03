% Calculate potential energy for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:52
% EndTime: 2019-12-31 18:25:52
% DurationCPUTime: 0.30s
% Computational Cost: add. (138->63), mult. (138->52), div. (0->0), fcn. (120->8), ass. (0->29)
t40 = -m(1) - m(2);
t39 = -m(3) - m(4);
t38 = -m(5) - m(6);
t17 = sin(qJ(5));
t20 = cos(qJ(5));
t37 = m(6) * pkin(4) + mrSges(6,1) * t20 - mrSges(6,2) * t17 + mrSges(5,1);
t36 = mrSges(2,2) - mrSges(3,3);
t35 = -m(4) * pkin(2) - mrSges(2,1) - mrSges(3,1);
t34 = m(6) * pkin(7) - mrSges(5,2) + mrSges(6,3);
t18 = sin(qJ(3));
t19 = sin(qJ(1));
t33 = t19 * t18;
t22 = cos(qJ(1));
t32 = t22 * t18;
t16 = pkin(5) + r_base(3);
t31 = qJ(3) + pkin(8);
t30 = t19 * pkin(1) + r_base(2);
t29 = -pkin(6) + t16;
t28 = t22 * pkin(1) + t19 * qJ(2) + r_base(1);
t27 = cos(t31);
t26 = sin(t31);
t21 = cos(qJ(3));
t11 = pkin(3) * t21 + pkin(2);
t5 = t19 * t11;
t4 = t19 * t21 - t32;
t3 = -t21 * t22 - t33;
t2 = -t19 * t27 + t22 * t26;
t1 = -t19 * t26 - t22 * t27;
t6 = (-m(1) * r_base(3) - m(4) * t29 + t17 * mrSges(6,1) + t20 * mrSges(6,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + mrSges(5,3) + (-m(2) - m(3)) * t16 + t38 * (-qJ(4) + t29)) * g(3) + (-mrSges(1,2) - t4 * mrSges(4,1) - t3 * mrSges(4,2) - m(5) * (t5 + t30) - m(6) * (-pkin(3) * t32 + t5) + t40 * r_base(2) + (-m(6) + t39) * (-qJ(2) * t22 + t30) + (-m(5) * (-pkin(3) * t18 - qJ(2)) - t36) * t22 + t37 * t2 + t35 * t19 + t34 * t1) * g(2) + (t3 * mrSges(4,1) - t4 * mrSges(4,2) - mrSges(1,1) + t40 * r_base(1) + t39 * t28 + t38 * (pkin(3) * t33 + t22 * t11 + t28) + t35 * t22 - t34 * t2 + t36 * t19 + t37 * t1) * g(1);
U = t6;
