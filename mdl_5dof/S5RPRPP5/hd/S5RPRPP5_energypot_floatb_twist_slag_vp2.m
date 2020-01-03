% Calculate potential energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:04
% EndTime: 2019-12-31 18:16:04
% DurationCPUTime: 0.35s
% Computational Cost: add. (101->62), mult. (129->46), div. (0->0), fcn. (93->4), ass. (0->19)
t33 = -mrSges(4,1) - mrSges(5,1);
t32 = -mrSges(4,2) + mrSges(6,2);
t31 = -m(1) - m(2);
t30 = -m(5) - m(6);
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t29 = t33 * t11 + t32 * t13 + mrSges(2,2) - mrSges(3,3);
t28 = m(6) * qJ(5) - mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t27 = pkin(3) * t11;
t24 = qJ(4) * t13;
t10 = pkin(5) + r_base(3);
t12 = sin(qJ(1));
t23 = t12 * pkin(1) + r_base(2);
t22 = pkin(2) + t10;
t21 = t12 * pkin(6) + t23;
t14 = cos(qJ(1));
t20 = t14 * pkin(1) + t12 * qJ(2) + r_base(1);
t18 = t14 * pkin(6) + t20;
t1 = (-m(1) * r_base(3) - m(4) * t22 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t30 * (t13 * pkin(3) + t11 * qJ(4) + t22) + (-m(2) - m(3)) * t10 + (-m(6) * pkin(4) - mrSges(6,1) + t33) * t13 + (-mrSges(5,3) - t32) * t11) * g(3) + (-m(3) * t23 - m(4) * t21 - mrSges(1,2) + t31 * r_base(2) + t30 * (t14 * t24 + t21) + (m(5) * t27 - t13 * mrSges(5,3) - (m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1)) * t11 + (m(3) + m(4) - t30) * qJ(2) - t29) * t14 + t28 * t12) * g(2) + (-m(3) * t20 - m(4) * t18 - mrSges(1,1) + t31 * r_base(1) + t30 * (t12 * t27 + t18) + t28 * t14 + (-(-m(5) * qJ(4) - mrSges(5,3)) * t13 - m(6) * (pkin(4) * t11 - t24) - t11 * mrSges(6,1) + t29) * t12) * g(1);
U = t1;
