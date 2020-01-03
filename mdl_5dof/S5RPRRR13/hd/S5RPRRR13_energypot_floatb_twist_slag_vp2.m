% Calculate potential energy for
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:27
% EndTime: 2019-12-31 19:14:28
% DurationCPUTime: 0.38s
% Computational Cost: add. (124->54), mult. (155->40), div. (0->0), fcn. (127->8), ass. (0->19)
t11 = sin(qJ(4));
t14 = cos(qJ(4));
t10 = qJ(4) + qJ(5);
t2 = sin(t10);
t3 = cos(t10);
t36 = mrSges(4,1) + m(5) * pkin(3) + t14 * mrSges(5,1) - t11 * mrSges(5,2) + m(6) * (pkin(4) * t14 + pkin(3)) + t3 * mrSges(6,1) - t2 * mrSges(6,2);
t35 = mrSges(4,2) - m(5) * pkin(7) + m(6) * (-pkin(8) - pkin(7)) - mrSges(5,3) - mrSges(6,3);
t12 = sin(qJ(3));
t15 = cos(qJ(3));
t34 = t36 * t12 + t35 * t15 - mrSges(2,2) + mrSges(3,3);
t33 = -m(1) - m(2);
t31 = -m(4) - m(5) - m(6);
t29 = -t2 * mrSges(6,1) - t14 * mrSges(5,2) - t3 * mrSges(6,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t11;
t9 = pkin(5) + r_base(3);
t13 = sin(qJ(1));
t27 = t13 * pkin(1) + r_base(2);
t16 = cos(qJ(1));
t24 = t16 * pkin(1) + t13 * qJ(2) + r_base(1);
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t9 + t31 * (pkin(2) + t9) - t36 * t15 + t35 * t12) * g(3) + (-m(3) * t27 - mrSges(1,2) + t33 * r_base(2) + t31 * (t13 * pkin(6) + t27) + ((m(3) - t31) * qJ(2) + t34) * t16 + t29 * t13) * g(2) + (-m(3) * t24 - mrSges(1,1) + t33 * r_base(1) + t31 * (t16 * pkin(6) + t24) + t29 * t16 - t34 * t13) * g(1);
U = t1;
