% Calculate potential energy for
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:41
% EndTime: 2019-12-31 17:50:41
% DurationCPUTime: 0.29s
% Computational Cost: add. (129->51), mult. (101->34), div. (0->0), fcn. (65->6), ass. (0->19)
t30 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t29 = mrSges(5,2) + mrSges(6,2);
t28 = -m(1) - m(2);
t27 = -m(5) - m(6);
t26 = m(4) - t27;
t24 = -m(5) * pkin(6) + m(6) * (-qJ(5) - pkin(6)) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t11 = sin(qJ(4));
t13 = cos(qJ(4));
t23 = t11 * t30 - t13 * t29 + mrSges(3,2) - mrSges(4,3);
t21 = pkin(5) + r_base(3);
t12 = sin(qJ(1));
t20 = t12 * pkin(1) + r_base(2);
t14 = cos(qJ(1));
t19 = t14 * pkin(1) + r_base(1);
t6 = qJ(2) + t21;
t9 = qJ(1) + pkin(7);
t5 = cos(t9);
t4 = sin(t9);
t1 = (-m(1) * r_base(3) - m(2) * t21 - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t6 + t27 * (pkin(3) + t6) + t30 * t13 + t29 * t11) * g(3) + (-m(3) * t20 - t12 * mrSges(2,1) - mrSges(2,2) * t14 - mrSges(1,2) + t28 * r_base(2) - t26 * (t4 * pkin(2) + t20) + (qJ(3) * t26 - t23) * t5 + t24 * t4) * g(2) + (-m(3) * t19 - mrSges(2,1) * t14 + t12 * mrSges(2,2) - mrSges(1,1) + t28 * r_base(1) - t26 * (t5 * pkin(2) + t4 * qJ(3) + t19) + t24 * t5 + t23 * t4) * g(1);
U = t1;
