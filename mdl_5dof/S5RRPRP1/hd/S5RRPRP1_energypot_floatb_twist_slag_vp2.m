% Calculate potential energy for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:50
% EndTime: 2020-01-03 11:58:51
% DurationCPUTime: 0.31s
% Computational Cost: add. (148->52), mult. (97->37), div. (0->0), fcn. (61->8), ass. (0->22)
t30 = m(5) + m(6);
t29 = -mrSges(5,2) - mrSges(6,2);
t28 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t27 = -m(1) - m(2);
t26 = -m(4) - t30;
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t25 = t30 * pkin(3) + t29 * t12 + t28 * t14 + mrSges(4,1);
t24 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t10 = qJ(1) + qJ(2);
t23 = pkin(5) + r_base(1);
t13 = sin(qJ(1));
t22 = t13 * pkin(1) + r_base(2);
t21 = pkin(6) + t23;
t15 = cos(qJ(1));
t19 = -pkin(1) * t15 + r_base(3);
t8 = cos(t10);
t7 = sin(t10);
t6 = pkin(8) + t10;
t2 = cos(t6);
t1 = sin(t6);
t3 = (-m(3) * t19 + mrSges(2,1) * t15 + t8 * mrSges(3,1) - t13 * mrSges(2,2) - t7 * mrSges(3,2) - mrSges(1,3) + t27 * r_base(3) + t25 * t2 + t26 * (-pkin(2) * t8 + t19) + t24 * t1) * g(3) + (-m(3) * t22 - t13 * mrSges(2,1) - t7 * mrSges(3,1) - mrSges(2,2) * t15 - t8 * mrSges(3,2) - mrSges(1,2) + t27 * r_base(2) + t26 * (pkin(2) * t7 + t22) + t24 * t2 - t25 * t1) * g(2) + (-m(1) * r_base(1) - m(2) * t23 - m(3) * t21 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t29 * t14 + t26 * (qJ(3) + t21) - t28 * t12) * g(1);
U = t3;
