% Calculate potential energy for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:08
% EndTime: 2020-01-03 11:49:09
% DurationCPUTime: 0.57s
% Computational Cost: add. (160->74), mult. (192->67), div. (0->0), fcn. (168->8), ass. (0->29)
t17 = sin(qJ(3));
t19 = cos(qJ(3));
t43 = -m(4) * pkin(2) - t19 * mrSges(4,1) + t17 * mrSges(4,2) - mrSges(3,1);
t42 = mrSges(3,2) - mrSges(5,3) - mrSges(6,3);
t41 = -m(2) - m(4);
t40 = -mrSges(5,1) - mrSges(6,1);
t39 = mrSges(5,2) + mrSges(6,2);
t38 = -m(3) - m(5) - m(6);
t37 = m(4) * pkin(6) + mrSges(4,3);
t15 = sin(pkin(8));
t16 = cos(pkin(8));
t36 = t42 * t15 + t43 * t16 - mrSges(2,1);
t34 = t17 * pkin(3);
t14 = qJ(3) + qJ(4);
t8 = sin(t14);
t35 = m(5) * t34 + m(6) * (pkin(4) * t8 + t34) - mrSges(2,2) + mrSges(3,3) + m(4) * qJ(2) + t17 * mrSges(4,1) + t19 * mrSges(4,2);
t21 = -pkin(7) - pkin(6);
t7 = t19 * pkin(3) + pkin(2);
t18 = sin(qJ(1));
t31 = t16 * t18;
t20 = cos(qJ(1));
t30 = t16 * t20;
t29 = t18 * pkin(1) + r_base(2);
t12 = -qJ(5) + t21;
t9 = cos(t14);
t5 = pkin(4) * t9 + t7;
t27 = -t12 * t15 + t16 * t5;
t26 = -t15 * t21 + t16 * t7;
t1 = (-mrSges(1,3) + t40 * (-t18 * t8 - t9 * t30) - t39 * (-t18 * t9 + t8 * t30) + (-m(1) + t41) * r_base(3) + t38 * (-t18 * qJ(2) + r_base(3)) + t35 * t18 + (m(3) * pkin(1) - m(4) * (-pkin(6) * t15 - pkin(1)) + t15 * mrSges(4,3) - m(5) * (-pkin(1) - t26) - m(6) * (-pkin(1) - t27) - t36) * t20) * g(3) + (-m(4) * t29 - mrSges(1,2) + (-m(1) - m(2)) * r_base(2) + t40 * (-t20 * t8 + t9 * t31) - t39 * (-t20 * t9 - t8 * t31) + t38 * (-t20 * qJ(2) + t29) + t35 * t20 + (-m(5) * t26 - m(6) * t27 - t37 * t15 + t36) * t18) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - mrSges(2,3) + (t38 + t41) * (pkin(5) + r_base(1)) + (-m(5) * t21 - m(6) * t12 + t37 - t42) * t16 + (-m(5) * t7 - m(6) * t5 + t39 * t8 + t40 * t9 + t43) * t15) * g(1);
U = t1;
