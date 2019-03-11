% Calculate potential energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRPR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:13
% EndTime: 2019-03-09 01:55:13
% DurationCPUTime: 0.53s
% Computational Cost: add. (176->80), mult. (181->59), div. (0->0), fcn. (143->8), ass. (0->31)
t44 = -mrSges(5,1) + mrSges(6,2);
t43 = -m(1) - m(2);
t42 = m(3) + m(4);
t41 = m(6) + m(7);
t17 = sin(qJ(6));
t19 = cos(qJ(6));
t40 = -t17 * mrSges(7,1) - t19 * mrSges(7,2) - mrSges(6,3);
t39 = m(7) * pkin(8) + mrSges(7,3);
t14 = sin(pkin(9));
t15 = cos(pkin(9));
t12 = pkin(9) + qJ(4);
t6 = sin(t12);
t7 = cos(t12);
t38 = -mrSges(4,1) * t14 - mrSges(4,2) * t15 - t7 * mrSges(5,2) + t44 * t6 + mrSges(2,2) - mrSges(3,3);
t16 = -pkin(7) - qJ(3);
t37 = -m(4) * qJ(3) - mrSges(2,1) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - m(7) * (pkin(5) - t16) - t19 * mrSges(7,1) + t17 * mrSges(7,2);
t36 = pkin(4) * t6;
t35 = pkin(3) * t14;
t13 = pkin(6) + r_base(3);
t18 = sin(qJ(1));
t33 = t18 * pkin(1) + r_base(2);
t32 = pkin(2) + t13;
t20 = cos(qJ(1));
t31 = t20 * pkin(1) + t18 * qJ(2) + r_base(1);
t29 = t15 * pkin(3) + t32;
t28 = t18 * t35 + t31;
t24 = -t18 * t16 + t33;
t22 = -t16 * t20 + t28;
t2 = t18 * t36;
t1 = t20 * t7 * qJ(5);
t3 = (-m(1) * r_base(3) - m(4) * t32 - m(5) * t29 - t15 * mrSges(4,1) + t14 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t41 * (t7 * pkin(4) + t6 * qJ(5) + t29) + (-m(2) - m(3)) * t13 + (-t39 + t44) * t7 + (mrSges(5,2) + t40) * t6) * g(3) + (-mrSges(1,2) - m(5) * t24 - m(6) * (t1 + t24) - m(7) * t1 + t43 * r_base(2) + (-m(7) - t42) * t33 + (m(6) * t36 - (m(7) * (-pkin(4) - pkin(8)) - mrSges(7,3)) * t6 + t40 * t7 + (-m(5) - t41) * (-qJ(2) - t35) + t42 * qJ(2) - t38) * t20 + t37 * t18) * g(2) + (-mrSges(1,1) - m(5) * t22 - m(6) * (t2 + t22) - m(7) * (t2 + t28) + t43 * r_base(1) - t42 * t31 + t37 * t20 + (-t39 * t6 + (t41 * qJ(5) - t40) * t7 + t38) * t18) * g(1);
U  = t3;
