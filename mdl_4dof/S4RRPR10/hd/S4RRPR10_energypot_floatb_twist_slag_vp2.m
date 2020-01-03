% Calculate potential energy for
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRPR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:05
% EndTime: 2019-12-31 17:11:06
% DurationCPUTime: 0.34s
% Computational Cost: add. (84->59), mult. (122->56), div. (0->0), fcn. (96->6), ass. (0->26)
t36 = -mrSges(3,1) + mrSges(4,2);
t35 = mrSges(3,2) - mrSges(4,3);
t34 = -m(1) - m(2);
t33 = m(4) + m(5);
t12 = sin(qJ(1));
t11 = sin(qJ(2));
t23 = qJ(3) * t11;
t14 = cos(qJ(2));
t26 = t12 * t14;
t32 = pkin(2) * t26 + t12 * t23;
t31 = t35 * t11 + t36 * t14 - mrSges(2,1);
t30 = mrSges(4,1) + mrSges(3,3) - mrSges(2,2);
t10 = sin(qJ(4));
t15 = cos(qJ(1));
t29 = t10 * t15;
t28 = t12 * t10;
t13 = cos(qJ(4));
t27 = t12 * t13;
t25 = t13 * t15;
t24 = t15 * t14;
t9 = pkin(4) + r_base(3);
t22 = t12 * pkin(1) + r_base(2);
t20 = t15 * pkin(1) + t12 * pkin(5) + r_base(1);
t19 = pkin(2) * t24 + t15 * t23 + t20;
t18 = -pkin(5) * t15 + t22;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t9 - t33 * (t11 * pkin(2) + t9) + (t10 * mrSges(5,1) + t13 * mrSges(5,2) + t33 * qJ(3) - t35) * t14 + (-m(5) * pkin(6) - mrSges(5,3) + t36) * t11) * g(3) + (-mrSges(1,2) - m(3) * t18 - m(4) * (t18 + t32) - m(5) * (pkin(6) * t26 + t22 + t32) - (t11 * t28 - t25) * mrSges(5,1) - (t11 * t27 + t29) * mrSges(5,2) - mrSges(5,3) * t26 + t34 * r_base(2) + (-m(5) * (-pkin(3) - pkin(5)) + t30) * t15 + t31 * t12) * g(2) + (-mrSges(1,1) - m(3) * t20 - m(4) * t19 - m(5) * (pkin(6) * t24 + t19) - (t11 * t29 + t27) * mrSges(5,1) - (t11 * t25 - t28) * mrSges(5,2) - mrSges(5,3) * t24 + t34 * r_base(1) + t31 * t15 + (-m(5) * pkin(3) - t30) * t12) * g(1);
U = t1;
