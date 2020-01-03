% Calculate potential energy for
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:02
% EndTime: 2019-12-31 19:38:02
% DurationCPUTime: 0.50s
% Computational Cost: add. (138->67), mult. (191->58), div. (0->0), fcn. (167->8), ass. (0->29)
t48 = mrSges(3,2) - mrSges(4,3);
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t47 = pkin(2) * t19 + qJ(3) * t17;
t15 = cos(pkin(8));
t46 = -m(6) * (pkin(4) * t15 + pkin(3)) - mrSges(3,1) - mrSges(4,1) - m(5) * pkin(3);
t42 = -m(5) - m(6);
t45 = -m(4) + t42;
t43 = -m(1) - m(2);
t18 = sin(qJ(1));
t41 = t47 * t18;
t14 = sin(pkin(8));
t37 = t14 * t17;
t25 = t15 * t19 + t37;
t26 = -t14 * t19 + t15 * t17;
t12 = pkin(8) + qJ(5);
t6 = sin(t12);
t7 = cos(t12);
t30 = t17 * t6 + t19 * t7;
t31 = t17 * t7 - t19 * t6;
t40 = -m(6) * pkin(4) * t37 - mrSges(5,1) * t25 - t30 * mrSges(6,1) - mrSges(5,2) * t26 - t31 * mrSges(6,2) + t48 * t17 + t46 * t19 - mrSges(2,1);
t39 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3) - mrSges(6,3);
t13 = pkin(5) + r_base(3);
t35 = t18 * pkin(1) + r_base(2);
t34 = t17 * pkin(2) + t13;
t20 = cos(qJ(1));
t24 = -pkin(6) * t20 + t35;
t16 = -pkin(7) - qJ(4);
t1 = (-m(1) * r_base(3) - m(6) * t34 - t26 * mrSges(5,1) - t31 * mrSges(6,1) + t25 * mrSges(5,2) + t30 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) + (-m(4) - m(5)) * (-qJ(3) * t19 + t34) + (-m(6) * (-pkin(4) * t14 - qJ(3)) - t48) * t19 + t46 * t17 + (-m(2) - m(3)) * t13) * g(3) + (-mrSges(1,2) - m(3) * t24 - m(4) * (t24 + t41) + t43 * r_base(2) + t42 * (t35 + t41) + (-m(5) * (-pkin(6) + qJ(4)) - m(6) * (-pkin(6) - t16) + t39) * t20 + t40 * t18) * g(2) + (-mrSges(1,1) + t43 * r_base(1) + (m(5) * qJ(4) - m(6) * t16 - t39) * t18 + (-m(3) + t45) * (t20 * pkin(1) + t18 * pkin(6) + r_base(1)) + (t45 * t47 + t40) * t20) * g(1);
U = t1;
