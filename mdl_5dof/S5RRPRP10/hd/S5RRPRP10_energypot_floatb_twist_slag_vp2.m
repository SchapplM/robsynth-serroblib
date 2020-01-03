% Calculate potential energy for
% S5RRPRP10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:19
% EndTime: 2019-12-31 20:09:19
% DurationCPUTime: 0.47s
% Computational Cost: add. (124->67), mult. (184->59), div. (0->0), fcn. (156->6), ass. (0->29)
t16 = sin(qJ(4));
t49 = -m(6) * pkin(4) * t16 + mrSges(3,2) - mrSges(4,3);
t48 = m(6) * (-qJ(5) - pkin(7)) - mrSges(3,1) + mrSges(4,2) - mrSges(6,3);
t47 = -m(1) - m(2);
t46 = -m(4) - m(6);
t18 = sin(qJ(1));
t17 = sin(qJ(2));
t31 = qJ(3) * t17;
t20 = cos(qJ(2));
t35 = t18 * t20;
t45 = pkin(2) * t35 + t18 * t31;
t44 = mrSges(5,1) + mrSges(6,1);
t43 = mrSges(5,2) + mrSges(6,2);
t42 = m(5) - t46;
t41 = -m(5) * pkin(7) - mrSges(5,3);
t40 = t49 * t17 + t48 * t20 - mrSges(2,1);
t19 = cos(qJ(4));
t39 = m(6) * (pkin(4) * t19 + pkin(3)) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t37 = t18 * t16;
t36 = t18 * t19;
t21 = cos(qJ(1));
t34 = t21 * t16;
t33 = t21 * t19;
t32 = t21 * t20;
t14 = pkin(5) + r_base(3);
t30 = t18 * pkin(1) + r_base(2);
t28 = t21 * pkin(1) + t18 * pkin(6) + r_base(1);
t25 = -t21 * pkin(6) + t30;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t14 - t42 * (t17 * pkin(2) + t14) + (t42 * qJ(3) + t44 * t16 + t43 * t19 - t49) * t20 + (t41 + t48) * t17) * g(3) + (-mrSges(1,2) - m(3) * t25 - m(5) * (pkin(7) * t35 + t30 + t45) - mrSges(5,3) * t35 + t47 * r_base(2) - t44 * (t17 * t37 - t33) - t43 * (t17 * t36 + t34) + t46 * (t25 + t45) + (-m(5) * (-pkin(3) - pkin(6)) + t39) * t21 + t40 * t18) * g(2) + (-m(3) * t28 - mrSges(1,1) + t47 * r_base(1) + t41 * t32 - t42 * (pkin(2) * t32 + t21 * t31 + t28) - t44 * (t17 * t34 + t36) - t43 * (t17 * t33 - t37) + (-m(5) * pkin(3) - t39) * t18 + t40 * t21) * g(1);
U = t1;
