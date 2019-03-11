% Calculate potential energy for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:05
% EndTime: 2019-03-09 04:18:06
% DurationCPUTime: 0.50s
% Computational Cost: add. (159->75), mult. (206->57), div. (0->0), fcn. (172->8), ass. (0->26)
t46 = -mrSges(4,1) + mrSges(5,2);
t45 = -m(1) - m(2);
t44 = mrSges(6,3) + mrSges(7,3);
t43 = m(5) + m(6) + m(7);
t15 = sin(qJ(5));
t18 = cos(qJ(5));
t14 = qJ(5) + qJ(6);
t4 = sin(t14);
t5 = cos(t14);
t42 = -t4 * mrSges(7,1) - t18 * mrSges(6,2) - t5 * mrSges(7,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t15;
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t41 = -t19 * mrSges(4,2) + t46 * t16 + mrSges(2,2) - mrSges(3,3);
t40 = -m(6) * pkin(4) - t18 * mrSges(6,1) + t15 * mrSges(6,2) - mrSges(2,1) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t18 + pkin(4)) - t5 * mrSges(7,1) + t4 * mrSges(7,2);
t21 = -pkin(9) - pkin(8);
t39 = -m(6) * pkin(8) + m(7) * t21 - t44;
t38 = pkin(3) * t16;
t13 = pkin(6) + r_base(3);
t17 = sin(qJ(1));
t35 = t17 * pkin(1) + r_base(2);
t34 = pkin(2) + t13;
t33 = t17 * pkin(7) + t35;
t20 = cos(qJ(1));
t32 = t20 * pkin(1) + t17 * qJ(2) + r_base(1);
t30 = t20 * pkin(7) + t32;
t1 = (-m(1) * r_base(3) - m(4) * t34 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t13 - t43 * (t19 * pkin(3) + t16 * qJ(4) + t34) + (t39 + t46) * t19 + (mrSges(4,2) + t42) * t16) * g(3) + (-m(3) * t35 - m(4) * t33 - mrSges(1,2) + t45 * r_base(2) - t43 * (t20 * t19 * qJ(4) + t33) + (m(5) * t38 + (-m(6) * (-pkin(3) - pkin(8)) - m(7) * (-pkin(3) + t21) + t44) * t16 + t42 * t19 + (m(3) + m(4) + t43) * qJ(2) - t41) * t20 + t40 * t17) * g(2) + (-m(3) * t32 - m(4) * t30 - mrSges(1,1) + t45 * r_base(1) - t43 * (t17 * t38 + t30) + t40 * t20 + (t39 * t16 + (t43 * qJ(4) - t42) * t19 + t41) * t17) * g(1);
U  = t1;
