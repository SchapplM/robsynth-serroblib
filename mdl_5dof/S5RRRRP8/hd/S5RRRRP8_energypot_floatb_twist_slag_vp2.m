% Calculate potential energy for
% S5RRRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:34
% EndTime: 2019-12-31 21:59:34
% DurationCPUTime: 0.43s
% Computational Cost: add. (160->62), mult. (192->55), div. (0->0), fcn. (168->8), ass. (0->25)
t17 = sin(qJ(3));
t20 = cos(qJ(3));
t7 = t20 * pkin(3) + pkin(2);
t16 = qJ(3) + qJ(4);
t9 = cos(t16);
t45 = -m(5) * t7 - m(6) * (pkin(4) * t9 + t7) - mrSges(3,1) - m(4) * pkin(2) - t20 * mrSges(4,1) + t17 * mrSges(4,2);
t23 = -pkin(8) - pkin(7);
t44 = m(5) * t23 + m(6) * (-qJ(5) + t23) + mrSges(3,2) - mrSges(5,3) - mrSges(6,3) - m(4) * pkin(7) - mrSges(4,3);
t43 = -m(1) - m(2);
t42 = -mrSges(5,1) - mrSges(6,1);
t41 = mrSges(5,2) + mrSges(6,2);
t40 = -m(3) - m(5) - m(6);
t39 = -m(4) + t40;
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t37 = t44 * t18 + t45 * t21 - mrSges(2,1);
t35 = t17 * pkin(3);
t8 = sin(t16);
t36 = -m(5) * t35 - m(6) * (pkin(4) * t8 + t35) + mrSges(2,2) - mrSges(3,3) - t17 * mrSges(4,1) - t20 * mrSges(4,2);
t19 = sin(qJ(1));
t34 = t19 * t21;
t22 = cos(qJ(1));
t33 = t21 * t22;
t32 = t19 * pkin(1) + r_base(2);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) + t39) * (pkin(5) + r_base(3)) - t44 * t21 + (t41 * t8 + t42 * t9 + t45) * t18) * g(3) + (-m(4) * t32 - mrSges(1,2) + t43 * r_base(2) + t42 * (-t22 * t8 + t9 * t34) - t41 * (-t22 * t9 - t8 * t34) + t40 * (-t22 * pkin(6) + t32) + (m(4) * pkin(6) - t36) * t22 + t37 * t19) * g(2) + (-mrSges(1,1) + t43 * r_base(1) + t42 * (t19 * t8 + t9 * t33) - t41 * (t19 * t9 - t8 * t33) + t39 * (t22 * pkin(1) + t19 * pkin(6) + r_base(1)) + t36 * t19 + t37 * t22) * g(1);
U = t1;
