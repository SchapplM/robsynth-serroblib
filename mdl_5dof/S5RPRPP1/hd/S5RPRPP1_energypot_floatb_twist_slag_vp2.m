% Calculate potential energy for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:35
% EndTime: 2019-12-31 18:08:35
% DurationCPUTime: 0.31s
% Computational Cost: add. (155->55), mult. (116->41), div. (0->0), fcn. (80->8), ass. (0->24)
t35 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t34 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t33 = -m(1) - m(2);
t32 = -m(3) - m(4);
t31 = -m(5) - m(6);
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t13 = qJ(3) + pkin(8);
t5 = sin(t13);
t7 = cos(t13);
t30 = -m(4) * pkin(2) - t18 * mrSges(4,1) + t16 * mrSges(4,2) - t34 * t5 + t35 * t7 - mrSges(3,1);
t29 = m(4) * pkin(6) - mrSges(3,2) + mrSges(6,2) + mrSges(4,3) + mrSges(5,3);
t28 = pkin(5) + r_base(3);
t17 = sin(qJ(1));
t27 = t17 * pkin(1) + r_base(2);
t19 = cos(qJ(1));
t26 = t19 * pkin(1) + r_base(1);
t9 = qJ(2) + t28;
t15 = -qJ(4) - pkin(6);
t14 = qJ(1) + pkin(7);
t8 = cos(t14);
t6 = sin(t14);
t4 = t18 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t28 - t16 * mrSges(4,1) - t18 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t32 * t9 + t31 * (t16 * pkin(3) + t9) + t34 * t7 + t35 * t5) * g(3) + (-t17 * mrSges(2,1) - t19 * mrSges(2,2) - mrSges(1,2) + t33 * r_base(2) + t32 * t27 + t31 * (t8 * t15 + t6 * t4 + t27) + t29 * t8 + t30 * t6) * g(2) + (-t19 * mrSges(2,1) + t17 * mrSges(2,2) - mrSges(1,1) + t33 * r_base(1) + t32 * t26 + t31 * (-t6 * t15 + t8 * t4 + t26) + t30 * t8 - t29 * t6) * g(1);
U = t1;
