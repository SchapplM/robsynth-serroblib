% Calculate potential energy for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:27
% EndTime: 2019-12-31 17:59:28
% DurationCPUTime: 0.32s
% Computational Cost: add. (141->55), mult. (119->41), div. (0->0), fcn. (87->8), ass. (0->22)
t12 = sin(qJ(5));
t15 = cos(qJ(5));
t35 = -m(6) * pkin(4) - t15 * mrSges(6,1) + t12 * mrSges(6,2) - mrSges(5,1);
t34 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t33 = -m(1) - m(2);
t32 = -m(5) - m(6);
t30 = -t12 * mrSges(6,1) - t15 * mrSges(6,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t13 = sin(qJ(4));
t16 = cos(qJ(4));
t29 = t35 * t13 - t34 * t16 + mrSges(3,2) - mrSges(4,3);
t28 = pkin(5) + r_base(3);
t14 = sin(qJ(1));
t27 = t14 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t26 = t17 * pkin(1) + r_base(1);
t11 = qJ(1) + pkin(8);
t6 = sin(t11);
t25 = t6 * pkin(2) + t27;
t8 = qJ(2) + t28;
t7 = cos(t11);
t23 = t7 * pkin(2) + t6 * qJ(3) + t26;
t1 = (-m(1) * r_base(3) - m(2) * t28 - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t8 + t32 * (pkin(3) + t8) + t35 * t16 + t34 * t13) * g(3) + (-m(3) * t27 - m(4) * t25 - t14 * mrSges(2,1) - mrSges(2,2) * t17 - mrSges(1,2) + t33 * r_base(2) + t32 * (t6 * pkin(6) + t25) + ((m(4) - t32) * qJ(3) - t29) * t7 + t30 * t6) * g(2) + (-m(3) * t26 - m(4) * t23 - mrSges(2,1) * t17 + t14 * mrSges(2,2) - mrSges(1,1) + t33 * r_base(1) + t32 * (t7 * pkin(6) + t23) + t30 * t7 + t29 * t6) * g(1);
U = t1;
