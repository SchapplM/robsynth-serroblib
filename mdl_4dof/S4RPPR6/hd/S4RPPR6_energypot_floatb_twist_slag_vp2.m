% Calculate potential energy for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:34
% EndTime: 2019-12-31 16:40:35
% DurationCPUTime: 0.32s
% Computational Cost: add. (85->48), mult. (125->41), div. (0->0), fcn. (101->6), ass. (0->20)
t36 = mrSges(3,2) - mrSges(4,3);
t10 = sin(pkin(6));
t11 = cos(pkin(6));
t35 = pkin(2) * t11 + qJ(3) * t10;
t34 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t32 = -m(1) - m(2);
t31 = -m(4) - m(5);
t13 = sin(qJ(1));
t30 = t35 * t13;
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t19 = t10 * t12 + t11 * t14;
t20 = t10 * t14 - t11 * t12;
t29 = -t19 * mrSges(5,1) - t20 * mrSges(5,2) + t36 * t10 + t34 * t11 - mrSges(2,1);
t28 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3);
t9 = pkin(4) + r_base(3);
t25 = t13 * pkin(1) + r_base(2);
t15 = cos(qJ(1));
t18 = -qJ(2) * t15 + t25;
t1 = (-m(1) * r_base(3) - t20 * mrSges(5,1) + t19 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t9 + t31 * (t10 * pkin(2) - qJ(3) * t11 + t9) - t36 * t11 + t34 * t10) * g(3) + (-mrSges(1,2) - m(3) * t18 - m(4) * (t18 + t30) - m(5) * (t25 + t30) + t32 * r_base(2) + (-m(5) * (pkin(5) - qJ(2)) + t28) * t15 + t29 * t13) * g(2) + (-mrSges(1,1) + t32 * r_base(1) + (m(5) * pkin(5) - t28) * t13 + (-m(3) + t31) * (t15 * pkin(1) + t13 * qJ(2) + r_base(1)) + (t31 * t35 + t29) * t15) * g(1);
U = t1;
