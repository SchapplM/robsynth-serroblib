% Calculate potential energy for
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:16
% EndTime: 2019-12-31 18:12:16
% DurationCPUTime: 0.34s
% Computational Cost: add. (143->57), mult. (140->42), div. (0->0), fcn. (104->6), ass. (0->23)
t13 = pkin(7) + qJ(3);
t10 = sin(t13);
t11 = cos(t13);
t41 = pkin(3) * t11 + qJ(4) * t10;
t40 = -m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t39 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t38 = -m(2) - m(3);
t37 = -m(5) - m(6);
t19 = cos(qJ(1));
t36 = t41 * t19;
t35 = -m(1) + t38;
t15 = sin(pkin(7));
t16 = cos(pkin(7));
t33 = -m(3) * pkin(1) - t16 * mrSges(3,1) + t15 * mrSges(3,2) + t39 * t10 + t40 * t11 - mrSges(2,1);
t32 = -m(3) * qJ(2) - mrSges(5,1) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t14 = pkin(5) + r_base(3);
t8 = t16 * pkin(2) + pkin(1);
t29 = t19 * t8 + r_base(1);
t27 = t15 * pkin(2) + t14;
t17 = -pkin(6) - qJ(2);
t18 = sin(qJ(1));
t23 = -t18 * t17 + t29;
t1 = (-m(1) * r_base(3) - m(4) * t27 - t15 * mrSges(3,1) - t16 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t37 * (t10 * pkin(3) - t11 * qJ(4) + t27) + t38 * t14 - t39 * t11 + t40 * t10) * g(3) + (-mrSges(1,2) + t35 * r_base(2) + (m(6) * pkin(4) - t32) * t19 + (-m(4) + t37) * (t19 * t17 + t18 * t8 + r_base(2)) + (t37 * t41 + t33) * t18) * g(2) + (-mrSges(1,1) - m(4) * t23 - m(5) * (t23 + t36) - m(6) * (t29 + t36) + t35 * r_base(1) + t33 * t19 + (-m(6) * (pkin(4) - t17) + t32) * t18) * g(1);
U = t1;
