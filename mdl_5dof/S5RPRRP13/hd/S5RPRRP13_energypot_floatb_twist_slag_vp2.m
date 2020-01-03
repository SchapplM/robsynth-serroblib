% Calculate potential energy for
% S5RPRRP13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:27
% EndTime: 2019-12-31 18:58:27
% DurationCPUTime: 0.42s
% Computational Cost: add. (116->60), mult. (168->54), div. (0->0), fcn. (144->6), ass. (0->25)
t45 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t44 = -m(1) - m(2);
t43 = -m(5) - m(6);
t42 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t41 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t40 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t17 = sin(qJ(3));
t20 = cos(qJ(3));
t39 = -t17 * mrSges(4,1) - t20 * t45 + mrSges(2,2) - mrSges(3,3);
t38 = pkin(7) * t20;
t16 = sin(qJ(4));
t21 = cos(qJ(1));
t37 = t16 * t21;
t18 = sin(qJ(1));
t36 = t18 * t17;
t19 = cos(qJ(4));
t35 = t18 * t19;
t32 = t21 * t19;
t15 = pkin(5) + r_base(3);
t31 = t18 * pkin(1) + r_base(2);
t30 = pkin(2) + t15;
t29 = t21 * pkin(1) + t18 * qJ(2) + r_base(1);
t27 = t18 * pkin(6) + t31;
t25 = t21 * pkin(6) + t29;
t1 = (-m(1) * r_base(3) - m(4) * t30 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t43 * (t20 * pkin(3) + t17 * pkin(7) + t30) + (-m(2) - m(3)) * t15 + (-t16 * t40 + t19 * t41 - mrSges(4,1)) * t20 + t45 * t17) * g(3) + (-m(3) * t31 - m(4) * t27 - mrSges(1,2) + t44 * r_base(2) + t43 * (t21 * t38 + t27) + t41 * (t16 * t18 - t17 * t32) + t40 * (t17 * t37 + t35) + t42 * t18 + (t43 * (-pkin(3) * t17 - qJ(2)) + (m(3) + m(4)) * qJ(2) - t39) * t21) * g(2) + (-m(3) * t29 - m(4) * t25 - mrSges(1,1) + t44 * r_base(1) + t43 * (pkin(3) * t36 - t18 * t38 + t25) + t41 * (t17 * t35 + t37) - t40 * (t16 * t36 - t32) + t42 * t21 + t39 * t18) * g(1);
U = t1;
