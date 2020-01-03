% Calculate potential energy for
% S5RPRRP12
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:11
% EndTime: 2019-12-31 18:56:12
% DurationCPUTime: 0.39s
% Computational Cost: add. (114->57), mult. (155->48), div. (0->0), fcn. (127->6), ass. (0->23)
t16 = cos(qJ(4));
t41 = m(5) * pkin(3) + m(6) * (t16 * pkin(4) + pkin(3)) + mrSges(4,1);
t40 = -m(5) * pkin(7) + m(6) * (-qJ(5) - pkin(7)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t39 = m(6) * pkin(4);
t38 = -m(1) - m(2);
t37 = -mrSges(5,1) - mrSges(6,1);
t36 = mrSges(5,2) + mrSges(6,2);
t35 = -m(4) - m(5) - m(6);
t34 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t14 = sin(qJ(3));
t17 = cos(qJ(3));
t33 = t41 * t14 + t40 * t17 - mrSges(2,2) + mrSges(3,3);
t13 = sin(qJ(4));
t15 = sin(qJ(1));
t32 = t15 * t13;
t31 = t15 * t16;
t18 = cos(qJ(1));
t28 = t18 * t13;
t27 = t18 * t16;
t11 = pkin(5) + r_base(3);
t26 = t15 * pkin(1) + r_base(2);
t23 = t18 * pkin(1) + t15 * qJ(2) + r_base(1);
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t11 + t35 * (pkin(2) + t11) + (t36 * t13 + t37 * t16 - t41) * t17 + t40 * t14) * g(3) + (-t32 * t39 - m(3) * t26 - mrSges(1,2) + t38 * r_base(2) + t37 * (-t14 * t27 + t32) - t36 * (t14 * t28 + t31) + t35 * (t15 * pkin(6) + t26) + t34 * t15 + ((m(3) - t35) * qJ(2) + t33) * t18) * g(2) + (-t28 * t39 - m(3) * t23 - mrSges(1,1) + t38 * r_base(1) + t35 * (t18 * pkin(6) + t23) + t37 * (t14 * t31 + t28) - t36 * (-t14 * t32 + t27) + t34 * t18 - t33 * t15) * g(1);
U = t1;
