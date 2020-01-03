% Calculate potential energy for
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:42:57
% EndTime: 2019-12-31 19:42:57
% DurationCPUTime: 0.54s
% Computational Cost: add. (149->67), mult. (254->65), div. (0->0), fcn. (250->8), ass. (0->31)
t51 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t50 = -t25 * mrSges(3,1) - mrSges(2,1) + (m(6) * pkin(7) - t51) * t22;
t49 = -m(1) - m(2);
t48 = -m(5) - m(6);
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t47 = (pkin(3) * t20 + qJ(4) * t19) * t22;
t45 = mrSges(2,2) - mrSges(3,3);
t21 = sin(qJ(5));
t24 = cos(qJ(5));
t42 = -t21 * mrSges(6,1) - t24 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t41 = -m(6) * pkin(4) - t24 * mrSges(6,1) + t21 * mrSges(6,2) - mrSges(4,1) - mrSges(5,1);
t23 = sin(qJ(1));
t39 = t23 * t25;
t26 = cos(qJ(1));
t37 = t26 * t25;
t36 = qJ(3) * t22;
t18 = pkin(5) + r_base(3);
t35 = t22 * pkin(2) + t18;
t34 = t26 * pkin(1) + t23 * pkin(6) + r_base(1);
t32 = t23 * pkin(1) - t26 * pkin(6) + r_base(2);
t31 = pkin(2) * t37 + t26 * t36 + t34;
t30 = -t25 * qJ(3) + t35;
t29 = pkin(2) * t39 + t23 * t36 + t32;
t6 = t23 * t19 + t20 * t37;
t5 = t19 * t37 - t23 * t20;
t4 = -t26 * t19 + t20 * t39;
t3 = t19 * t39 + t26 * t20;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t30 - m(5) * (t30 + t47) - m(6) * (t35 + t47) + (-m(2) - m(3)) * t18 + (-m(6) * (pkin(7) - qJ(3)) + t51) * t25 + (t42 * t19 + t41 * t20 - mrSges(3,1)) * t22) * g(3) + (-m(3) * t32 - m(4) * t29 - mrSges(1,2) + t49 * r_base(2) + t48 * (t4 * pkin(3) + t3 * qJ(4) + t29) + t41 * t4 + t42 * t3 - t45 * t26 + t50 * t23) * g(2) + (-m(3) * t34 - m(4) * t31 - mrSges(1,1) + t49 * r_base(1) + t48 * (t6 * pkin(3) + t5 * qJ(4) + t31) + t41 * t6 + t42 * t5 + t45 * t23 + t50 * t26) * g(1);
U = t1;
