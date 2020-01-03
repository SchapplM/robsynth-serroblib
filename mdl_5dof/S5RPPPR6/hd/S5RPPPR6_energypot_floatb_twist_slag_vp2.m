% Calculate potential energy for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:32
% EndTime: 2019-12-31 17:47:33
% DurationCPUTime: 0.56s
% Computational Cost: add. (134->80), mult. (215->80), div. (0->0), fcn. (199->8), ass. (0->35)
t52 = -mrSges(3,1) + mrSges(4,2);
t51 = mrSges(3,2) - mrSges(4,3);
t50 = -m(1) - m(2);
t49 = -m(5) - m(6);
t24 = sin(qJ(1));
t20 = sin(pkin(7));
t38 = qJ(3) * t20;
t22 = cos(pkin(7));
t44 = t22 * t24;
t48 = pkin(2) * t44 + t24 * t38;
t47 = t51 * t20 + t52 * t22 - mrSges(2,1);
t46 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t45 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t26 = cos(qJ(1));
t43 = t22 * t26;
t19 = sin(pkin(8));
t42 = t24 * t19;
t21 = cos(pkin(8));
t41 = t24 * t21;
t40 = t26 * t19;
t39 = t26 * t21;
t37 = qJ(4) * t22;
t18 = pkin(5) + r_base(3);
t36 = t24 * pkin(1) + r_base(2);
t35 = t20 * pkin(2) + t18;
t34 = t26 * pkin(1) + t24 * qJ(2) + r_base(1);
t30 = pkin(2) * t43 + t26 * t38 + t34;
t29 = -t26 * qJ(2) + t36;
t28 = t24 * pkin(3) + t26 * t37 + t30;
t27 = t24 * t37 + (-pkin(3) - qJ(2)) * t26 + t36 + t48;
t25 = cos(qJ(5));
t23 = sin(qJ(5));
t4 = t20 * t42 - t39;
t2 = t20 * t40 + t41;
t1 = (-m(1) * r_base(3) - m(4) * t35 - mrSges(1,3) - mrSges(2,3) + t49 * (t20 * qJ(4) + t35) + (-m(2) - m(3)) * t18 + ((m(6) * pkin(4) + t25 * mrSges(6,1) - t23 * mrSges(6,2) + mrSges(5,1)) * t19 - t45 * t21 + (m(4) - t49) * qJ(3) - t51) * t22 + (-t23 * mrSges(6,1) - t25 * mrSges(6,2) - mrSges(5,3) + t52) * t20) * g(3) + (-mrSges(1,2) - m(3) * t29 - m(4) * (t29 + t48) - m(5) * t27 - t4 * mrSges(5,1) - mrSges(5,3) * t44 - m(6) * (t4 * pkin(4) + t27) - (t23 * t44 + t4 * t25) * mrSges(6,1) - (-t4 * t23 + t25 * t44) * mrSges(6,2) + t50 * r_base(2) + t45 * (t20 * t41 + t40) - t46 * t26 + t47 * t24) * g(2) + (-mrSges(1,1) - m(3) * t34 - m(4) * t30 - m(5) * t28 - t2 * mrSges(5,1) - mrSges(5,3) * t43 - m(6) * (t2 * pkin(4) + t28) - (t2 * t25 + t23 * t43) * mrSges(6,1) - (-t2 * t23 + t25 * t43) * mrSges(6,2) + t50 * r_base(1) - t45 * (-t20 * t39 + t42) + t47 * t26 + t46 * t24) * g(1);
U = t1;
