% Calculate potential energy for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:16
% EndTime: 2019-12-05 15:40:16
% DurationCPUTime: 0.43s
% Computational Cost: add. (126->62), mult. (194->59), div. (0->0), fcn. (170->6), ass. (0->27)
t50 = mrSges(3,2) - mrSges(4,3);
t49 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t48 = -m(1) - m(2);
t47 = -m(5) - m(6);
t19 = sin(pkin(7));
t22 = sin(qJ(2));
t35 = qJ(3) * t22;
t24 = cos(qJ(2));
t41 = t19 * t24;
t46 = pkin(2) * t41 + t19 * t35;
t45 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t44 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t43 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t42 = t50 * t22 + t49 * t24 - mrSges(2,1);
t20 = cos(pkin(7));
t40 = t20 * t24;
t21 = sin(qJ(4));
t39 = t21 * t22;
t23 = cos(qJ(4));
t38 = t22 * t23;
t34 = t19 * pkin(1) + r_base(2);
t18 = qJ(1) + r_base(3);
t33 = t20 * pkin(1) + t19 * pkin(5) + r_base(1);
t32 = t22 * pkin(2) + t18;
t28 = -t20 * pkin(5) + t34;
t27 = pkin(2) * t40 + t20 * t35 + t33;
t1 = (-m(1) * r_base(3) - m(4) * t32 - mrSges(1,3) - mrSges(2,3) + t47 * (t22 * pkin(6) + t32) + (-m(2) - m(3)) * t18 + (-t43 * t23 + t44 * t21 + (m(4) - t47) * qJ(3) - t50) * t24 + t49 * t22) * g(3) + (-mrSges(1,2) - m(3) * t28 - m(4) * (t28 + t46) + t48 * r_base(2) + t47 * (pkin(6) * t41 + (-pkin(3) - pkin(5)) * t20 + t34 + t46) - t44 * (t19 * t39 - t20 * t23) + t43 * (t19 * t38 + t20 * t21) - t45 * t20 + t42 * t19) * g(2) + (-m(3) * t33 - m(4) * t27 - mrSges(1,1) + t48 * r_base(1) + t47 * (t19 * pkin(3) + pkin(6) * t40 + t27) - t44 * (t19 * t23 + t20 * t39) - t43 * (t19 * t21 - t20 * t38) + t45 * t19 + t42 * t20) * g(1);
U = t1;
