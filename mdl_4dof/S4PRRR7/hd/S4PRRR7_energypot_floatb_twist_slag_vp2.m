% Calculate potential energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:00
% EndTime: 2019-12-31 16:36:00
% DurationCPUTime: 0.40s
% Computational Cost: add. (145->71), mult. (282->86), div. (0->0), fcn. (311->10), ass. (0->33)
t49 = -m(1) - m(2);
t48 = -m(4) - m(5);
t47 = -m(5) * pkin(7) + mrSges(4,2) - mrSges(5,3);
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t46 = -t25 * mrSges(5,1) - t28 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t45 = -m(5) * pkin(3) - t28 * mrSges(5,1) + t25 * mrSges(5,2) - mrSges(4,1);
t21 = sin(pkin(8));
t22 = sin(pkin(4));
t44 = t21 * t22;
t27 = sin(qJ(2));
t43 = t22 * t27;
t29 = cos(qJ(3));
t42 = t22 * t29;
t30 = cos(qJ(2));
t41 = t22 * t30;
t23 = cos(pkin(8));
t40 = t23 * t22;
t24 = cos(pkin(4));
t39 = t24 * t27;
t38 = t24 * t30;
t37 = qJ(1) + r_base(3);
t36 = t23 * pkin(1) + pkin(5) * t44 + r_base(1);
t35 = t24 * pkin(5) + t37;
t34 = t21 * pkin(1) - pkin(5) * t40 + r_base(2);
t32 = pkin(2) * t43 - pkin(6) * t41 + t35;
t26 = sin(qJ(3));
t12 = t24 * t26 + t27 * t42;
t10 = -t21 * t39 + t23 * t30;
t9 = t21 * t38 + t23 * t27;
t8 = t21 * t30 + t23 * t39;
t7 = t21 * t27 - t23 * t38;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t37 - mrSges(2,3) - m(3) * t35 - t24 * mrSges(3,3) - (t27 * mrSges(3,1) + t30 * mrSges(3,2)) * t22 - m(4) * t32 - t12 * mrSges(4,1) + mrSges(4,3) * t41 - m(5) * (t12 * pkin(3) + t32) - (t12 * t28 - t25 * t41) * mrSges(5,1) - (-t12 * t25 - t28 * t41) * mrSges(5,2) + t47 * (-t24 * t29 + t26 * t43)) * g(3) + (-m(3) * t34 - t21 * mrSges(2,1) - t8 * mrSges(3,1) - t23 * mrSges(2,2) + mrSges(3,3) * t40 - mrSges(1,2) + t49 * r_base(2) + t48 * (t8 * pkin(2) + t7 * pkin(6) + t34) + t45 * (-t26 * t40 + t8 * t29) + t46 * t7 + t47 * (t8 * t26 + t29 * t40)) * g(2) + (-m(3) * t36 - t23 * mrSges(2,1) - t10 * mrSges(3,1) + t21 * mrSges(2,2) - mrSges(3,3) * t44 - mrSges(1,1) + t49 * r_base(1) + t48 * (t10 * pkin(2) + t9 * pkin(6) + t36) + t45 * (t10 * t29 + t26 * t44) + t46 * t9 + t47 * (t10 * t26 - t21 * t42)) * g(1);
U = t1;
