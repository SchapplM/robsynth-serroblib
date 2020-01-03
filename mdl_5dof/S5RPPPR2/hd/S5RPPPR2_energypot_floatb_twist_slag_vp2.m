% Calculate potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:15
% EndTime: 2020-01-03 11:22:15
% DurationCPUTime: 0.71s
% Computational Cost: add. (176->85), mult. (325->93), div. (0->0), fcn. (348->10), ass. (0->38)
t59 = -m(1) - m(2);
t58 = -m(5) - m(6);
t24 = sin(pkin(7));
t27 = cos(pkin(7));
t57 = -t27 * mrSges(3,1) + t24 * mrSges(3,2) - mrSges(2,1);
t56 = -mrSges(2,2) + mrSges(3,3);
t55 = pkin(2) * t27 + qJ(3) * t24;
t54 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t53 = t28 * mrSges(6,1) + t30 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t52 = -m(6) * pkin(4) - t30 * mrSges(6,1) + mrSges(6,2) * t28 - mrSges(5,1);
t23 = sin(pkin(8));
t49 = t23 * t24;
t26 = cos(pkin(8));
t48 = t24 * t26;
t29 = sin(qJ(1));
t47 = t24 * t29;
t31 = cos(qJ(1));
t46 = t24 * t31;
t45 = t27 * t31;
t44 = t29 * t23;
t43 = t29 * t26;
t21 = pkin(5) + r_base(1);
t41 = -t29 * qJ(2) + r_base(3);
t39 = t29 * pkin(1) - t31 * qJ(2) + r_base(2);
t38 = -pkin(1) - t55;
t37 = t24 * pkin(2) - qJ(3) * t27 + t21;
t36 = t29 * t55 + t39;
t33 = pkin(3) * t48 + qJ(4) * t49 + t37;
t25 = cos(pkin(9));
t22 = sin(pkin(9));
t12 = -t26 * t45 - t44;
t11 = t23 * t45 - t43;
t10 = -t23 * t31 + t27 * t43;
t9 = t26 * t31 + t27 * t44;
t8 = -t22 * t27 + t25 * t48;
t1 = (-t12 * mrSges(4,1) - mrSges(1,3) + t59 * r_base(3) + t58 * (t12 * pkin(3) - t11 * qJ(4) + t38 * t31 + t41) + (-m(3) - m(4)) * t41 + t54 * (t12 * t22 + t25 * t46) + t56 * t29 + t52 * (t12 * t25 - t22 * t46) + t53 * t11 + (m(3) * pkin(1) - m(4) * t38 + t24 * mrSges(4,3) - t57) * t31) * g(3) + (-m(3) * t39 - m(4) * t36 - t10 * mrSges(4,1) - mrSges(4,3) * t47 - mrSges(1,2) + t59 * r_base(2) + t58 * (t10 * pkin(3) + t9 * qJ(4) + t36) + t52 * (t10 * t25 + t22 * t47) - t53 * t9 + t56 * t31 + t57 * t29 + t54 * (t10 * t22 - t25 * t47)) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - mrSges(2,3) - m(4) * t37 - m(5) * t33 - t8 * mrSges(5,1) - mrSges(5,3) * t49 - m(6) * (pkin(4) * t8 + t33) - (t28 * t49 + t30 * t8) * mrSges(6,1) - (-t28 * t8 + t30 * t49) * mrSges(6,2) + t54 * (t22 * t48 + t25 * t27) + (-mrSges(3,2) + mrSges(4,3)) * t27 + (-t26 * mrSges(4,1) + t23 * mrSges(4,2) - mrSges(3,1)) * t24 + (-m(2) - m(3)) * t21) * g(1);
U = t1;
