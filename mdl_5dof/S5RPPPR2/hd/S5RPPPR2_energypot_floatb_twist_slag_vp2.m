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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:30:53
% EndTime: 2019-12-05 17:30:53
% DurationCPUTime: 0.75s
% Computational Cost: add. (176->86), mult. (325->94), div. (0->0), fcn. (348->10), ass. (0->38)
t59 = -m(1) - m(2);
t58 = -m(5) - m(6);
t26 = sin(pkin(7));
t29 = cos(pkin(7));
t57 = -t29 * mrSges(3,1) + t26 * mrSges(3,2) - mrSges(2,1);
t56 = mrSges(2,2) - mrSges(3,3);
t55 = -pkin(6) * m(6) + mrSges(5,2) - mrSges(6,3);
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t54 = t30 * mrSges(6,1) + t32 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t53 = -pkin(4) * m(6) - t32 * mrSges(6,1) + t30 * mrSges(6,2) - mrSges(5,1);
t25 = sin(pkin(8));
t51 = t25 * t26;
t28 = cos(pkin(8));
t50 = t26 * t28;
t31 = sin(qJ(1));
t49 = t26 * t31;
t33 = cos(qJ(1));
t48 = t26 * t33;
t47 = t29 * t33;
t46 = t31 * t25;
t45 = t31 * t28;
t44 = qJ(3) * t26;
t23 = pkin(5) + r_base(1);
t43 = t33 * qJ(2) + r_base(2);
t42 = t33 * pkin(1) + t31 * qJ(2) + r_base(3);
t40 = pkin(2) * t47 + t33 * t44 + t42;
t39 = -pkin(2) * t29 - pkin(1) - t44;
t37 = t26 * pkin(2) - qJ(3) * t29 + t23;
t35 = pkin(3) * t50 + qJ(4) * t51 + t37;
t27 = cos(pkin(9));
t24 = sin(pkin(9));
t12 = t28 * t47 + t46;
t11 = t25 * t47 - t45;
t10 = t25 * t33 - t29 * t45;
t9 = t28 * t33 + t29 * t46;
t8 = -t24 * t29 + t27 * t50;
t1 = (-m(3) * t42 - m(4) * t40 - t12 * mrSges(4,1) - mrSges(4,3) * t48 - mrSges(1,3) + t59 * r_base(3) + t58 * (t12 * pkin(3) + qJ(4) * t11 + t40) + t57 * t33 + t56 * t31 + t55 * (t12 * t24 - t27 * t48) + t53 * (t12 * t27 + t24 * t48) - t54 * t11) * g(3) + (-t10 * mrSges(4,1) - mrSges(1,2) + t59 * r_base(2) + t58 * (t10 * pkin(3) - qJ(4) * t9 + t39 * t31 + t43) + t53 * (t10 * t27 - t24 * t49) + t54 * t9 + (-m(3) - m(4)) * t43 + t56 * t33 + t55 * (t10 * t24 + t27 * t49) + (m(3) * pkin(1) - m(4) * t39 + t26 * mrSges(4,3) - t57) * t31) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - mrSges(2,3) - m(4) * t37 - m(5) * t35 - t8 * mrSges(5,1) - mrSges(5,3) * t51 - m(6) * (pkin(4) * t8 + t35) - (t30 * t51 + t32 * t8) * mrSges(6,1) - (-t30 * t8 + t32 * t51) * mrSges(6,2) + t55 * (t24 * t50 + t27 * t29) + (-mrSges(3,2) + mrSges(4,3)) * t29 + (-t28 * mrSges(4,1) + t25 * mrSges(4,2) - mrSges(3,1)) * t26 + (-m(2) - m(3)) * t23) * g(1);
U = t1;
