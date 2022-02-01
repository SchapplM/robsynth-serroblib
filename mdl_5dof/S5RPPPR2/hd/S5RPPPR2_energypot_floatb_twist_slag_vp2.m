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
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 08:58:59
% EndTime: 2022-01-23 08:58:59
% DurationCPUTime: 0.60s
% Computational Cost: add. (176->96), mult. (275->109), div. (0->0), fcn. (280->10), ass. (0->39)
t58 = -m(2) - m(6);
t57 = -m(3) - m(4);
t31 = cos(qJ(5));
t56 = t31 * mrSges(6,1);
t55 = mrSges(4,2) - mrSges(5,3);
t54 = mrSges(5,2) - mrSges(6,3);
t53 = -m(1) - m(5) + t58;
t23 = sin(pkin(9));
t26 = cos(pkin(9));
t18 = t26 * pkin(4) + t23 * pkin(6) + pkin(3);
t24 = sin(pkin(8));
t27 = cos(pkin(8));
t38 = qJ(4) * t27 - qJ(2);
t29 = sin(qJ(5));
t48 = t24 * t29;
t52 = m(5) * (-t24 * pkin(3) + t38) + m(6) * (-t18 * t24 + t38) + (t26 * t48 + t31 * t27) * mrSges(6,2) + mrSges(2,2) - mrSges(3,3);
t49 = t24 * qJ(4) + pkin(2);
t17 = t27 * pkin(3) + t49;
t25 = sin(pkin(7));
t28 = cos(pkin(7));
t35 = t23 * pkin(4) - t26 * pkin(6) + qJ(3);
t47 = t24 * t31;
t50 = t25 * qJ(3) + pkin(1);
t6 = t18 * t27 + t49;
t45 = t28 * t26;
t9 = t25 * t23 + t27 * t45;
t51 = -m(4) * (pkin(2) * t28 + t50) - m(5) * (t17 * t28 + t50) - m(6) * (t25 * t35 + t6 * t28 + pkin(1)) - (t28 * t47 - t9 * t29) * mrSges(6,2) - mrSges(2,1) - m(3) * pkin(1) - t28 * mrSges(3,1) + t25 * mrSges(3,2);
t46 = t25 * t27;
t30 = sin(qJ(1));
t44 = t30 * t25;
t43 = t30 * t28;
t32 = cos(qJ(1));
t42 = t32 * t24;
t41 = t32 * t25;
t40 = t32 * t28;
t22 = pkin(5) + r_base(3);
t14 = t30 * t24 + t27 * t40;
t12 = t27 * t43 - t42;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t54 * (t23 * t46 + t45) + (-m(4) - m(5)) * (-t28 * qJ(3) + t22) + (m(6) * t35 - mrSges(3,2) + mrSges(4,3)) * t28 + (-m(3) + t58) * t22 + (t29 * mrSges(6,2) - mrSges(5,1) - t56) * (-t28 * t23 + t26 * t46) + (-m(4) * pkin(2) - m(5) * t17 - m(6) * t6 - t27 * mrSges(4,1) - t48 * mrSges(6,1) - t47 * mrSges(6,2) + t55 * t24 - mrSges(3,1)) * t25) * g(3) + (-mrSges(1,2) - t12 * mrSges(4,1) - mrSges(4,3) * t44 - (t12 * t26 + t23 * t44) * mrSges(5,1) + t26 * t42 * t56 + t57 * (-t32 * qJ(2) + r_base(2)) + t53 * r_base(2) + (-t9 * t56 + t51) * t30 - t52 * t32 + t54 * (t12 * t23 - t26 * t44) + (-t29 * mrSges(6,1) + t55) * (t24 * t43 + t32 * t27)) * g(2) + (-mrSges(1,1) - t14 * mrSges(4,1) - mrSges(4,3) * t41 - (t14 * t26 + t23 * t41) * mrSges(5,1) + t54 * (t14 * t23 - t26 * t41) + t57 * (t30 * qJ(2) + r_base(1)) + t53 * r_base(1) + (-(t28 * t48 + t9 * t31) * mrSges(6,1) + t51) * t32 + (-(t26 * t47 - t29 * t27) * mrSges(6,1) + t52) * t30 + t55 * (t24 * t40 - t30 * t27)) * g(1);
U = t1;
