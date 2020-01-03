% Calculate potential energy for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:12
% EndTime: 2019-12-31 21:43:13
% DurationCPUTime: 0.52s
% Computational Cost: add. (214->83), mult. (425->88), div. (0->0), fcn. (481->10), ass. (0->48)
t62 = -m(1) - m(2);
t61 = -m(6) * pkin(9) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t27 = sin(qJ(5));
t31 = cos(qJ(5));
t60 = -t27 * mrSges(6,1) - t31 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t59 = m(6) * (pkin(4) + pkin(8)) + t31 * mrSges(6,1) - t27 * mrSges(6,2) + mrSges(5,1) + mrSges(4,3);
t58 = mrSges(3,2) - t59;
t26 = cos(pkin(5));
t33 = cos(qJ(2));
t34 = cos(qJ(1));
t46 = t34 * t33;
t29 = sin(qJ(2));
t30 = sin(qJ(1));
t49 = t30 * t29;
t12 = -t26 * t46 + t49;
t56 = t12 * pkin(8);
t47 = t34 * t29;
t48 = t30 * t33;
t14 = t26 * t48 + t47;
t55 = t14 * pkin(8);
t25 = sin(pkin(5));
t54 = t25 * t29;
t53 = t25 * t30;
t32 = cos(qJ(3));
t52 = t25 * t32;
t51 = t25 * t33;
t50 = t25 * t34;
t45 = pkin(6) + r_base(3);
t44 = pkin(8) * t51;
t43 = t26 * pkin(7) + t45;
t42 = t34 * pkin(1) + pkin(7) * t53 + r_base(1);
t41 = pkin(2) * t54 + t43;
t15 = -t26 * t49 + t46;
t40 = t15 * pkin(2) + t42;
t39 = t30 * pkin(1) - pkin(7) * t50 + r_base(2);
t13 = t26 * t47 + t48;
t38 = t13 * pkin(2) + t39;
t28 = sin(qJ(3));
t5 = t15 * t28 - t30 * t52;
t6 = t15 * t32 + t28 * t53;
t37 = t6 * pkin(3) + t5 * qJ(4) + t40;
t10 = -t26 * t32 + t28 * t54;
t11 = t26 * t28 + t29 * t52;
t36 = t11 * pkin(3) + t10 * qJ(4) + t41;
t3 = t13 * t28 + t32 * t50;
t4 = t13 * t32 - t28 * t50;
t35 = t4 * pkin(3) + t3 * qJ(4) + t38;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t45 - mrSges(2,3) - m(3) * t43 - t26 * mrSges(3,3) - (t29 * mrSges(3,1) + t33 * mrSges(3,2)) * t25 - m(4) * (t41 - t44) - m(5) * (t36 - t44) - m(6) * t36 + t59 * t51 + t60 * t10 + t61 * t11) * g(3) + (-mrSges(1,2) - t30 * mrSges(2,1) - t34 * mrSges(2,2) - m(3) * t39 - t13 * mrSges(3,1) + mrSges(3,3) * t50 - m(4) * (t38 + t56) - m(5) * (t35 + t56) - m(6) * t35 + t62 * r_base(2) + t61 * t4 + t60 * t3 + t58 * t12) * g(2) + (-mrSges(1,1) - t34 * mrSges(2,1) + t30 * mrSges(2,2) - m(3) * t42 - t15 * mrSges(3,1) - mrSges(3,3) * t53 - m(4) * (t40 + t55) - m(5) * (t37 + t55) - m(6) * t37 + t62 * r_base(1) + t61 * t6 + t60 * t5 + t58 * t14) * g(1);
U = t1;
