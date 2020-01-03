% Calculate potential energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:36:58
% EndTime: 2019-12-31 21:36:59
% DurationCPUTime: 0.64s
% Computational Cost: add. (237->92), mult. (447->104), div. (0->0), fcn. (511->12), ass. (0->43)
t60 = -m(1) - m(2);
t59 = -m(4) - m(5);
t58 = -m(5) * qJ(4) + m(6) * (-pkin(9) - qJ(4)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t24 = pkin(10) + qJ(5);
t19 = sin(t24);
t20 = cos(t24);
t25 = sin(pkin(10));
t27 = cos(pkin(10));
t45 = pkin(4) * t25 + pkin(8);
t57 = -m(6) * t45 - t25 * mrSges(5,1) - t19 * mrSges(6,1) - t27 * mrSges(5,2) - t20 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t18 = pkin(4) * t27 + pkin(3);
t56 = -m(5) * pkin(3) - m(6) * t18 - t27 * mrSges(5,1) - t20 * mrSges(6,1) + t25 * mrSges(5,2) + t19 * mrSges(6,2) - mrSges(4,1);
t26 = sin(pkin(5));
t31 = sin(qJ(2));
t55 = t26 * t31;
t32 = sin(qJ(1));
t54 = t26 * t32;
t33 = cos(qJ(3));
t53 = t26 * t33;
t34 = cos(qJ(2));
t52 = t26 * t34;
t35 = cos(qJ(1));
t51 = t26 * t35;
t50 = t31 * t32;
t49 = t31 * t35;
t48 = t32 * t34;
t47 = t34 * t35;
t46 = pkin(6) + r_base(3);
t28 = cos(pkin(5));
t44 = t28 * pkin(7) + t46;
t43 = t35 * pkin(1) + pkin(7) * t54 + r_base(1);
t42 = pkin(2) * t55 + t44;
t12 = -t28 * t50 + t47;
t41 = t12 * pkin(2) + t43;
t40 = t32 * pkin(1) - pkin(7) * t51 + r_base(2);
t10 = t28 * t49 + t48;
t39 = t10 * pkin(2) + t40;
t37 = -pkin(8) * t52 + t42;
t30 = sin(qJ(3));
t11 = t28 * t48 + t49;
t9 = -t28 * t47 + t50;
t8 = t28 * t30 + t31 * t53;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t46 - mrSges(2,3) - m(3) * t44 - t28 * mrSges(3,3) - (t31 * mrSges(3,1) + t34 * mrSges(3,2)) * t26 - m(4) * t37 - t8 * mrSges(4,1) + mrSges(4,3) * t52 - m(5) * (pkin(3) * t8 + t37) - (-t25 * t52 + t27 * t8) * mrSges(5,1) - (-t25 * t8 - t27 * t52) * mrSges(5,2) - m(6) * (t18 * t8 - t45 * t52 + t42) - (-t19 * t52 + t20 * t8) * mrSges(6,1) - (-t19 * t8 - t20 * t52) * mrSges(6,2) + t58 * (-t28 * t33 + t30 * t55)) * g(3) + (-m(3) * t40 - m(6) * t39 - t32 * mrSges(2,1) - t10 * mrSges(3,1) - t35 * mrSges(2,2) + mrSges(3,3) * t51 - mrSges(1,2) + t60 * r_base(2) + t59 * (t9 * pkin(8) + t39) + t56 * (t10 * t33 - t30 * t51) + t57 * t9 + t58 * (t10 * t30 + t33 * t51)) * g(2) + (-m(3) * t43 - m(6) * t41 - t35 * mrSges(2,1) - t12 * mrSges(3,1) + t32 * mrSges(2,2) - mrSges(3,3) * t54 - mrSges(1,1) + t60 * r_base(1) + t59 * (pkin(8) * t11 + t41) + t56 * (t12 * t33 + t30 * t54) + t57 * t11 + t58 * (t12 * t30 - t32 * t53)) * g(1);
U = t1;
