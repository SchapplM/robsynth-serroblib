% Calculate potential energy for
% S5PRPRR4
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:30
% EndTime: 2019-12-05 15:49:31
% DurationCPUTime: 0.70s
% Computational Cost: add. (271->75), mult. (561->84), div. (0->0), fcn. (668->12), ass. (0->39)
t42 = cos(qJ(2));
t73 = mrSges(3,2) * t42;
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t37 = sin(qJ(5));
t40 = cos(qJ(5));
t61 = -m(6) * pkin(4) - t40 * mrSges(6,1) + t37 * mrSges(6,2) - mrSges(5,1);
t65 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t70 = -m(5) - m(6);
t72 = t70 * pkin(3) + t65 * t38 + t61 * t41 - mrSges(4,1);
t33 = sin(pkin(5));
t36 = cos(pkin(5));
t39 = sin(qJ(2));
t55 = t36 * t39;
t69 = -mrSges(3,3) - mrSges(4,3);
t71 = -t55 * mrSges(3,1) - t36 * t73 - mrSges(2,2) + (m(3) * pkin(6) - t61 * t38 + t65 * t41 - t69) * t33;
t31 = sin(pkin(10));
t34 = cos(pkin(10));
t68 = t39 * t31 - t34 * t42;
t67 = -m(1) - m(2) - m(3);
t63 = t37 * mrSges(6,1) + t40 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t62 = -m(3) * pkin(1) - t42 * mrSges(3,1) + t39 * mrSges(3,2) - mrSges(2,1);
t52 = qJ(1) + r_base(3);
t19 = pkin(2) * t55 + (-pkin(6) - qJ(3)) * t33;
t28 = pkin(2) * t42 + pkin(1);
t32 = sin(pkin(9));
t35 = cos(pkin(9));
t51 = t35 * t19 + t32 * t28 + r_base(2);
t50 = t36 * pkin(6) + t52;
t49 = t31 * t42 + t39 * t34;
t48 = -t19 * t32 + t35 * t28 + r_base(1);
t47 = t33 * t39 * pkin(2) + t36 * qJ(3) + t50;
t46 = t68 * t36;
t18 = t49 * t36;
t17 = t49 * t33;
t16 = t68 * t33;
t9 = t32 * t46 - t35 * t49;
t7 = -t32 * t49 - t35 * t46;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t52 - mrSges(2,3) - m(3) * t50 - (mrSges(3,1) * t39 + t73) * t33 - m(4) * t47 - t17 * mrSges(4,1) + t70 * (t17 * pkin(3) + pkin(7) * t16 + t47) + t69 * t36 + t61 * (t17 * t41 + t36 * t38) - t63 * t16 + t65 * (t17 * t38 - t36 * t41)) * g(3) + (-m(4) * t51 - mrSges(1,2) + t62 * t32 + t67 * r_base(2) + t70 * (-pkin(7) * t7 + t51) + t63 * t7 + t72 * (t18 * t35 - t32 * t68) + t71 * t35) * g(2) + (-m(4) * t48 - mrSges(1,1) + t62 * t35 + t67 * r_base(1) + t70 * (-pkin(7) * t9 + t48) + t63 * t9 + t72 * (-t18 * t32 - t35 * t68) - t71 * t32) * g(1);
U = t1;
