% Calculate potential energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:09
% EndTime: 2019-12-05 17:23:10
% DurationCPUTime: 0.66s
% Computational Cost: add. (374->97), mult. (901->127), div. (0->0), fcn. (1108->14), ass. (0->50)
t41 = sin(pkin(6));
t44 = cos(pkin(6));
t45 = cos(pkin(5));
t42 = sin(pkin(5));
t52 = cos(qJ(2));
t71 = t42 * t52;
t25 = -t41 * t71 + t45 * t44;
t40 = sin(pkin(11));
t43 = cos(pkin(11));
t49 = sin(qJ(2));
t68 = t45 * t52;
t28 = -t40 * t68 - t43 * t49;
t73 = t42 * t44;
t19 = -t28 * t41 + t40 * t73;
t82 = -m(1) - m(2);
t81 = -m(5) - m(6);
t80 = -m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t46 = sin(qJ(5));
t50 = cos(qJ(5));
t79 = -mrSges(6,1) * t46 - mrSges(6,2) * t50 + mrSges(4,2) - mrSges(5,3);
t78 = -m(6) * pkin(4) - t50 * mrSges(6,1) + t46 * mrSges(6,2) - mrSges(5,1);
t77 = cos(qJ(3));
t75 = t40 * t42;
t74 = t42 * t43;
t72 = t42 * t49;
t69 = t45 * t49;
t65 = qJ(1) + r_base(3);
t64 = t41 * t77;
t63 = t44 * t77;
t62 = t43 * pkin(1) + pkin(7) * t75 + r_base(1);
t61 = t45 * pkin(7) + t65;
t60 = t42 * t64;
t26 = -t40 * t49 + t43 * t68;
t18 = -t26 * t41 - t43 * t73;
t59 = t40 * pkin(1) - pkin(7) * t74 + r_base(2);
t29 = -t40 * t69 + t43 * t52;
t58 = t29 * pkin(2) + t19 * pkin(8) + t62;
t57 = pkin(2) * t72 + t25 * pkin(8) + t61;
t27 = t40 * t52 + t43 * t69;
t54 = t27 * pkin(2) + pkin(8) * t18 + t59;
t51 = cos(qJ(4));
t48 = sin(qJ(3));
t47 = sin(qJ(4));
t17 = t45 * t41 * t48 + (t44 * t48 * t52 + t49 * t77) * t42;
t16 = -t45 * t64 + t48 * t72 - t63 * t71;
t10 = t29 * t77 + (t28 * t44 + t41 * t75) * t48;
t9 = -t28 * t63 + t29 * t48 - t40 * t60;
t8 = t27 * t77 + (t26 * t44 - t41 * t74) * t48;
t7 = -t26 * t63 + t27 * t48 + t43 * t60;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t65 - mrSges(2,3) - m(3) * t61 - t45 * mrSges(3,3) - (t49 * mrSges(3,1) + t52 * mrSges(3,2)) * t42 - m(4) * t57 - t17 * mrSges(4,1) - t25 * mrSges(4,3) + t81 * (t17 * pkin(3) + t16 * pkin(9) + t57) + t78 * (t17 * t51 + t25 * t47) + t79 * t16 + t80 * (t17 * t47 - t25 * t51)) * g(3) + (-m(3) * t59 - m(4) * t54 - t40 * mrSges(2,1) - t27 * mrSges(3,1) - t8 * mrSges(4,1) - t43 * mrSges(2,2) - t26 * mrSges(3,2) + mrSges(3,3) * t74 - t18 * mrSges(4,3) - mrSges(1,2) + t82 * r_base(2) + t81 * (t8 * pkin(3) + pkin(9) * t7 + t54) + t78 * (t18 * t47 + t51 * t8) + t79 * t7 + t80 * (-t18 * t51 + t47 * t8)) * g(2) + (-m(3) * t62 - m(4) * t58 - t43 * mrSges(2,1) - t29 * mrSges(3,1) - t10 * mrSges(4,1) + t40 * mrSges(2,2) - t28 * mrSges(3,2) - mrSges(3,3) * t75 - t19 * mrSges(4,3) - mrSges(1,1) + t82 * r_base(1) + t81 * (t10 * pkin(3) + pkin(9) * t9 + t58) + t78 * (t10 * t51 + t19 * t47) + t79 * t9 + t80 * (t10 * t47 - t19 * t51)) * g(1);
U = t1;
