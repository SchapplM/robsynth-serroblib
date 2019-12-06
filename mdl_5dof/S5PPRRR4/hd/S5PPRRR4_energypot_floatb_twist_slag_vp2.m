% Calculate potential energy for
% S5PPRRR4
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:31
% EndTime: 2019-12-05 15:18:31
% DurationCPUTime: 0.67s
% Computational Cost: add. (374->97), mult. (901->128), div. (0->0), fcn. (1108->14), ass. (0->51)
t42 = sin(pkin(6));
t46 = cos(pkin(6));
t47 = cos(pkin(5));
t43 = sin(pkin(5));
t44 = cos(pkin(11));
t73 = t43 * t44;
t25 = -t42 * t73 + t46 * t47;
t40 = sin(pkin(11));
t45 = cos(pkin(10));
t41 = sin(pkin(10));
t74 = t41 * t47;
t28 = -t40 * t45 - t44 * t74;
t71 = t43 * t46;
t19 = -t28 * t42 + t41 * t71;
t83 = -m(1) - m(2);
t82 = -m(5) - m(6);
t81 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t48 = sin(qJ(5));
t51 = cos(qJ(5));
t80 = -mrSges(6,1) * t48 - mrSges(6,2) * t51 + mrSges(4,2) - mrSges(5,3);
t79 = -m(6) * pkin(4) - t51 * mrSges(6,1) + t48 * mrSges(6,2) - mrSges(5,1);
t78 = cos(qJ(3));
t76 = t40 * t43;
t75 = t41 * t43;
t72 = t43 * t45;
t70 = t45 * t47;
t68 = qJ(2) * t43;
t65 = qJ(1) + r_base(3);
t64 = t42 * t78;
t63 = t46 * t78;
t62 = t45 * pkin(1) + t41 * t68 + r_base(1);
t61 = t47 * qJ(2) + t65;
t60 = t43 * t64;
t26 = -t40 * t41 + t44 * t70;
t18 = -t26 * t42 - t45 * t71;
t59 = t41 * pkin(1) - t45 * t68 + r_base(2);
t29 = -t40 * t74 + t44 * t45;
t58 = t29 * pkin(2) + t19 * pkin(7) + t62;
t57 = pkin(2) * t76 + t25 * pkin(7) + t61;
t27 = t40 * t70 + t41 * t44;
t54 = t27 * pkin(2) + pkin(7) * t18 + t59;
t52 = cos(qJ(4));
t50 = sin(qJ(3));
t49 = sin(qJ(4));
t17 = t47 * t42 * t50 + (t44 * t46 * t50 + t40 * t78) * t43;
t16 = -t47 * t64 + t50 * t76 - t63 * t73;
t10 = t29 * t78 + (t28 * t46 + t42 * t75) * t50;
t9 = -t28 * t63 + t29 * t50 - t41 * t60;
t8 = t27 * t78 + (t26 * t46 - t42 * t72) * t50;
t7 = -t26 * t63 + t27 * t50 + t45 * t60;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t65 - mrSges(2,3) - m(3) * t61 - t47 * mrSges(3,3) - (t40 * mrSges(3,1) + t44 * mrSges(3,2)) * t43 - m(4) * t57 - t17 * mrSges(4,1) - t25 * mrSges(4,3) + t82 * (t17 * pkin(3) + pkin(8) * t16 + t57) + t79 * (t17 * t52 + t25 * t49) + t80 * t16 + t81 * (t17 * t49 - t25 * t52)) * g(3) + (-m(3) * t59 - m(4) * t54 - t41 * mrSges(2,1) - t27 * mrSges(3,1) - t8 * mrSges(4,1) - t45 * mrSges(2,2) - t26 * mrSges(3,2) + mrSges(3,3) * t72 - t18 * mrSges(4,3) - mrSges(1,2) + t83 * r_base(2) + t82 * (t8 * pkin(3) + pkin(8) * t7 + t54) + t79 * (t18 * t49 + t52 * t8) + t80 * t7 + t81 * (-t18 * t52 + t49 * t8)) * g(2) + (-m(3) * t62 - m(4) * t58 - t45 * mrSges(2,1) - t29 * mrSges(3,1) - t10 * mrSges(4,1) + t41 * mrSges(2,2) - t28 * mrSges(3,2) - mrSges(3,3) * t75 - t19 * mrSges(4,3) - mrSges(1,1) + t83 * r_base(1) + t82 * (t10 * pkin(3) + pkin(8) * t9 + t58) + t79 * (t10 * t52 + t19 * t49) + t80 * t9 + t81 * (t10 * t49 - t19 * t52)) * g(1);
U = t1;
