% Calculate potential energy for
% S5RRRRR12
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:29
% EndTime: 2019-12-31 22:46:30
% DurationCPUTime: 0.63s
% Computational Cost: add. (374->97), mult. (901->124), div. (0->0), fcn. (1108->14), ass. (0->51)
t43 = cos(pkin(5));
t48 = sin(qJ(1));
t51 = cos(qJ(2));
t69 = t48 * t51;
t47 = sin(qJ(2));
t52 = cos(qJ(1));
t71 = t47 * t52;
t28 = -t43 * t69 - t71;
t40 = sin(pkin(6));
t42 = cos(pkin(6));
t41 = sin(pkin(5));
t75 = t41 * t48;
t19 = -t28 * t40 + t42 * t75;
t74 = t41 * t51;
t25 = -t40 * t74 + t42 * t43;
t83 = -m(1) - m(2);
t82 = -m(5) - m(6);
t81 = -m(6) * pkin(11) + mrSges(5,2) - mrSges(6,3);
t44 = sin(qJ(5));
t49 = cos(qJ(5));
t80 = -mrSges(6,1) * t44 - mrSges(6,2) * t49 + mrSges(4,2) - mrSges(5,3);
t79 = -m(6) * pkin(4) - t49 * mrSges(6,1) + t44 * mrSges(6,2) - mrSges(5,1);
t78 = cos(qJ(3));
t76 = t41 * t47;
t73 = t41 * t52;
t70 = t48 * t47;
t68 = t51 * t52;
t67 = pkin(7) + r_base(3);
t64 = t40 * t78;
t63 = t42 * t78;
t62 = t43 * pkin(8) + t67;
t61 = t52 * pkin(1) + pkin(8) * t75 + r_base(1);
t60 = t41 * t64;
t26 = t43 * t68 - t70;
t18 = -t26 * t40 - t42 * t73;
t59 = t48 * pkin(1) - pkin(8) * t73 + r_base(2);
t29 = -t43 * t70 + t68;
t58 = t29 * pkin(2) + t19 * pkin(9) + t61;
t57 = pkin(2) * t76 + t25 * pkin(9) + t62;
t27 = t43 * t71 + t69;
t54 = t27 * pkin(2) + pkin(9) * t18 + t59;
t50 = cos(qJ(4));
t46 = sin(qJ(3));
t45 = sin(qJ(4));
t17 = t43 * t40 * t46 + (t42 * t46 * t51 + t47 * t78) * t41;
t16 = -t43 * t64 + t46 * t76 - t63 * t74;
t12 = t29 * t78 + (t28 * t42 + t40 * t75) * t46;
t11 = -t28 * t63 + t29 * t46 - t48 * t60;
t10 = t27 * t78 + (t26 * t42 - t40 * t73) * t46;
t9 = -t26 * t63 + t27 * t46 + t52 * t60;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t67 - mrSges(2,3) - m(3) * t62 - t43 * mrSges(3,3) - (mrSges(3,1) * t47 + mrSges(3,2) * t51) * t41 - m(4) * t57 - t17 * mrSges(4,1) - t25 * mrSges(4,3) + t82 * (t17 * pkin(3) + pkin(10) * t16 + t57) + t81 * (t17 * t45 - t25 * t50) + t79 * (t17 * t50 + t25 * t45) + t80 * t16) * g(3) + (-m(3) * t59 - m(4) * t54 - t48 * mrSges(2,1) - t27 * mrSges(3,1) - t10 * mrSges(4,1) - t52 * mrSges(2,2) - t26 * mrSges(3,2) + mrSges(3,3) * t73 - t18 * mrSges(4,3) - mrSges(1,2) + t83 * r_base(2) + t82 * (t10 * pkin(3) + t9 * pkin(10) + t54) + t79 * (t10 * t50 + t18 * t45) + t80 * t9 + t81 * (t10 * t45 - t18 * t50)) * g(2) + (-m(3) * t61 - m(4) * t58 - t52 * mrSges(2,1) - t29 * mrSges(3,1) - t12 * mrSges(4,1) + t48 * mrSges(2,2) - t28 * mrSges(3,2) - mrSges(3,3) * t75 - t19 * mrSges(4,3) - mrSges(1,1) + t83 * r_base(1) + t82 * (t12 * pkin(3) + pkin(10) * t11 + t58) + t81 * (t12 * t45 - t19 * t50) + t79 * (t12 * t50 + t19 * t45) + t80 * t11) * g(1);
U = t1;
