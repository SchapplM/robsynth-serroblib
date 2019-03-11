% Calculate potential energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:23
% EndTime: 2019-03-09 15:14:23
% DurationCPUTime: 0.64s
% Computational Cost: add. (313->101), mult. (682->118), div. (0->0), fcn. (786->10), ass. (0->50)
t79 = -m(1) - m(2);
t78 = -m(6) - m(7);
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t77 = -t45 * mrSges(3,1) + mrSges(3,2) * t42 - mrSges(2,1);
t76 = mrSges(2,2) - mrSges(3,3);
t75 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t74 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1) - mrSges(5,3);
t73 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t39 = sin(pkin(6));
t72 = t39 * t45;
t41 = sin(qJ(3));
t71 = t41 * t42;
t43 = sin(qJ(1));
t70 = t42 * t43;
t46 = cos(qJ(1));
t69 = t42 * t46;
t68 = t43 * t45;
t67 = t45 * t46;
t66 = qJ(4) * t39;
t40 = cos(pkin(6));
t65 = qJ(4) * t40;
t64 = cos(pkin(10));
t37 = pkin(7) + r_base(3);
t63 = t39 * t71;
t62 = t42 * t65;
t61 = pkin(2) * t42 + t37;
t60 = t40 * t64;
t59 = t42 * t64;
t58 = pkin(1) * t46 + pkin(8) * t43 + r_base(1);
t57 = t39 * t59;
t55 = pkin(1) * t43 - pkin(8) * t46 + r_base(2);
t54 = pkin(2) * t67 + pkin(9) * t69 + t58;
t53 = pkin(2) * t68 + pkin(9) * t70 + t55;
t44 = cos(qJ(3));
t52 = qJ(4) * t63 + t42 * t44 * pkin(3) + (-pkin(9) - t65) * t45 + t61;
t19 = -t41 * t67 + t43 * t44;
t20 = t41 * t43 + t44 * t67;
t51 = pkin(3) * t20 - t19 * t66 + t46 * t62 + t54;
t17 = -t41 * t68 - t44 * t46;
t18 = -t41 * t46 + t44 * t68;
t50 = pkin(3) * t18 - t17 * t66 + t43 * t62 + t53;
t38 = sin(pkin(10));
t9 = t44 * t59 + (-t40 * t71 - t72) * t38;
t8 = t64 * t72 + (t38 * t44 + t41 * t60) * t42;
t6 = t20 * t64 + (t19 * t40 + t39 * t69) * t38;
t5 = -t19 * t60 + t20 * t38 - t46 * t57;
t4 = t18 * t64 + (t17 * t40 + t39 * t70) * t38;
t3 = -t17 * t60 + t18 * t38 - t43 * t57;
t1 = (-m(1) * r_base(3) - m(4) * t61 - m(5) * t52 - mrSges(1,3) - mrSges(2,3) + t78 * (pkin(4) * t9 + qJ(5) * t8 + t52) + (m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3)) * t45 + (-mrSges(4,1) * t44 + mrSges(4,2) * t41 - mrSges(3,1)) * t42 + (-m(2) - m(3)) * t37 + t73 * t9 + t75 * t8 + t74 * (-t40 * t45 + t63)) * g(3) + (-m(3) * t55 - m(4) * t53 - m(5) * t50 - t18 * mrSges(4,1) - t17 * mrSges(4,2) - mrSges(4,3) * t70 - mrSges(1,2) + t79 * r_base(2) + t78 * (pkin(4) * t4 + qJ(5) * t3 + t50) - t76 * t46 + t77 * t43 + t73 * t4 + t75 * t3 + t74 * (-t17 * t39 + t40 * t70)) * g(2) + (-m(3) * t58 - m(4) * t54 - m(5) * t51 - t20 * mrSges(4,1) - t19 * mrSges(4,2) - mrSges(4,3) * t69 - mrSges(1,1) + t79 * r_base(1) + t78 * (pkin(4) * t6 + qJ(5) * t5 + t51) + t77 * t46 + t76 * t43 + t73 * t6 + t75 * t5 + t74 * (-t19 * t39 + t40 * t69)) * g(1);
U  = t1;
