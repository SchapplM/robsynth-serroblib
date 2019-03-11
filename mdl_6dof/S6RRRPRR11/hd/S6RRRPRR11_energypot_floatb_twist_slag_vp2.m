% Calculate potential energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:43
% EndTime: 2019-03-09 19:24:44
% DurationCPUTime: 0.74s
% Computational Cost: add. (355->101), mult. (751->113), div. (0->0), fcn. (900->12), ass. (0->50)
t79 = -m(6) - m(7);
t38 = sin(qJ(6));
t43 = cos(qJ(6));
t81 = -t38 * mrSges(7,1) - t43 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t80 = -m(1) - m(2);
t78 = -mrSges(4,1) - mrSges(5,1);
t77 = mrSges(4,2) - mrSges(5,3);
t76 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t74 = -m(7) * pkin(5) - t43 * mrSges(7,1) + t38 * mrSges(7,2) - mrSges(6,1);
t73 = mrSges(3,2) - t81 + t79 * (pkin(9) - pkin(10));
t41 = sin(qJ(2));
t45 = cos(qJ(2));
t46 = cos(qJ(1));
t42 = sin(qJ(1));
t64 = cos(pkin(6));
t60 = t42 * t64;
t26 = t46 * t41 + t45 * t60;
t71 = pkin(9) * t26;
t59 = t46 * t64;
t24 = t41 * t42 - t45 * t59;
t70 = t24 * pkin(9);
t69 = cos(qJ(3));
t37 = sin(pkin(6));
t68 = t37 * t41;
t67 = t37 * t42;
t66 = t37 * t45;
t65 = t37 * t46;
t63 = pkin(7) + r_base(3);
t62 = t37 * t69;
t61 = t64 * pkin(8) + t63;
t58 = t46 * pkin(1) + pkin(8) * t67 + r_base(1);
t27 = -t41 * t60 + t46 * t45;
t56 = t27 * pkin(2) + t58;
t55 = t42 * pkin(1) - pkin(8) * t65 + r_base(2);
t25 = t41 * t59 + t42 * t45;
t54 = t25 * pkin(2) + t55;
t53 = pkin(2) * t68 - pkin(9) * t66 + t61;
t40 = sin(qJ(3));
t15 = t27 * t40 - t42 * t62;
t16 = t27 * t69 + t40 * t67;
t52 = t16 * pkin(3) + qJ(4) * t15 + t56;
t13 = t25 * t40 + t46 * t62;
t14 = t25 * t69 - t40 * t65;
t50 = t14 * pkin(3) + t13 * qJ(4) + t54;
t22 = t40 * t68 - t64 * t69;
t23 = t64 * t40 + t41 * t62;
t49 = t23 * pkin(3) + qJ(4) * t22 + t53;
t44 = cos(qJ(5));
t39 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t63 - mrSges(2,3) - m(3) * t61 - t64 * mrSges(3,3) - (t41 * mrSges(3,1) + t45 * mrSges(3,2)) * t37 - m(4) * t53 - m(5) * t49 + t79 * (t23 * pkin(4) + pkin(10) * t66 + t49) + t76 * (-t22 * t44 + t23 * t39) + t78 * t23 + t77 * t22 + t74 * (t22 * t39 + t23 * t44) + t81 * t66) * g(3) + (-m(4) * (t54 + t70) - m(5) * (t50 + t70) + mrSges(3,3) * t65 - m(3) * t55 - mrSges(1,2) - t25 * mrSges(3,1) - t42 * mrSges(2,1) - t46 * mrSges(2,2) + t80 * r_base(2) + t79 * (t14 * pkin(4) + t50) + t78 * t14 + t77 * t13 + t76 * (-t13 * t44 + t14 * t39) + t74 * (t13 * t39 + t14 * t44) + t73 * t24) * g(2) + (-m(4) * (t56 + t71) - m(5) * (t52 + t71) - m(3) * t58 - mrSges(1,1) - mrSges(3,3) * t67 - t27 * mrSges(3,1) + t42 * mrSges(2,2) - t46 * mrSges(2,1) + t80 * r_base(1) + t79 * (t16 * pkin(4) + t52) + t76 * (-t15 * t44 + t16 * t39) + t78 * t16 + t77 * t15 + t74 * (t15 * t39 + t16 * t44) + t73 * t26) * g(1);
U  = t1;
