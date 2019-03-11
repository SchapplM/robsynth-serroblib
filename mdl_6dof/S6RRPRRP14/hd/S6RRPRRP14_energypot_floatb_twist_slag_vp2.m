% Calculate potential energy for
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP14_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:39
% EndTime: 2019-03-09 13:03:39
% DurationCPUTime: 0.65s
% Computational Cost: add. (296->99), mult. (596->108), div. (0->0), fcn. (686->10), ass. (0->49)
t76 = -m(1) - m(2);
t75 = -m(6) - m(7);
t74 = -mrSges(3,1) + mrSges(4,2);
t73 = mrSges(3,3) + mrSges(4,1);
t72 = -mrSges(4,3) + mrSges(3,2);
t71 = -mrSges(5,3) + t74;
t70 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t69 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t68 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t36 = sin(pkin(6));
t40 = sin(qJ(2));
t67 = t36 * t40;
t41 = sin(qJ(1));
t66 = t36 * t41;
t44 = cos(qJ(2));
t65 = t36 * t44;
t45 = cos(qJ(1));
t64 = t36 * t45;
t63 = t40 * t41;
t62 = t40 * t45;
t61 = t41 * t44;
t60 = t44 * t45;
t59 = pkin(7) + r_base(3);
t58 = pkin(8) * t64;
t57 = t41 * pkin(1) + r_base(2);
t37 = cos(pkin(6));
t56 = t37 * pkin(8) + t59;
t55 = t45 * pkin(1) + pkin(8) * t66 + r_base(1);
t54 = pkin(2) * t67 + t56;
t21 = -t37 * t60 + t63;
t22 = t37 * t62 + t61;
t53 = t22 * pkin(2) + t21 * qJ(3) + t57;
t23 = t37 * t61 + t62;
t24 = -t37 * t63 + t60;
t52 = t24 * pkin(2) + qJ(3) * t23 + t55;
t51 = t37 * pkin(3) + pkin(9) * t67 - qJ(3) * t65 + t54;
t50 = pkin(3) * t66 + pkin(9) * t24 + t52;
t48 = t22 * pkin(9) + (-pkin(3) - pkin(8)) * t64 + t53;
t43 = cos(qJ(4));
t42 = cos(qJ(5));
t39 = sin(qJ(4));
t38 = sin(qJ(5));
t20 = t37 * t43 - t39 * t65;
t19 = t37 * t39 + t43 * t65;
t12 = t21 * t39 - t43 * t64;
t11 = t21 * t43 + t39 * t64;
t10 = t23 * t39 + t43 * t66;
t9 = -t23 * t43 + t39 * t66;
t1 = (-m(1) * r_base(3) - m(2) * t59 - m(3) * t56 - m(4) * t54 - m(5) * t51 - t20 * mrSges(5,1) - mrSges(5,3) * t67 - mrSges(1,3) - mrSges(2,3) + t75 * (t20 * pkin(4) + pkin(10) * t19 + t51) + t69 * (t20 * t42 + t38 * t67) + t68 * (t20 * t38 - t42 * t67) - t73 * t37 + ((m(4) * qJ(3) - t72) * t44 + t74 * t40) * t36 + t70 * t19) * g(3) + (-m(3) * (t57 - t58) - m(4) * (t53 - t58) - m(5) * t48 - mrSges(1,2) - t12 * mrSges(5,1) - t41 * mrSges(2,1) - t45 * mrSges(2,2) + t76 * r_base(2) + t73 * t64 + t75 * (t12 * pkin(4) - t11 * pkin(10) + t48) + t69 * (t12 * t42 + t22 * t38) + t68 * (t12 * t38 - t22 * t42) + t72 * t21 + t71 * t22 - t70 * t11) * g(2) + (-m(3) * t55 - m(4) * t52 - m(5) * t50 - t45 * mrSges(2,1) - t10 * mrSges(5,1) + t41 * mrSges(2,2) - mrSges(1,1) + t76 * r_base(1) - t73 * t66 + t75 * (t10 * pkin(4) + pkin(10) * t9 + t50) + t72 * t23 + t69 * (t10 * t42 + t24 * t38) + t68 * (t10 * t38 - t24 * t42) + t70 * t9 + t71 * t24) * g(1);
U  = t1;
