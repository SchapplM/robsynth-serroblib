% Calculate Gravitation load on the joints for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:07
% EndTime: 2019-12-05 16:54:10
% DurationCPUTime: 0.69s
% Computational Cost: add. (311->77), mult. (774->119), div. (0->0), fcn. (897->10), ass. (0->46)
t38 = cos(qJ(4));
t84 = -m(6) * (t38 * pkin(4) + pkin(3)) - mrSges(4,1);
t83 = m(6) * (-qJ(5) - pkin(8)) + mrSges(4,2) - mrSges(6,3);
t82 = mrSges(5,1) + mrSges(6,1);
t77 = -mrSges(5,2) - mrSges(6,2);
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t81 = t83 * t36 + t84 * t39 - mrSges(3,1);
t70 = m(6) * pkin(4);
t78 = mrSges(3,2) - mrSges(4,3);
t76 = -m(4) - m(5) - m(6);
t35 = sin(qJ(4));
t75 = m(5) * pkin(3) + t77 * t35 + t82 * t38 - t84;
t74 = -t70 - t82;
t73 = -m(5) * pkin(8) - mrSges(5,3) + t83;
t72 = m(5) * (pkin(3) * t39 + pkin(8) * t36) + t36 * mrSges(5,3);
t71 = t72 - t81;
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t57 = cos(pkin(9));
t58 = cos(pkin(5));
t46 = t58 * t57;
t56 = sin(pkin(9));
t20 = t37 * t46 + t56 * t40;
t69 = t20 * t35;
t45 = t58 * t56;
t22 = -t37 * t45 + t57 * t40;
t68 = t22 * t35;
t33 = sin(pkin(5));
t67 = t33 * t37;
t66 = t33 * t40;
t65 = t35 * t37;
t64 = t35 * t39;
t61 = t38 * t39;
t60 = t39 * t40;
t51 = t33 * t57;
t50 = t33 * t56;
t24 = t58 * t36 + t39 * t67;
t23 = t36 * t67 - t58 * t39;
t21 = t57 * t37 + t40 * t45;
t19 = t56 * t37 - t40 * t46;
t12 = t22 * t39 + t36 * t50;
t11 = t22 * t36 - t39 * t50;
t10 = t20 * t39 - t36 * t51;
t9 = t20 * t36 + t39 * t51;
t1 = [(-m(2) - m(3) + t76) * g(3), (-t69 * t70 - t82 * (-t19 * t61 + t69) + t76 * (-t19 * pkin(2) + t20 * pkin(7)) + t77 * (t19 * t64 + t20 * t38) + t78 * t20 + t71 * t19) * g(2) + (-t68 * t70 - t82 * (-t21 * t61 + t68) + t77 * (t21 * t64 + t22 * t38) + t76 * (-t21 * pkin(2) + t22 * pkin(7)) + t78 * t22 + t71 * t21) * g(1) + (-t72 * t66 + t76 * (pkin(2) * t66 + pkin(7) * t67) + (-t82 * (t38 * t60 + t65) + t77 * (-t35 * t60 + t37 * t38) - t65 * t70 + t78 * t37 + t81 * t40) * t33) * g(3), (t75 * t23 + t73 * t24) * g(3) + (t73 * t10 + t75 * t9) * g(2) + (t75 * t11 + t73 * t12) * g(1), (t77 * (-t24 * t38 + t35 * t66) + t74 * (-t24 * t35 - t38 * t66)) * g(3) + (t77 * (-t10 * t38 - t19 * t35) + t74 * (-t10 * t35 + t19 * t38)) * g(2) + (t77 * (-t12 * t38 - t21 * t35) + t74 * (-t12 * t35 + t21 * t38)) * g(1), (-g(1) * t11 - g(2) * t9 - g(3) * t23) * m(6)];
taug = t1(:);
