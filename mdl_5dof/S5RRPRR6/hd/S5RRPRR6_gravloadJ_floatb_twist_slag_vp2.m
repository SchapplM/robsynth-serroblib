% Calculate Gravitation load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:41
% EndTime: 2019-12-05 18:35:43
% DurationCPUTime: 0.40s
% Computational Cost: add. (334->81), mult. (316->96), div. (0->0), fcn. (294->10), ass. (0->46)
t81 = mrSges(3,2) - mrSges(4,3);
t36 = sin(pkin(9));
t37 = cos(pkin(9));
t40 = cos(qJ(4));
t80 = mrSges(3,1) - m(5) * (-pkin(3) * t37 - pkin(2)) - m(6) * (-(t40 * pkin(4) + pkin(3)) * t37 - pkin(2)) + t37 * mrSges(4,1) + (m(5) * pkin(7) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3)) * t36;
t78 = m(6) * pkin(4) + mrSges(5,1);
t35 = qJ(1) + qJ(2);
t31 = sin(t35);
t63 = t36 * (-pkin(8) - pkin(7));
t33 = cos(t35);
t38 = sin(qJ(4));
t64 = t33 * t38;
t77 = pkin(4) * t64 + t31 * t63;
t34 = qJ(4) + qJ(5);
t30 = sin(t34);
t32 = cos(t34);
t65 = t33 * t37;
t11 = t30 * t65 - t31 * t32;
t12 = -t31 * t30 - t32 * t65;
t61 = t37 * t38;
t19 = -t31 * t40 + t33 * t61;
t60 = t37 * t40;
t66 = t31 * t38;
t20 = -t33 * t60 - t66;
t76 = -t20 * mrSges(5,1) - t12 * mrSges(6,1) - t19 * mrSges(5,2) - t11 * mrSges(6,2) - t81 * t31 + t80 * t33;
t67 = t31 * t37;
t10 = -t33 * t30 + t32 * t67;
t17 = t31 * t61 + t33 * t40;
t18 = t31 * t60 - t64;
t9 = t30 * t67 + t33 * t32;
t75 = t18 * mrSges(5,1) + t10 * mrSges(6,1) - t17 * mrSges(5,2) - t9 * mrSges(6,2) + t80 * t31 + t81 * t33;
t73 = t9 * mrSges(6,1) + t10 * mrSges(6,2);
t72 = -t11 * mrSges(6,1) + t12 * mrSges(6,2);
t71 = g(1) * t36;
t39 = sin(qJ(1));
t70 = t39 * pkin(1);
t41 = cos(qJ(1));
t69 = t41 * pkin(1);
t59 = t31 * qJ(3);
t26 = t33 * qJ(3);
t57 = t26 - t70;
t56 = -t31 * pkin(2) + t26;
t53 = -mrSges(6,1) * t30 - mrSges(6,2) * t32;
t52 = -t33 * pkin(2) - t59;
t50 = -pkin(4) * t66 + t33 * t63 - t59;
t1 = [(t39 * mrSges(2,1) + t41 * mrSges(2,2) + m(3) * t70 - m(4) * (t56 - t70) - m(5) * t57 - m(6) * (t57 + t77) + t75) * g(3) + (t41 * mrSges(2,1) - t39 * mrSges(2,2) + m(3) * t69 - m(4) * (t52 - t69) - m(5) * (-t59 - t69) - m(6) * (t50 - t69) + t76) * g(2), (-m(4) * t56 - m(5) * t26 - m(6) * (t26 + t77) + t75) * g(3) + (-m(4) * t52 + m(5) * t59 - m(6) * t50 + t76) * g(2), (g(2) * t33 + g(3) * t31) * (-m(4) - m(5) - m(6)), (mrSges(5,2) * t40 + t78 * t38 - t53) * t71 + (-t20 * mrSges(5,2) + t78 * t19 - t72) * g(3) + (-t18 * mrSges(5,2) - t78 * t17 - t73) * g(2), -g(2) * t73 - g(3) * t72 - t53 * t71];
taug = t1(:);
