% Calculate Gravitation load on the joints for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:14
% EndTime: 2024-09-27 21:45:14
% DurationCPUTime: 0.52s
% Computational Cost: add. (542->99), mult. (370->128), div. (0->0), fcn. (322->20), ass. (0->64)
t87 = -m(5) * pkin(3) - mrSges(4,1);
t52 = qJ(3) + qJ(4);
t46 = pkin(5) + t52;
t36 = cos(t46) / 0.2e1;
t67 = pkin(5) - t52;
t38 = cos(t67);
t65 = sin(t67);
t68 = sin(t46) / 0.2e1;
t39 = qJ(5) + t46;
t33 = cos(t39) / 0.2e1;
t66 = -qJ(5) + t67;
t37 = cos(t66);
t62 = sin(t66);
t69 = sin(t39) / 0.2e1;
t71 = (t69 + t62 / 0.2e1) * mrSges(6,1) + (t33 - t37 / 0.2e1) * mrSges(6,2);
t86 = -(t68 + t65 / 0.2e1) * mrSges(5,1) - (t36 - t38 / 0.2e1) * mrSges(5,2) - t71;
t27 = t38 / 0.2e1 + t36;
t51 = pkin(10) + qJ(2);
t44 = sin(t51);
t45 = cos(t51);
t47 = sin(t52);
t75 = t45 * t47;
t12 = -t44 * t27 - t75;
t26 = t68 - t65 / 0.2e1;
t48 = cos(t52);
t22 = t69 - t62 / 0.2e1;
t49 = qJ(5) + t52;
t41 = cos(t49);
t23 = t37 / 0.2e1 + t33;
t40 = sin(t49);
t6 = -t44 * t23 - t45 * t40;
t80 = t6 * mrSges(6,1) + (t44 * t22 - t45 * t41) * mrSges(6,2);
t85 = -t12 * mrSges(5,1) - (t44 * t26 - t45 * t48) * mrSges(5,2) - t80;
t76 = t44 * t47;
t11 = -t45 * t27 + t76;
t5 = -t45 * t23 + t44 * t40;
t81 = -t5 * mrSges(6,1) + (-t45 * t22 - t44 * t41) * mrSges(6,2);
t84 = t11 * mrSges(5,1) - (-t45 * t26 - t44 * t48) * mrSges(5,2) - t81;
t82 = pkin(7) + pkin(8);
t53 = sin(pkin(5));
t78 = m(6) * t53;
t58 = cos(qJ(3));
t50 = t58 * pkin(3);
t43 = t50 + pkin(2);
t54 = cos(pkin(5));
t56 = sin(qJ(3));
t74 = t54 * t56;
t73 = t54 * t58;
t55 = sin(qJ(4));
t72 = t55 * t56;
t57 = cos(qJ(4));
t42 = t57 * pkin(4) + pkin(3);
t64 = -pkin(4) * t72 + t42 * t58;
t63 = -t44 * t56 + t45 * t73;
t18 = -t44 * t73 - t45 * t56;
t61 = pkin(4) * (t57 * t58 - t72);
t60 = m(4) * pkin(2) + m(5) * t43 + m(6) * (pkin(4) * t48 + t43) + t48 * mrSges(5,1) + t41 * mrSges(6,1) + mrSges(3,1);
t59 = m(5) * (pkin(3) * t74 - t53 * t82) + m(6) * (-t53 * (pkin(9) + t82) + (t58 * t55 * pkin(4) + t56 * t42) * t54) + t26 * mrSges(5,1) + t22 * mrSges(6,1) + mrSges(3,2) + (-m(4) * pkin(7) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3)) * t53;
t30 = -t56 * pkin(3) - pkin(4) * t47;
t19 = -t44 * t74 + t45 * t58;
t17 = -t44 * t58 - t45 * t74;
t15 = t54 * t61;
t14 = t64 * t54;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), (-t19 * mrSges(4,1) - t18 * mrSges(4,2) - t12 * mrSges(5,2) - t6 * mrSges(6,2) + t59 * t44 - t60 * t45) * g(2) + (-t17 * mrSges(4,1) + mrSges(4,2) * t63 - t11 * mrSges(5,2) - t5 * mrSges(6,2) + t60 * t44 + t59 * t45) * g(1), (-t64 * t78 + (-m(5) * t50 - mrSges(4,1) * t58 + mrSges(4,2) * t56) * t53 + t86) * g(3) + (-t17 * mrSges(4,2) - m(6) * (t45 * t14 + t44 * t30) + t84 + t87 * t63) * g(2) + (t19 * mrSges(4,2) - m(6) * (-t44 * t14 + t45 * t30) + t87 * t18 + t85) * g(1), (-t61 * t78 + t86) * g(3) + (-m(6) * (-pkin(4) * t76 + t45 * t15) + t84) * g(2) + (-m(6) * (-pkin(4) * t75 - t44 * t15) + t85) * g(1), -g(1) * t80 - g(2) * t81 - g(3) * t71];
taug = t1(:);
