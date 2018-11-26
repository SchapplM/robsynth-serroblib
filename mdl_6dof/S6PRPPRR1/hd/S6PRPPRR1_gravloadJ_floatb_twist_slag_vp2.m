% Calculate Gravitation load on the joints for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:53:28
% EndTime: 2018-11-23 14:53:29
% DurationCPUTime: 0.69s
% Computational Cost: add. (980->96), mult. (762->128), div. (0->0), fcn. (674->22), ass. (0->60)
t56 = sin(qJ(6));
t58 = cos(qJ(6));
t97 = m(7) * pkin(5) + mrSges(7,1) * t58 - mrSges(7,2) * t56 + mrSges(6,1);
t92 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t48 = pkin(6) - qJ(2);
t34 = sin(t48) / 0.2e1;
t47 = pkin(6) + qJ(2);
t39 = sin(t47);
t96 = t39 / 0.2e1 + t34;
t95 = -m(4) - m(5);
t94 = -m(6) - m(7);
t35 = cos(t47) / 0.2e1;
t44 = cos(t48);
t27 = t44 / 0.2e1 + t35;
t91 = m(5) - t94;
t45 = pkin(12) + qJ(5);
t37 = sin(t45);
t41 = cos(t45);
t52 = cos(pkin(12));
t90 = m(5) * pkin(3) + t52 * mrSges(5,1) - sin(pkin(12)) * mrSges(5,2) + mrSges(4,1) + t97 * t41 - t92 * t37;
t89 = -m(5) * qJ(4) - t56 * mrSges(7,1) - t58 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t50 = sin(pkin(10));
t51 = sin(pkin(6));
t86 = t50 * t51;
t57 = sin(qJ(2));
t85 = t50 * t57;
t53 = cos(pkin(10));
t84 = t51 * t53;
t83 = t53 * t57;
t82 = t96 * pkin(2);
t46 = qJ(2) + pkin(11);
t81 = m(4) + t91;
t78 = pkin(6) - t46;
t77 = pkin(6) + t46;
t25 = t27 * pkin(2);
t76 = -pkin(2) * t85 + t53 * t25;
t66 = sin(t77) / 0.2e1;
t70 = sin(t78);
t22 = t66 - t70 / 0.2e1;
t42 = cos(t46);
t73 = t53 * t22 + t42 * t50;
t72 = -t22 * t50 + t42 * t53;
t71 = cos(t77);
t69 = -pkin(2) * t83 - t25 * t50;
t67 = cos(t78) / 0.2e1;
t61 = t71 / 0.2e1 + t67;
t59 = cos(qJ(2));
t55 = -pkin(8) - qJ(4);
t54 = cos(pkin(6));
t38 = sin(t46);
t36 = pkin(4) * t52 + pkin(3);
t26 = t34 - t39 / 0.2e1;
t24 = t67 - t71 / 0.2e1;
t23 = t70 / 0.2e1 + t66;
t14 = t53 * t38 + t50 * t61;
t11 = t38 * t50 - t53 * t61;
t10 = t24 * t41 + t37 * t54;
t4 = t37 * t86 + t41 * t72;
t2 = -t37 * t84 + t41 * t73;
t1 = [(-m(2) - m(3) - t81) * g(3) (-t96 * mrSges(3,1) - (t35 - t44 / 0.2e1) * mrSges(3,2) + t95 * t82 + t94 * (t23 * t36 - t24 * t55 + t82) + t89 * t24 - t90 * t23) * g(3) + (-(t27 * t53 - t85) * mrSges(3,1) - (t53 * t26 - t50 * t59) * mrSges(3,2) + t95 * t76 + t94 * (-t11 * t36 - t55 * t73 + t76) + t89 * t73 + t90 * t11) * g(2) + (-(-t27 * t50 - t83) * mrSges(3,1) - (-t50 * t26 - t53 * t59) * mrSges(3,2) + t95 * t69 + t94 * (-t14 * t36 - t55 * t72 + t69) + t89 * t72 + t90 * t14) * g(1) (-g(3) * t54 + (-g(1) * t50 + g(2) * t53) * t51) * t81, t91 * (-g(1) * t14 - g(2) * t11 + g(3) * t23) (-t97 * (-t24 * t37 + t41 * t54) + t92 * t10) * g(3) + (t92 * t2 - t97 * (-t37 * t73 - t41 * t84)) * g(2) + (t92 * t4 - t97 * (-t37 * t72 + t41 * t86)) * g(1), -g(1) * ((t14 * t58 - t4 * t56) * mrSges(7,1) + (-t14 * t56 - t4 * t58) * mrSges(7,2)) - g(2) * ((t11 * t58 - t2 * t56) * mrSges(7,1) + (-t11 * t56 - t2 * t58) * mrSges(7,2)) - g(3) * ((-t10 * t56 - t23 * t58) * mrSges(7,1) + (-t10 * t58 + t23 * t56) * mrSges(7,2))];
taug  = t1(:);
