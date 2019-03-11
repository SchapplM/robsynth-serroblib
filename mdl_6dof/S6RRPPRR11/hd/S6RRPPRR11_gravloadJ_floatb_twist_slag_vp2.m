% Calculate Gravitation load on the joints for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:04
% EndTime: 2019-03-09 09:39:06
% DurationCPUTime: 0.93s
% Computational Cost: add. (550->110), mult. (1045->145), div. (0->0), fcn. (1180->12), ass. (0->55)
t70 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t49 = sin(qJ(6));
t52 = cos(qJ(6));
t96 = m(7) * pkin(5) + t52 * mrSges(7,1) - t49 * mrSges(7,2) + mrSges(6,1);
t44 = pkin(11) + qJ(5);
t41 = sin(t44);
t42 = cos(t44);
t45 = sin(pkin(11));
t47 = cos(pkin(11));
t63 = -t45 * mrSges(5,1) - t47 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t102 = -t96 * t41 - t70 * t42 + t63;
t98 = m(4) + m(5);
t97 = m(6) + m(7);
t64 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t89 = t49 * mrSges(7,1) + t52 * mrSges(7,2) - t64;
t94 = m(5) + t97;
t76 = m(4) + t94;
t92 = t76 * qJ(3) - t63;
t90 = (-t98 - t97) * qJ(3) + t102;
t88 = pkin(4) * t45;
t46 = sin(pkin(6));
t50 = sin(qJ(2));
t87 = t46 * t50;
t51 = sin(qJ(1));
t86 = t46 * t51;
t53 = cos(qJ(2));
t85 = t46 * t53;
t54 = cos(qJ(1));
t84 = t46 * t54;
t83 = pkin(2) * t85 + qJ(3) * t87;
t82 = t54 * pkin(1) + pkin(8) * t86;
t81 = cos(pkin(6));
t72 = t51 * t81;
t29 = -t50 * t72 + t53 * t54;
t78 = t29 * pkin(2) + t82;
t74 = -t51 * pkin(1) + pkin(8) * t84;
t71 = t54 * t81;
t27 = t50 * t71 + t51 * t53;
t69 = t27 * pkin(2) - t74;
t28 = t54 * t50 + t53 * t72;
t40 = pkin(4) * t47 + pkin(3);
t48 = -pkin(9) - qJ(4);
t65 = t28 * t88 - t29 * t48 + t40 * t86 + t78;
t26 = t50 * t51 - t53 * t71;
t7 = -t26 * t41 + t42 * t84;
t5 = t26 * t42 + t41 * t84;
t56 = mrSges(2,2) + (-m(5) * pkin(3) - mrSges(5,1) * t47 + mrSges(5,2) * t45 - mrSges(4,1) - mrSges(3,3)) * t46;
t24 = t28 * pkin(2);
t22 = t26 * pkin(2);
t13 = -t41 * t85 + t81 * t42;
t4 = t28 * t41 + t42 * t86;
t3 = -t28 * t42 + t41 * t86;
t2 = t29 * t49 + t4 * t52;
t1 = t29 * t52 - t4 * t49;
t6 = [(-t54 * mrSges(2,1) - m(3) * t82 - m(6) * t65 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t65) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t70 * t3 - t92 * t28 + t56 * t51 + t64 * t29 - t98 * t78) * g(2) + (t51 * mrSges(2,1) - m(3) * t74 + t70 * t5 - t96 * t7 + t92 * t26 + t56 * t54 + t89 * t27 + t98 * t69 + t97 * (t26 * t88 - t27 * t48 - t40 * t84 + t69)) * g(1) (-t98 * t83 - t97 * (t87 * t88 + t83) + ((t97 * t48 - t89) * t53 + t102 * t50) * t46) * g(3) + (-t97 * (t26 * t48 + t27 * t88 - t22) + t98 * t22 + t90 * t27 + t89 * t26) * g(2) + (-t97 * (t28 * t48 + t29 * t88 - t24) + t98 * t24 + t90 * t29 + t89 * t28) * g(1) (-g(1) * t28 - g(2) * t26 + g(3) * t85) * t76, t94 * (-g(1) * t29 - g(2) * t27 - g(3) * t87) (t70 * t13 - t96 * (-t81 * t41 - t42 * t85)) * g(3) + (-t96 * t5 - t70 * t7) * g(2) + (t96 * t3 + t70 * t4) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t27 * t52 + t49 * t7) * mrSges(7,1) + (-t27 * t49 + t52 * t7) * mrSges(7,2)) - g(3) * ((-t13 * t49 + t52 * t87) * mrSges(7,1) + (-t13 * t52 - t49 * t87) * mrSges(7,2))];
taug  = t6(:);
