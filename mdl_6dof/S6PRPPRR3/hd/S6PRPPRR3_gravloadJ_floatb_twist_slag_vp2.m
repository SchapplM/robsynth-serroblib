% Calculate Gravitation load on the joints for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 14:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:54:38
% EndTime: 2018-11-23 14:54:38
% DurationCPUTime: 0.60s
% Computational Cost: add. (1058->96), mult. (1222->132), div. (0->0), fcn. (1203->16), ass. (0->66)
t100 = -m(6) - m(7);
t52 = sin(qJ(6));
t55 = cos(qJ(6));
t98 = -m(7) * pkin(5) - mrSges(7,1) * t55 + mrSges(7,2) * t52 - mrSges(6,1);
t96 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t99 = mrSges(3,1) + mrSges(4,1);
t97 = mrSges(3,2) - mrSges(4,3);
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t95 = t96 * t53 + t98 * t56 - mrSges(5,1);
t94 = -t52 * mrSges(7,1) - t55 * mrSges(7,2) + t100 * pkin(8) + mrSges(5,2) - mrSges(6,3);
t47 = sin(pkin(10));
t50 = cos(pkin(10));
t54 = sin(qJ(2));
t83 = pkin(6) + qJ(2);
t70 = cos(t83) / 0.2e1;
t84 = pkin(6) - qJ(2);
t76 = cos(t84);
t59 = t76 / 0.2e1 + t70;
t28 = t47 * t54 - t50 * t59;
t46 = sin(pkin(11));
t93 = t28 * t46;
t49 = cos(pkin(11));
t92 = t28 * t49;
t31 = t47 * t59 + t50 * t54;
t91 = t31 * t46;
t90 = t31 * t49;
t57 = cos(qJ(2));
t89 = t47 * t57;
t48 = sin(pkin(6));
t88 = t48 * t53;
t87 = t48 * t56;
t86 = t50 * t57;
t85 = m(5) - t100;
t82 = m(4) + t85;
t74 = sin(t83);
t68 = t74 / 0.2e1;
t75 = sin(t84);
t69 = t75 / 0.2e1;
t40 = t68 + t69;
t41 = t70 - t76 / 0.2e1;
t80 = t40 * t49 - t41 * t46;
t18 = t40 * t46 + t41 * t49;
t63 = t69 - t74 / 0.2e1;
t30 = -t50 * t63 + t89;
t79 = -t28 * pkin(2) + t30 * qJ(3);
t33 = t47 * t63 + t86;
t78 = -t31 * pkin(2) + t33 * qJ(3);
t77 = t40 * pkin(2) - t41 * qJ(3);
t73 = -t28 * pkin(3) + t79;
t72 = -t31 * pkin(3) + t78;
t71 = t40 * pkin(3) + t77;
t62 = t68 - t75 / 0.2e1;
t51 = cos(pkin(6));
t32 = -t47 * t62 + t86;
t29 = t50 * t62 + t89;
t16 = -t18 * t56 - t51 * t53;
t14 = t32 * t49 + t91;
t13 = t32 * t46 - t90;
t12 = t33 * t46 - t90;
t10 = t29 * t49 + t93;
t9 = t29 * t46 - t92;
t8 = t30 * t46 - t92;
t4 = t14 * t56 - t47 * t88;
t2 = t10 * t56 + t50 * t88;
t1 = [(-m(2) - m(3) - t82) * g(3) (-m(4) * t77 - m(5) * t71 + t100 * (pkin(4) * t80 + t71) - t97 * t41 - t99 * t40 + t95 * t80 + t94 * t18) * g(3) + (-m(4) * t79 - m(5) * t73 + t100 * (t8 * pkin(4) + t73) + t97 * t30 + t99 * t28 + t95 * t8 + t94 * (-t30 * t49 - t93)) * g(2) + (-m(4) * t78 - m(5) * t72 + t100 * (t12 * pkin(4) + t72) + t97 * t33 + t99 * t31 + t95 * t12 + t94 * (-t33 * t49 - t91)) * g(1) (-g(1) * t31 - g(2) * t28 + g(3) * t40) * t82 (t51 * g(3) + (g(1) * t47 - g(2) * t50) * t48) * t85 (t96 * t16 + t98 * (t18 * t53 - t51 * t56)) * g(3) + (t96 * t2 + t98 * (-t10 * t53 + t50 * t87)) * g(2) + (t96 * t4 + t98 * (-t14 * t53 - t47 * t87)) * g(1), -g(1) * ((t13 * t55 - t4 * t52) * mrSges(7,1) + (-t13 * t52 - t4 * t55) * mrSges(7,2)) - g(2) * ((-t2 * t52 + t55 * t9) * mrSges(7,1) + (-t2 * t55 - t52 * t9) * mrSges(7,2)) - g(3) * ((-t16 * t52 + t55 * t80) * mrSges(7,1) + (-t16 * t55 - t52 * t80) * mrSges(7,2))];
taug  = t1(:);
