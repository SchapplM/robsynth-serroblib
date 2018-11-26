% Calculate Gravitation load on the joints for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:55:52
% EndTime: 2018-11-23 14:55:53
% DurationCPUTime: 0.90s
% Computational Cost: add. (1154->99), mult. (962->132), div. (0->0), fcn. (889->22), ass. (0->66)
t43 = pkin(12) + qJ(6);
t35 = sin(t43);
t39 = cos(t43);
t47 = sin(pkin(12));
t50 = cos(pkin(12));
t96 = mrSges(5,1) + m(7) * (pkin(5) * t50 + pkin(4)) + t39 * mrSges(7,1) - t35 * mrSges(7,2) + m(6) * pkin(4) + t50 * mrSges(6,1) - t47 * mrSges(6,2);
t95 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3);
t53 = sin(qJ(4));
t55 = cos(qJ(4));
t94 = -t53 * t95 + t55 * t96 + mrSges(4,1);
t101 = m(6) + m(7);
t97 = -m(5) - t101;
t46 = pkin(6) - qJ(2);
t31 = sin(t46) / 0.2e1;
t45 = pkin(6) + qJ(2);
t37 = sin(t45);
t102 = t37 / 0.2e1 + t31;
t32 = cos(t45) / 0.2e1;
t42 = cos(t46);
t23 = t42 / 0.2e1 + t32;
t79 = m(4) - t97;
t93 = -t35 * mrSges(7,1) - t50 * mrSges(6,2) - t39 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t47 + t97 * pkin(8);
t44 = qJ(2) + pkin(11);
t75 = pkin(6) + t44;
t30 = sin(t75);
t92 = t30 / 0.2e1;
t40 = cos(t44);
t48 = sin(pkin(10));
t89 = t48 * t40;
t54 = sin(qJ(2));
t88 = t48 * t54;
t49 = sin(pkin(6));
t87 = t49 * t53;
t86 = t49 * t55;
t51 = cos(pkin(10));
t85 = t51 * t40;
t84 = t51 * t54;
t83 = t102 * pkin(2);
t82 = cos(pkin(6));
t76 = pkin(6) - t44;
t70 = sin(t76);
t81 = t92 - t70 / 0.2e1;
t21 = t23 * pkin(2);
t74 = -pkin(2) * t88 + t21 * t51;
t71 = cos(t75);
t69 = -pkin(2) * t84 - t21 * t48;
t68 = cos(t76) / 0.2e1;
t67 = t70 / 0.2e1;
t61 = -t30 / 0.2e1 + t67;
t59 = t71 / 0.2e1 + t68;
t56 = cos(qJ(2));
t36 = sin(t44);
t22 = t31 - t37 / 0.2e1;
t20 = t68 - t71 / 0.2e1;
t19 = t67 + t92;
t14 = t20 * t55 + t53 * t82;
t13 = t20 * t53 - t55 * t82;
t11 = -t48 * t81 + t85;
t10 = t36 * t51 + t48 * t59;
t8 = t51 * t81 + t89;
t7 = t36 * t48 - t51 * t59;
t4 = t11 * t55 + t48 * t87;
t3 = t11 * t53 - t48 * t86;
t2 = -t51 * t87 + t55 * t8;
t1 = t51 * t86 + t53 * t8;
t5 = [(-m(2) - m(3) - t79) * g(3) (-t102 * mrSges(3,1) - (t32 - t42 / 0.2e1) * mrSges(3,2) - m(4) * t83 + t97 * (pkin(3) * t19 + t83) + t93 * t20 - t94 * t19) * g(3) + (-(t23 * t51 - t88) * mrSges(3,1) - (t22 * t51 - t48 * t56) * mrSges(3,2) - m(4) * t74 + t97 * (-pkin(3) * t7 + t74) + t93 * (-t51 * t61 + t89) + t94 * t7) * g(2) + (-(-t23 * t48 - t84) * mrSges(3,1) - (-t22 * t48 - t51 * t56) * mrSges(3,2) - m(4) * t69 + t97 * (-pkin(3) * t10 + t69) + t93 * (t48 * t61 + t85) + t94 * t10) * g(1) ((-g(1) * t48 + g(2) * t51) * t49 - g(3) * t82) * t79 (t13 * t96 + t14 * t95) * g(3) + (t1 * t96 + t2 * t95) * g(2) + (t3 * t96 + t4 * t95) * g(1), t101 * (-g(1) * t3 - g(2) * t1 - g(3) * t13) -g(1) * ((t10 * t39 - t35 * t4) * mrSges(7,1) + (-t10 * t35 - t39 * t4) * mrSges(7,2)) - g(2) * ((-t2 * t35 + t39 * t7) * mrSges(7,1) + (-t2 * t39 - t35 * t7) * mrSges(7,2)) - g(3) * ((-t14 * t35 - t19 * t39) * mrSges(7,1) + (-t14 * t39 + t19 * t35) * mrSges(7,2))];
taug  = t5(:);
