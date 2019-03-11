% Calculate Gravitation load on the joints for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:16
% EndTime: 2019-03-08 19:55:18
% DurationCPUTime: 0.98s
% Computational Cost: add. (598->98), mult. (1518->149), div. (0->0), fcn. (1871->12), ass. (0->59)
t52 = cos(qJ(5));
t117 = m(6) * pkin(4) + m(7) * (pkin(5) * t52 + pkin(4)) + mrSges(5,1);
t116 = m(6) * pkin(9) - m(7) * (-qJ(6) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t115 = mrSges(6,1) + mrSges(7,1);
t107 = -mrSges(6,2) - mrSges(7,2);
t44 = sin(pkin(11));
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t84 = cos(pkin(11));
t33 = -t54 * t44 - t51 * t84;
t45 = sin(pkin(10));
t46 = cos(pkin(10));
t47 = cos(pkin(6));
t64 = -t51 * t44 + t54 * t84;
t59 = t47 * t64;
t23 = t33 * t46 - t45 * t59;
t85 = t33 * t47;
t24 = t45 * t85 + t46 * t64;
t91 = t47 * t54;
t65 = -t45 * t91 - t46 * t51;
t62 = t65 * pkin(2);
t114 = t23 * pkin(3) + pkin(8) * t24 + t62;
t113 = -t45 * t51 + t46 * t91;
t50 = sin(qJ(4));
t53 = cos(qJ(4));
t112 = t116 * t50 + t117 * t53 + mrSges(4,1);
t100 = m(7) * pkin(5);
t111 = -m(5) - m(6);
t108 = mrSges(4,2) - mrSges(5,3);
t106 = -m(7) + t111;
t49 = sin(qJ(5));
t105 = t107 * t49 + t115 * t52 + t117;
t103 = -t100 - t115;
t102 = m(4) - t106;
t19 = -t45 * t64 + t46 * t85;
t99 = t19 * t49;
t98 = t24 * t49;
t83 = sin(pkin(6));
t66 = t84 * t83;
t72 = t54 * t83;
t31 = t44 * t72 + t51 * t66;
t97 = t31 * t49;
t92 = t47 * t51;
t90 = t49 * t53;
t87 = t52 * t53;
t76 = t50 * t83;
t74 = t51 * t83;
t73 = t53 * t83;
t71 = t113 * pkin(2);
t40 = pkin(2) * t72;
t30 = t44 * t74 - t54 * t66;
t26 = t31 * t53 + t47 * t50;
t25 = t31 * t50 - t47 * t53;
t20 = t45 * t33 + t46 * t59;
t14 = t24 * t53 + t45 * t76;
t13 = t24 * t50 - t45 * t73;
t12 = -t19 * t53 - t46 * t76;
t11 = -t19 * t50 + t46 * t73;
t1 = [(-m(2) - m(3) - t102) * g(3) (-t97 * t100 - m(4) * t40 - mrSges(3,1) * t72 + mrSges(3,2) * t74 + t106 * (-t30 * pkin(3) + pkin(8) * t31 + t40) + t108 * t31 - t115 * (-t30 * t87 + t97) + t107 * (t30 * t90 + t31 * t52) + t112 * t30) * g(3) + (-t113 * mrSges(3,1) - (-t45 * t54 - t46 * t92) * mrSges(3,2) - m(4) * t71 + t99 * t100 + t106 * (t20 * pkin(3) - pkin(8) * t19 + t71) - t115 * (t20 * t87 - t99) + t107 * (-t19 * t52 - t20 * t90) - t108 * t19 - t112 * t20) * g(2) + (-t65 * mrSges(3,1) - (t45 * t92 - t46 * t54) * mrSges(3,2) - m(4) * t62 - m(7) * (pkin(5) * t98 + t114) - t115 * (t23 * t87 + t98) + t107 * (-t23 * t90 + t24 * t52) + t111 * t114 + t108 * t24 - t112 * t23) * g(1) (-g(3) * t47 + (-g(1) * t45 + g(2) * t46) * t83) * t102 (t105 * t25 - t116 * t26) * g(3) + (t105 * t11 - t116 * t12) * g(2) + (t105 * t13 - t116 * t14) * g(1) (t107 * (-t26 * t52 - t30 * t49) + t103 * (-t26 * t49 + t30 * t52)) * g(3) + (t107 * (-t12 * t52 + t20 * t49) + t103 * (-t12 * t49 - t20 * t52)) * g(2) + (t107 * (-t14 * t52 + t23 * t49) + t103 * (-t14 * t49 - t23 * t52)) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t25) * m(7)];
taug  = t1(:);
