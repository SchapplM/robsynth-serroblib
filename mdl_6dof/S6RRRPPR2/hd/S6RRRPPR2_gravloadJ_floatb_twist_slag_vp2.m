% Calculate Gravitation load on the joints for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2018-11-23 17:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:33:28
% EndTime: 2018-11-23 17:33:29
% DurationCPUTime: 0.71s
% Computational Cost: add. (490->119), mult. (447->130), div. (0->0), fcn. (370->10), ass. (0->64)
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t117 = mrSges(7,1) * t41 + mrSges(7,2) * t44;
t40 = qJ(2) + qJ(3);
t35 = pkin(10) + t40;
t31 = sin(t35);
t32 = cos(t35);
t36 = sin(t40);
t37 = cos(t40);
t116 = mrSges(4,1) * t36 + mrSges(5,1) * t31 + mrSges(4,2) * t37 + mrSges(5,2) * t32;
t115 = t117 * t32;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t109 = g(1) * t46 + g(2) * t43;
t114 = -t37 * mrSges(4,1) + t36 * mrSges(4,2) + (-mrSges(5,1) + mrSges(6,2)) * t32 + (mrSges(5,2) - mrSges(6,3)) * t31;
t103 = m(6) + m(7);
t24 = t31 * qJ(5);
t87 = t32 * t46;
t113 = pkin(4) * t87 + t46 * t24;
t89 = t31 * t46;
t110 = -mrSges(6,2) * t89 - mrSges(6,3) * t87 - t115 * t46;
t90 = t31 * t43;
t108 = -mrSges(6,2) * t90 + (-mrSges(6,3) * t32 - t115) * t43;
t106 = -t32 * mrSges(7,3) - t117 * t31 + t114;
t45 = cos(qJ(2));
t38 = t45 * pkin(2);
t42 = sin(qJ(2));
t60 = t45 * mrSges(3,1) - t42 * mrSges(3,2);
t105 = -mrSges(2,1) - m(3) * pkin(1) - t60 - m(4) * (t38 + pkin(1)) + t114;
t47 = -pkin(8) - pkin(7);
t39 = -qJ(4) + t47;
t104 = -m(3) * pkin(7) + m(4) * t47 - m(7) * (pkin(5) - t39) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t101 = pkin(2) * t42;
t100 = pkin(3) * t36;
t33 = pkin(3) * t37;
t97 = g(3) * t32;
t29 = t32 * pkin(4);
t13 = -t100 - t101;
t77 = qJ(5) * t32;
t14 = t43 * t77;
t95 = t43 * t13 + t14;
t16 = t46 * t77;
t94 = t46 * t13 + t16;
t85 = t41 * t43;
t84 = t41 * t46;
t83 = t43 * t44;
t82 = t44 * t46;
t78 = t33 + t38;
t72 = t29 + t24 + t33;
t12 = pkin(1) + t78;
t5 = t46 * t12;
t69 = -t39 * t43 + t5;
t65 = -t12 - t24;
t64 = t38 + t72;
t63 = m(7) * (-pkin(4) - pkin(9)) - mrSges(7,3);
t62 = -pkin(4) * t31 - t100;
t54 = t63 * t31;
t48 = -m(7) * t100 + t54;
t28 = t32 * pkin(9);
t4 = -t31 * t85 + t82;
t3 = t31 * t83 + t84;
t2 = t31 * t84 + t83;
t1 = t31 * t82 - t85;
t6 = [(-m(5) * t69 - m(6) * (t69 + t113) - m(7) * (pkin(9) * t87 + t113 + t5) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - mrSges(7,3) * t87 + t105 * t46 + t104 * t43) * g(2) + (-t4 * mrSges(7,1) + t3 * mrSges(7,2) + ((m(5) + m(6)) * t39 + t104) * t46 + (m(5) * t12 - m(6) * (t65 - t29) - m(7) * t65 - t63 * t32 - t105) * t43) * g(1) (-m(6) * (-pkin(4) * t90 + t95) - m(7) * t95 - t43 * t54 + t108) * g(2) + (-m(6) * (-pkin(4) * t89 + t94) - m(7) * t94 - t46 * t54 + t110) * g(1) + (-t60 - m(4) * t38 - m(5) * t78 - m(6) * t64 - m(7) * (t28 + t64) + t106) * g(3) + t109 * (m(4) * t101 - m(5) * t13 + mrSges(3,1) * t42 + mrSges(3,2) * t45 + t116) (-m(6) * (t62 * t43 + t14) - m(7) * t14 - t48 * t43 + t108) * g(2) + (-m(6) * (t62 * t46 + t16) - m(7) * t16 - t48 * t46 + t110) * g(1) + (-m(5) * t33 - m(6) * t72 - m(7) * (t28 + t72) + t106) * g(3) + (m(5) * t100 + t116) * t109 (-g(1) * t43 + g(2) * t46) * (m(5) + t103) (-t109 * t31 + t97) * t103, -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t3 + mrSges(7,2) * t4) - (-mrSges(7,1) * t44 + mrSges(7,2) * t41) * t97];
taug  = t6(:);
