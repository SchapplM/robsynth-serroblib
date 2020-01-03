% Calculate Gravitation load on the joints for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:10
% EndTime: 2019-12-31 22:08:13
% DurationCPUTime: 1.00s
% Computational Cost: add. (460->109), mult. (1111->158), div. (0->0), fcn. (1305->10), ass. (0->57)
t48 = cos(qJ(4));
t40 = pkin(4) * t48 + pkin(3);
t105 = -m(6) * t40 - mrSges(4,1);
t104 = m(6) * (-qJ(5) - pkin(9)) + mrSges(4,2) - mrSges(6,3);
t103 = mrSges(5,1) + mrSges(6,1);
t97 = mrSges(5,2) + mrSges(6,2);
t45 = sin(qJ(3));
t49 = cos(qJ(3));
t102 = t104 * t45 + t105 * t49 - mrSges(3,1);
t46 = sin(qJ(2));
t47 = sin(qJ(1));
t50 = cos(qJ(2));
t72 = cos(pkin(5));
t90 = cos(qJ(1));
t59 = t72 * t90;
t28 = t46 * t59 + t47 * t50;
t42 = sin(pkin(5));
t67 = t42 * t90;
t14 = t28 * t49 - t45 * t67;
t27 = t46 * t47 - t50 * t59;
t44 = sin(qJ(4));
t101 = t14 * t44 - t27 * t48;
t100 = -t14 * t48 - t27 * t44;
t91 = m(6) * pkin(4);
t98 = mrSges(3,2) - mrSges(4,3);
t71 = m(4) + m(5) + m(6);
t96 = m(5) * pkin(3) + t103 * t48 - t97 * t44 - t105;
t95 = -t91 - t103;
t56 = -t44 * t91 + t98;
t53 = -m(5) * pkin(9) - mrSges(5,3) + t104;
t94 = m(5) * (pkin(3) * t49 + pkin(9) * t45) + t45 * mrSges(5,3);
t93 = pkin(8) * t71 - t56;
t92 = t94 - t102;
t87 = t28 * t44;
t62 = t47 * t72;
t30 = -t46 * t62 + t50 * t90;
t86 = t30 * t44;
t85 = t42 * t46;
t84 = t42 * t47;
t83 = t42 * t49;
t82 = t42 * t50;
t81 = t44 * t49;
t78 = t48 * t49;
t77 = t49 * t50;
t73 = t90 * pkin(1) + pkin(7) * t84;
t69 = t30 * pkin(2) + t73;
t65 = -pkin(1) * t47 + pkin(7) * t67;
t61 = -t28 * pkin(2) + t65;
t18 = t30 * t49 + t45 * t84;
t29 = t46 * t90 + t50 * t62;
t5 = -t18 * t44 + t29 * t48;
t13 = t28 * t45 + t49 * t67;
t26 = t45 * t72 + t46 * t83;
t25 = t45 * t85 - t49 * t72;
t17 = t30 * t45 - t47 * t83;
t6 = t18 * t48 + t29 * t44;
t1 = [(-t90 * mrSges(2,1) - m(3) * t73 - t30 * mrSges(3,1) - m(4) * t69 - t18 * mrSges(4,1) - m(5) * (pkin(3) * t18 + t69) - m(6) * (t18 * t40 + t69) - t103 * t6 - t97 * t5 + (-mrSges(3,3) * t42 + mrSges(2,2)) * t47 - t93 * t29 + t53 * t17) * g(2) + (t47 * mrSges(2,1) + t90 * mrSges(2,2) - m(3) * t65 + t28 * mrSges(3,1) - mrSges(3,3) * t67 - m(4) * t61 + t14 * mrSges(4,1) - m(5) * (-pkin(3) * t14 + t61) - m(6) * (-t14 * t40 + t61) - t103 * t100 - t97 * t101 + t93 * t27 - t53 * t13) * g(1), (-t87 * t91 - t103 * (-t27 * t78 + t87) - t97 * (t27 * t81 + t28 * t48) - t71 * (-t27 * pkin(2) + pkin(8) * t28) + t98 * t28 + t92 * t27) * g(2) + (-t86 * t91 - t97 * (t29 * t81 + t30 * t48) - t71 * (-t29 * pkin(2) + pkin(8) * t30) + t98 * t30 - t103 * (-t29 * t78 + t86) + t92 * t29) * g(1) + (-t94 * t82 - t71 * (pkin(2) * t82 + pkin(8) * t85) + (-t103 * (t44 * t46 + t48 * t77) - t97 * (-t44 * t77 + t46 * t48) + t102 * t50 + t56 * t46) * t42) * g(3), (t96 * t25 + t53 * t26) * g(3) + (t96 * t13 + t53 * t14) * g(2) + (t96 * t17 + t53 * t18) * g(1), (-t97 * (-t26 * t48 + t44 * t82) + t95 * (-t26 * t44 - t48 * t82)) * g(3) + (-t97 * t100 - t95 * t101) * g(2) + (t95 * t5 + t97 * t6) * g(1), (-g(1) * t17 - g(2) * t13 - g(3) * t25) * m(6)];
taug = t1(:);
