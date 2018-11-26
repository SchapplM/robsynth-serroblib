% Calculate Gravitation load on the joints for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:03:53
% EndTime: 2018-11-23 15:03:54
% DurationCPUTime: 1.06s
% Computational Cost: add. (1351->110), mult. (1111->149), div. (0->0), fcn. (1048->22), ass. (0->66)
t50 = qJ(5) + qJ(6);
t45 = sin(t50);
t46 = cos(t50);
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t116 = mrSges(5,1) + m(7) * (pkin(5) * t58 + pkin(4)) + t46 * mrSges(7,1) - t45 * mrSges(7,2) + m(6) * pkin(4) + t58 * mrSges(6,1) - t55 * mrSges(6,2);
t106 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) - mrSges(6,3);
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t105 = -t106 * t56 + t116 * t59 + mrSges(4,1);
t108 = -m(5) - m(6) - m(7);
t109 = -m(7) * pkin(5) - mrSges(6,1);
t49 = pkin(6) - qJ(2);
t36 = sin(t49) / 0.2e1;
t48 = pkin(6) + qJ(2);
t40 = sin(t48);
t113 = t40 / 0.2e1 + t36;
t37 = cos(t48) / 0.2e1;
t44 = cos(t49);
t29 = t44 / 0.2e1 + t37;
t104 = -t45 * mrSges(7,1) - t58 * mrSges(6,2) - t46 * mrSges(7,2) + t108 * pkin(8) + t109 * t55 + mrSges(4,2) - mrSges(5,3);
t47 = qJ(2) + pkin(12);
t39 = sin(t47);
t51 = sin(pkin(11));
t53 = cos(pkin(11));
t86 = pkin(6) - t47;
t74 = cos(t86) / 0.2e1;
t85 = pkin(6) + t47;
t78 = cos(t85);
t64 = t78 / 0.2e1 + t74;
t13 = t39 * t51 - t53 * t64;
t76 = sin(t85);
t72 = t76 / 0.2e1;
t77 = sin(t86);
t67 = t72 - t77 / 0.2e1;
t42 = cos(t47);
t97 = t51 * t42;
t14 = t53 * t67 + t97;
t52 = sin(pkin(6));
t95 = t52 * t56;
t8 = t14 * t59 - t53 * t95;
t100 = (t13 * t46 - t45 * t8) * mrSges(7,1) + (-t13 * t45 - t46 * t8) * mrSges(7,2);
t93 = t53 * t42;
t17 = -t51 * t67 + t93;
t10 = t17 * t59 + t51 * t95;
t16 = t53 * t39 + t51 * t64;
t99 = (-t10 * t45 + t16 * t46) * mrSges(7,1) + (-t10 * t46 - t16 * t45) * mrSges(7,2);
t26 = t74 - t78 / 0.2e1;
t54 = cos(pkin(6));
t20 = t26 * t59 + t54 * t56;
t73 = t77 / 0.2e1;
t25 = t73 + t72;
t98 = (-t20 * t45 - t25 * t46) * mrSges(7,1) + (-t20 * t46 + t25 * t45) * mrSges(7,2);
t57 = sin(qJ(2));
t96 = t51 * t57;
t94 = t52 * t59;
t92 = t53 * t57;
t91 = t113 * pkin(2);
t89 = m(4) - t108;
t27 = t29 * pkin(2);
t84 = -pkin(2) * t96 + t53 * t27;
t75 = -pkin(2) * t92 - t51 * t27;
t68 = -t76 / 0.2e1 + t73;
t60 = cos(qJ(2));
t28 = t36 - t40 / 0.2e1;
t1 = [(-m(2) - m(3) - t89) * g(3) (-t113 * mrSges(3,1) - (t37 - t44 / 0.2e1) * mrSges(3,2) - m(4) * t91 + t108 * (t25 * pkin(3) + t91) + t104 * t26 - t105 * t25) * g(3) + (-(t29 * t53 - t96) * mrSges(3,1) - (t28 * t53 - t51 * t60) * mrSges(3,2) - m(4) * t84 + t108 * (-t13 * pkin(3) + t84) + t104 * (-t53 * t68 + t97) + t105 * t13) * g(2) + (-(-t29 * t51 - t92) * mrSges(3,1) - (-t28 * t51 - t53 * t60) * mrSges(3,2) - m(4) * t75 + t108 * (-t16 * pkin(3) + t75) + t104 * (t51 * t68 + t93) + t105 * t16) * g(1) (-g(3) * t54 + (-g(1) * t51 + g(2) * t53) * t52) * t89 (t106 * t20 - t116 * (-t26 * t56 + t54 * t59)) * g(3) + (t106 * t8 - t116 * (-t14 * t56 - t53 * t94)) * g(2) + (-t116 * (-t17 * t56 + t51 * t94) + t106 * t10) * g(1) (-(-t20 * t58 + t25 * t55) * mrSges(6,2) - t98 + t109 * (-t20 * t55 - t25 * t58)) * g(3) + (-(-t13 * t55 - t58 * t8) * mrSges(6,2) - t100 + t109 * (t13 * t58 - t55 * t8)) * g(2) + (-(-t10 * t58 - t16 * t55) * mrSges(6,2) - t99 + t109 * (-t10 * t55 + t16 * t58)) * g(1), -g(1) * t99 - g(2) * t100 - g(3) * t98];
taug  = t1(:);
