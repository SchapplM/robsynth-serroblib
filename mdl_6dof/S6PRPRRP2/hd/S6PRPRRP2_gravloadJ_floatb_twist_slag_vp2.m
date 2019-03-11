% Calculate Gravitation load on the joints for
% S6PRPRRP2
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:59
% EndTime: 2019-03-08 20:00:01
% DurationCPUTime: 0.81s
% Computational Cost: add. (675->84), mult. (1738->131), div. (0->0), fcn. (2171->12), ass. (0->50)
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t77 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t80 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t102 = t77 * t60 - t80 * t63 - mrSges(5,1);
t101 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t104 = -m(6) - m(7);
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t105 = t104 * (pkin(4) * t64 + pkin(9) * t61) + t101 * t61 - mrSges(4,1) + t102 * t64;
t110 = -m(5) + t104;
t56 = sin(pkin(10));
t58 = cos(pkin(10));
t62 = sin(qJ(2));
t59 = cos(pkin(6));
t65 = cos(qJ(2));
t93 = t59 * t65;
t109 = -t56 * t62 + t58 * t93;
t106 = t80 * t60 + t77 * t63 - mrSges(4,2) + mrSges(5,3);
t55 = sin(pkin(11));
t87 = cos(pkin(11));
t72 = -t62 * t55 + t65 * t87;
t57 = sin(pkin(6));
t96 = t57 * t61;
t95 = t57 * t64;
t94 = t59 * t62;
t53 = t57 * t65 * pkin(2);
t85 = m(4) - t110;
t79 = t109 * pkin(2);
t49 = -t65 * t55 - t62 * t87;
t47 = t49 * t59;
t31 = -t47 * t58 + t56 * t72;
t32 = -t56 * t47 - t58 * t72;
t74 = -t56 * t93 - t58 * t62;
t71 = t74 * pkin(2);
t46 = t72 * t59;
t45 = t49 * t57;
t44 = t72 * t57;
t37 = -t45 * t64 + t59 * t61;
t36 = t45 * t61 + t59 * t64;
t33 = -t46 * t56 + t49 * t58;
t30 = t46 * t58 + t49 * t56;
t16 = -t32 * t64 + t56 * t96;
t15 = t32 * t61 + t56 * t95;
t14 = t31 * t64 - t58 * t96;
t13 = -t31 * t61 - t58 * t95;
t9 = t37 * t60 + t44 * t63;
t3 = t16 * t60 + t33 * t63;
t1 = t14 * t60 + t30 * t63;
t2 = [(-m(2) - m(3) - t85) * g(3) (-(mrSges(3,1) * t65 - mrSges(3,2) * t62) * t57 - m(4) * t53 + t110 * (t44 * pkin(3) - pkin(8) * t45 + t53) + t106 * t45 + t105 * t44) * g(3) + (-t109 * mrSges(3,1) - (-t56 * t65 - t58 * t94) * mrSges(3,2) - m(4) * t79 + t110 * (t30 * pkin(3) + pkin(8) * t31 + t79) - t106 * t31 + t105 * t30) * g(2) + (-t74 * mrSges(3,1) - (t56 * t94 - t58 * t65) * mrSges(3,2) - m(4) * t71 + t110 * (t33 * pkin(3) - t32 * pkin(8) + t71) + t106 * t32 + t105 * t33) * g(1) (-g(3) * t59 + (-g(1) * t56 + g(2) * t58) * t57) * t85 (t104 * (t36 * pkin(4) + pkin(9) * t37) + t101 * t37 + t102 * t36) * g(3) + (t104 * (t13 * pkin(4) + pkin(9) * t14) + t101 * t14 + t102 * t13) * g(2) + (t104 * (t15 * pkin(4) + pkin(9) * t16) + t101 * t16 + t102 * t15) * g(1) (t80 * t9 + t77 * (t37 * t63 - t44 * t60)) * g(3) + (t77 * (t14 * t63 - t30 * t60) + t80 * t1) * g(2) + (t77 * (t16 * t63 - t33 * t60) + t80 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
