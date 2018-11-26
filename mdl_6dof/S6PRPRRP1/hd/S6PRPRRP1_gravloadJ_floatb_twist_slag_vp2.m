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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:59:30
% EndTime: 2018-11-23 14:59:31
% DurationCPUTime: 0.96s
% Computational Cost: add. (1240->110), mult. (1058->151), div. (0->0), fcn. (989->20), ass. (0->74)
t63 = cos(qJ(5));
t122 = -m(6) * pkin(4) - m(7) * (pkin(5) * t63 + pkin(4)) - mrSges(5,1);
t112 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t121 = mrSges(6,1) + mrSges(7,1);
t115 = -mrSges(6,2) - mrSges(7,2);
t55 = pkin(6) - qJ(2);
t43 = sin(t55) / 0.2e1;
t54 = pkin(6) + qJ(2);
t48 = sin(t54);
t120 = t48 / 0.2e1 + t43;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t119 = t112 * t61 + t122 * t64 - mrSges(4,1);
t109 = m(7) * pkin(5);
t116 = mrSges(4,2) - mrSges(5,3);
t44 = cos(t54) / 0.2e1;
t52 = cos(t55);
t35 = t52 / 0.2e1 + t44;
t114 = -m(5) - m(6) - m(7);
t60 = sin(qJ(5));
t113 = t115 * t60 + t121 * t63 - t122;
t111 = -t109 - t121;
t91 = m(4) - t114;
t53 = qJ(2) + pkin(11);
t50 = cos(t53);
t56 = sin(pkin(10));
t103 = t56 * t50;
t58 = cos(pkin(10));
t88 = pkin(6) - t53;
t80 = sin(t88);
t76 = t80 / 0.2e1;
t87 = pkin(6) + t53;
t79 = sin(t87);
t73 = -t79 / 0.2e1 + t76;
t21 = -t58 * t73 + t103;
t106 = t21 * t60;
t99 = t58 * t50;
t24 = t56 * t73 + t99;
t105 = t24 * t60;
t77 = cos(t88) / 0.2e1;
t81 = cos(t87);
t32 = t77 - t81 / 0.2e1;
t104 = t32 * t60;
t62 = sin(qJ(2));
t102 = t56 * t62;
t57 = sin(pkin(6));
t101 = t57 * t61;
t100 = t57 * t64;
t98 = t58 * t62;
t97 = t60 * t64;
t94 = t63 * t64;
t93 = t120 * pkin(2);
t92 = cos(pkin(6));
t33 = t35 * pkin(2);
t86 = -pkin(2) * t102 + t58 * t33;
t78 = -pkin(2) * t98 - t33 * t56;
t75 = t79 / 0.2e1;
t72 = t75 - t80 / 0.2e1;
t66 = t81 / 0.2e1 + t77;
t65 = cos(qJ(2));
t47 = sin(t53);
t34 = t43 - t48 / 0.2e1;
t31 = t76 + t75;
t26 = t32 * t64 + t61 * t92;
t25 = t32 * t61 - t64 * t92;
t23 = -t56 * t72 + t99;
t22 = t58 * t47 + t56 * t66;
t20 = t58 * t72 + t103;
t19 = t47 * t56 - t58 * t66;
t14 = t101 * t56 + t23 * t64;
t13 = -t100 * t56 + t23 * t61;
t12 = -t101 * t58 + t20 * t64;
t11 = t100 * t58 + t20 * t61;
t1 = [(-m(2) - m(3) - t91) * g(3) (-t120 * mrSges(3,1) - (t44 - t52 / 0.2e1) * mrSges(3,2) - m(4) * t93 - t104 * t109 + t114 * (t31 * pkin(3) + pkin(8) * t32 + t93) + t116 * t32 - t121 * (t31 * t94 + t104) + t115 * (-t31 * t97 + t32 * t63) + t119 * t31) * g(3) + (-(t35 * t58 - t102) * mrSges(3,1) - (t58 * t34 - t56 * t65) * mrSges(3,2) - m(4) * t86 - t106 * t109 + t114 * (-t19 * pkin(3) + pkin(8) * t21 + t86) - t121 * (-t19 * t94 + t106) + t115 * (t19 * t97 + t21 * t63) + t116 * t21 - t119 * t19) * g(2) + (-(-t35 * t56 - t98) * mrSges(3,1) - (-t56 * t34 - t58 * t65) * mrSges(3,2) - m(4) * t78 - t105 * t109 - t121 * (-t22 * t94 + t105) + t115 * (t22 * t97 + t24 * t63) + t114 * (-t22 * pkin(3) + pkin(8) * t24 + t78) + t116 * t24 - t119 * t22) * g(1) ((-g(1) * t56 + g(2) * t58) * t57 - g(3) * t92) * t91 (t112 * t26 + t113 * t25) * g(3) + (t113 * t11 + t112 * t12) * g(2) + (t112 * t14 + t113 * t13) * g(1) (t115 * (-t26 * t63 + t31 * t60) + t111 * (-t26 * t60 - t31 * t63)) * g(3) + (t115 * (-t12 * t63 - t19 * t60) + t111 * (-t12 * t60 + t19 * t63)) * g(2) + (t115 * (-t14 * t63 - t22 * t60) + t111 * (-t14 * t60 + t22 * t63)) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t25) * m(7)];
taug  = t1(:);
