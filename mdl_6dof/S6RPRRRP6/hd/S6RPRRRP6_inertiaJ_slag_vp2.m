% Calculate joint inertia matrix for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:27:01
% EndTime: 2018-11-23 16:27:02
% DurationCPUTime: 1.08s
% Computational Cost: add. (1967->260), mult. (3789->362), div. (0->0), fcn. (4188->8), ass. (0->102)
t148 = Ifges(6,5) + Ifges(7,5);
t147 = Ifges(6,6) + Ifges(7,6);
t132 = Ifges(6,3) + Ifges(7,3);
t106 = cos(qJ(4));
t100 = sin(pkin(10));
t101 = cos(pkin(10));
t104 = sin(qJ(3));
t135 = cos(qJ(3));
t75 = t135 * t100 + t104 * t101;
t119 = t106 * t75;
t74 = t100 * t104 - t135 * t101;
t146 = Ifges(5,5) * t119 + Ifges(5,3) * t74;
t102 = sin(qJ(5));
t103 = sin(qJ(4));
t105 = cos(qJ(5));
t76 = -t102 * t103 + t105 * t106;
t77 = t102 * t106 + t103 * t105;
t145 = t147 * t76 + t148 * t77;
t144 = (t105 * mrSges(6,1) + (-mrSges(6,2) - mrSges(7,2)) * t102) * pkin(4);
t131 = pkin(7) + qJ(2);
t81 = t131 * t100;
t82 = t131 * t101;
t51 = t104 * t82 + t135 * t81;
t143 = t51 ^ 2;
t142 = 2 * mrSges(7,1);
t141 = 0.2e1 * t51;
t89 = -pkin(2) * t101 - pkin(1);
t140 = 0.2e1 * t89;
t97 = t101 ^ 2;
t139 = m(7) * pkin(5);
t137 = -pkin(9) - pkin(8);
t136 = pkin(5) * t76;
t121 = t103 * t75;
t41 = pkin(3) * t74 - pkin(8) * t75 + t89;
t53 = -t104 * t81 + t135 * t82;
t20 = t103 * t41 + t106 * t53;
t15 = -pkin(9) * t121 + t20;
t19 = -t103 * t53 + t106 * t41;
t9 = pkin(4) * t74 - pkin(9) * t119 + t19;
t6 = t102 * t9 + t105 * t15;
t134 = m(7) * t102;
t133 = pkin(4) * t105;
t36 = t77 * t75;
t21 = -mrSges(7,2) * t74 - mrSges(7,3) * t36;
t22 = -mrSges(6,2) * t74 - mrSges(6,3) * t36;
t130 = t21 + t22;
t45 = -t76 * mrSges(7,1) + t77 * mrSges(7,2);
t86 = t137 * t103;
t87 = t137 * t106;
t57 = t102 * t86 - t105 * t87;
t127 = Ifges(5,5) * t103 + Ifges(5,6) * t106;
t126 = t100 ^ 2 + t97;
t125 = t103 ^ 2 + t106 ^ 2;
t124 = Ifges(5,4) * t103;
t123 = Ifges(5,4) * t106;
t122 = t102 * t76;
t37 = t76 * t75;
t5 = -t102 * t15 + t105 * t9;
t2 = pkin(5) * t74 - qJ(6) * t37 + t5;
t23 = mrSges(7,1) * t74 - mrSges(7,3) * t37;
t118 = m(7) * t2 + t23;
t91 = -pkin(4) * t106 - pkin(3);
t17 = t36 * mrSges(7,1) + t37 * mrSges(7,2);
t46 = -t76 * mrSges(6,1) + t77 * mrSges(6,2);
t116 = -t101 * mrSges(3,1) + t100 * mrSges(3,2);
t56 = t102 * t87 + t105 * t86;
t27 = pkin(4) * t121 + t51;
t34 = -qJ(6) * t77 + t56;
t115 = m(7) * t34 - t77 * mrSges(7,3);
t83 = -t106 * mrSges(5,1) + t103 * mrSges(5,2);
t114 = mrSges(5,1) * t103 + mrSges(5,2) * t106;
t113 = -t46 - t45;
t112 = t132 * t74 - t147 * t36 + t148 * t37;
t35 = qJ(6) * t76 + t57;
t111 = t56 * mrSges(6,1) + t34 * mrSges(7,1) - t57 * mrSges(6,2) - t35 * mrSges(7,2) + t145;
t3 = -qJ(6) * t36 + t6;
t110 = t5 * mrSges(6,1) + t2 * mrSges(7,1) - t6 * mrSges(6,2) - t3 * mrSges(7,2) + t112;
t109 = pkin(4) ^ 2;
t95 = t102 ^ 2 * t109;
t90 = pkin(5) + t133;
t85 = Ifges(5,1) * t103 + t123;
t84 = Ifges(5,2) * t106 + t124;
t62 = t75 * mrSges(4,2);
t61 = t102 * pkin(4) * t77;
t58 = t91 - t136;
t50 = Ifges(6,1) * t77 + Ifges(6,4) * t76;
t49 = Ifges(7,1) * t77 + Ifges(7,4) * t76;
t48 = Ifges(6,4) * t77 + Ifges(6,2) * t76;
t47 = Ifges(7,4) * t77 + Ifges(7,2) * t76;
t43 = mrSges(5,1) * t74 - mrSges(5,3) * t119;
t42 = -mrSges(5,2) * t74 - mrSges(5,3) * t121;
t40 = t114 * t75;
t26 = Ifges(5,5) * t74 + (Ifges(5,1) * t106 - t124) * t75;
t25 = Ifges(5,6) * t74 + (-Ifges(5,2) * t103 + t123) * t75;
t24 = mrSges(6,1) * t74 - mrSges(6,3) * t37;
t18 = mrSges(6,1) * t36 + mrSges(6,2) * t37;
t16 = pkin(5) * t36 + t27;
t13 = Ifges(6,1) * t37 - Ifges(6,4) * t36 + Ifges(6,5) * t74;
t12 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t74;
t11 = Ifges(6,4) * t37 - Ifges(6,2) * t36 + Ifges(6,6) * t74;
t10 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t74;
t1 = [-0.2e1 * pkin(1) * t116 + Ifges(3,2) * t97 + t62 * t140 + 0.2e1 * t20 * t42 + 0.2e1 * t19 * t43 + t40 * t141 + 0.2e1 * t3 * t21 + 0.2e1 * t6 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t5 * t24 + 0.2e1 * t27 * t18 + 0.2e1 * t16 * t17 + Ifges(2,3) + (Ifges(3,1) * t100 + 0.2e1 * Ifges(3,4) * t101) * t100 + (t12 + t13) * t37 - (t10 + t11) * t36 + 0.2e1 * t126 * qJ(2) * mrSges(3,3) + (mrSges(4,3) * t141 + Ifges(4,1) * t75 - t103 * t25 + t106 * t26) * t75 + (mrSges(4,1) * t140 - 0.2e1 * mrSges(4,3) * t53 + Ifges(4,2) * t74 + (-Ifges(5,6) * t103 - (2 * Ifges(4,4))) * t75 + t112 + t146) * t74 + m(3) * (t126 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t53 ^ 2 + t89 ^ 2 + t143) + m(5) * (t19 ^ 2 + t20 ^ 2 + t143) + m(7) * (t16 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t27 ^ 2 + t5 ^ 2 + t6 ^ 2); -m(3) * pkin(1) + t74 * mrSges(4,1) + t103 * t42 + t106 * t43 + t62 + t130 * t77 + (t23 + t24) * t76 + m(7) * (t2 * t76 + t3 * t77) + m(6) * (t5 * t76 + t6 * t77) + m(5) * (t103 * t20 + t106 * t19) + m(4) * t89 + t116; m(3) + m(4) + m(5) * t125 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t76 ^ 2 + t77 ^ 2); (t75 * t85 / 0.2e1 + pkin(8) * t42 + t20 * mrSges(5,3) + t25 / 0.2e1) * t106 + (-t19 * mrSges(5,3) - t75 * t84 / 0.2e1 - pkin(8) * t43 + t26 / 0.2e1) * t103 + m(6) * (t27 * t91 + t5 * t56 + t57 * t6) + m(7) * (t16 * t58 + t2 * t34 + t3 * t35) + t91 * t18 - Ifges(4,6) * t74 + Ifges(4,5) * t75 + t58 * t17 + t16 * t45 + t27 * t46 + (-t2 * mrSges(7,3) - t5 * mrSges(6,3) + t12 / 0.2e1 + t13 / 0.2e1) * t77 - t53 * mrSges(4,2) + t56 * t24 + t57 * t22 + t34 * t23 + t35 * t21 - pkin(3) * t40 + (t3 * mrSges(7,3) + t6 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1) * t76 + m(5) * (-pkin(3) * t51 + (-t19 * t103 + t20 * t106) * pkin(8)) + (t83 - mrSges(4,1)) * t51 - (t47 / 0.2e1 + t48 / 0.2e1) * t36 + (t49 / 0.2e1 + t50 / 0.2e1) * t37 + (t127 + t145) * t74 / 0.2e1; m(6) * (t56 * t76 + t57 * t77) + m(7) * (t34 * t76 + t35 * t77); -0.2e1 * pkin(3) * t83 + t103 * t85 + t106 * t84 + 0.2e1 * t58 * t45 + 0.2e1 * t91 * t46 + Ifges(4,3) + 0.2e1 * t125 * pkin(8) * mrSges(5,3) + m(7) * (t34 ^ 2 + t35 ^ 2 + t58 ^ 2) + m(6) * (t56 ^ 2 + t57 ^ 2 + t91 ^ 2) + m(5) * (t125 * pkin(8) ^ 2 + pkin(3) ^ 2) + (-0.2e1 * mrSges(6,3) * t56 - 0.2e1 * mrSges(7,3) * t34 + t49 + t50) * t77 + (0.2e1 * mrSges(6,3) * t57 + 0.2e1 * mrSges(7,3) * t35 + t47 + t48) * t76; t110 - Ifges(5,6) * t121 + t19 * mrSges(5,1) - t20 * mrSges(5,2) + t118 * t90 + (t105 * t24 + t130 * t102 + t3 * t134 + m(6) * (t102 * t6 + t105 * t5)) * pkin(4) + t146; m(6) * (t76 * t133 + t61) + m(7) * (t76 * t90 + t61) + t113 - t83; t115 * t90 - t114 * pkin(8) + (mrSges(7,3) * t122 + (-t105 * t77 + t122) * mrSges(6,3) + t35 * t134 + m(6) * (t102 * t57 + t105 * t56)) * pkin(4) + t111 + t127; t90 * t142 + Ifges(5,3) + m(7) * (t90 ^ 2 + t95) + m(6) * (t105 ^ 2 * t109 + t95) + 0.2e1 * t144 + t132; t118 * pkin(5) + t110; m(7) * t136 + t113; t115 * pkin(5) + t111; t90 * t139 + (pkin(5) + t90) * mrSges(7,1) + t144 + t132; (t142 + t139) * pkin(5) + t132; m(7) * t16 + t17; 0; m(7) * t58 + t45; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
