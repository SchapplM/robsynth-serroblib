% Calculate joint inertia matrix for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:52:35
% EndTime: 2018-11-23 16:52:36
% DurationCPUTime: 0.90s
% Computational Cost: add. (1541->224), mult. (2710->303), div. (0->0), fcn. (2884->8), ass. (0->82)
t64 = sin(qJ(2));
t66 = cos(qJ(2));
t110 = t64 ^ 2 + t66 ^ 2;
t62 = sin(qJ(6));
t65 = cos(qJ(6));
t83 = t62 ^ 2 + t65 ^ 2;
t79 = t83 * mrSges(7,3);
t44 = -t65 * mrSges(7,1) + mrSges(7,2) * t62;
t87 = mrSges(6,1) - t44;
t60 = sin(pkin(10));
t61 = cos(pkin(10));
t35 = -t60 * t64 - t61 * t66;
t36 = -t60 * t66 + t61 * t64;
t63 = sin(qJ(5));
t96 = cos(qJ(5));
t17 = -t35 * t96 + t36 * t63;
t18 = t63 * t35 + t36 * t96;
t93 = t18 * t62;
t10 = -mrSges(7,2) * t17 - mrSges(7,3) * t93;
t92 = t18 * t65;
t11 = mrSges(7,1) * t17 - mrSges(7,3) * t92;
t109 = t65 * t10 - t62 * t11;
t108 = -m(4) * pkin(2) - mrSges(4,1);
t107 = -m(7) * pkin(5) - t87;
t86 = pkin(7) - qJ(4);
t80 = t86 * t66;
t81 = t86 * t64;
t26 = t60 * t81 + t61 * t80;
t13 = pkin(8) * t35 + t26;
t25 = -t60 * t80 + t61 * t81;
t70 = -t36 * pkin(8) + t25;
t6 = t13 * t63 - t70 * t96;
t106 = t6 ^ 2;
t105 = 0.2e1 * t6;
t33 = t60 * t63 - t61 * t96;
t104 = t33 ^ 2;
t43 = -t66 * pkin(2) - t64 * qJ(3) - pkin(1);
t31 = t66 * pkin(3) - t43;
t24 = -pkin(4) * t35 + t31;
t103 = 0.2e1 * t24;
t102 = 0.2e1 * t35;
t101 = -0.2e1 * t43;
t100 = -0.2e1 * t44;
t99 = t65 / 0.2e1;
t98 = t33 * t6;
t8 = t13 * t96 + t63 * t70;
t97 = t8 * mrSges(6,2);
t95 = Ifges(7,4) * t62;
t94 = Ifges(7,4) * t65;
t67 = -pkin(2) - pkin(3);
t41 = -qJ(3) * t60 + t61 * t67;
t40 = -pkin(4) + t41;
t42 = qJ(3) * t61 + t60 * t67;
t22 = t40 * t96 - t63 * t42;
t91 = t22 * mrSges(6,1);
t23 = t63 * t40 + t96 * t42;
t90 = t23 * mrSges(6,2);
t85 = Ifges(7,5) * t92 + Ifges(7,3) * t17;
t84 = t110 * pkin(7) ^ 2;
t45 = Ifges(7,5) * t62 + Ifges(7,6) * t65;
t82 = t45 / 0.2e1 - Ifges(6,6);
t21 = -pkin(9) + t23;
t78 = t83 * t21;
t37 = t60 * t96 + t63 * t61;
t77 = t83 * t37;
t76 = -t35 * mrSges(5,1) + t36 * mrSges(5,2);
t5 = pkin(5) * t17 - pkin(9) * t18 + t24;
t1 = t5 * t65 - t62 * t8;
t2 = t5 * t62 + t65 * t8;
t74 = -t1 * t62 + t2 * t65;
t73 = mrSges(7,1) * t62 + mrSges(7,2) * t65;
t46 = Ifges(7,2) * t65 + t95;
t47 = Ifges(7,1) * t62 + t94;
t72 = t65 * t46 + t62 * t47 + Ifges(6,3);
t71 = Ifges(6,5) + t47 * t99 - t62 * t46 / 0.2e1;
t32 = t37 ^ 2;
t20 = pkin(5) - t22;
t15 = t18 * mrSges(6,2);
t9 = t73 * t18;
t4 = Ifges(7,5) * t17 + (Ifges(7,1) * t65 - t95) * t18;
t3 = Ifges(7,6) * t17 + (-Ifges(7,2) * t62 + t94) * t18;
t7 = [t26 * mrSges(5,3) * t102 + 0.2e1 * t31 * t76 + Ifges(5,2) * t35 ^ 2 + t15 * t103 + t9 * t105 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t101 + (Ifges(3,2) + Ifges(4,3)) * t66) * t66 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t101 + (Ifges(3,1) + Ifges(4,1)) * t64 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t66) * t64 + (mrSges(6,1) * t103 - 0.2e1 * t8 * mrSges(6,3) + Ifges(6,2) * t17 + t85) * t17 + (-0.2e1 * mrSges(5,3) * t25 + Ifges(5,1) * t36 + Ifges(5,4) * t102) * t36 + (mrSges(6,3) * t105 + Ifges(6,1) * t18 - t3 * t62 + t4 * t65 + (-Ifges(7,6) * t62 - (2 * Ifges(6,4))) * t17) * t18 + m(4) * (t43 ^ 2 + t84) + m(3) * (pkin(1) ^ 2 + t84) + m(6) * (t24 ^ 2 + t8 ^ 2 + t106) + m(5) * (t25 ^ 2 + t26 ^ 2 + t31 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t106) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(7) * t110; -t25 * mrSges(5,1) + t26 * mrSges(5,2) + t97 - Ifges(5,5) * t36 - Ifges(5,6) * t35 + t20 * t9 + t87 * t6 + (t42 * t35 - t41 * t36) * mrSges(5,3) + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t66 + (-t3 / 0.2e1 + t21 * t10 - t2 * mrSges(7,3)) * t65 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t64 + (-t4 / 0.2e1 + t1 * mrSges(7,3) - t21 * t11) * t62 + (-t23 * mrSges(6,3) - t82) * t17 + m(7) * (t20 * t6 + t21 * t74) + m(5) * (t25 * t41 + t26 * t42) + m(6) * (-t22 * t6 + t23 * t8) + (-t22 * mrSges(6,3) - t71) * t18 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t66 + (-mrSges(3,1) + t108) * t64) * pkin(7); 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t41 * mrSges(5,1) - 0.2e1 * t91 + 0.2e1 * t42 * mrSges(5,2) + 0.2e1 * t90 + 0.2e1 * qJ(3) * mrSges(4,3) + t20 * t100 + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) - 0.2e1 * t21 * t79 + m(7) * (t21 ^ 2 * t83 + t20 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t41 ^ 2 + t42 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t72; (m(4) * pkin(7) + mrSges(4,2)) * t64 + (t18 * mrSges(6,3) + t9) * t33 + (t35 * t60 - t36 * t61) * mrSges(5,3) + (-t17 * mrSges(6,3) + t109) * t37 + m(7) * (t37 * t74 + t98) + m(6) * (t37 * t8 + t98) + m(5) * (t25 * t61 + t26 * t60); -t61 * mrSges(5,1) + t60 * mrSges(5,2) + t87 * t33 + (mrSges(6,2) - t79) * t37 + m(7) * (t20 * t33 + t21 * t77) + m(6) * (-t22 * t33 + t23 * t37) + m(5) * (t41 * t61 + t42 * t60) + t108; m(4) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(6) * (t32 + t104) + m(7) * (t32 * t83 + t104); t17 * mrSges(6,1) + t62 * t10 + t65 * t11 + t15 + m(7) * (t1 * t65 + t2 * t62) + m(6) * t24 + m(5) * t31 + t76; 0; 0; m(7) * t83 + m(5) + m(6); t62 * t4 / 0.2e1 + t3 * t99 - pkin(5) * t9 - t97 + t82 * t17 + t74 * mrSges(7,3) + t71 * t18 + t107 * t6 + (m(7) * t74 + t109) * pkin(9); m(7) * (-pkin(5) * t20 + pkin(9) * t78) - t90 + t91 + (t20 + pkin(5)) * t44 + (-pkin(9) * t83 + t78) * mrSges(7,3) - t72; -t37 * mrSges(6,2) + (m(7) * pkin(9) + mrSges(7,3)) * t77 + t107 * t33; 0; pkin(5) * t100 + m(7) * (pkin(9) ^ 2 * t83 + pkin(5) ^ 2) + 0.2e1 * pkin(9) * t79 + t72; mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,6) * t93 + t85; -t21 * t73 - t45; -t73 * t37; -t44; -pkin(9) * t73 + t45; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
