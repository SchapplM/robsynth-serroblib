% Calculate joint inertia matrix for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2018-11-23 15:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:58:15
% EndTime: 2018-11-23 15:58:16
% DurationCPUTime: 0.88s
% Computational Cost: add. (1075->240), mult. (2129->320), div. (0->0), fcn. (2000->8), ass. (0->97)
t86 = cos(pkin(10));
t89 = sin(qJ(3));
t115 = t86 * t89;
t84 = sin(pkin(10));
t116 = t84 * t89;
t50 = mrSges(5,1) * t116 + mrSges(5,2) * t115;
t120 = cos(qJ(5));
t88 = sin(qJ(5));
t60 = t120 * t84 + t88 * t86;
t44 = t60 * t89;
t94 = t120 * t86 - t88 * t84;
t46 = t94 * t89;
t14 = t44 * mrSges(7,1) - t46 * mrSges(7,3);
t15 = t44 * mrSges(6,1) + t46 * mrSges(6,2);
t93 = -t14 - t15;
t129 = -t50 + t93;
t85 = sin(pkin(9));
t75 = pkin(1) * t85 + pkin(7);
t128 = 0.2e1 * t75;
t127 = m(6) + m(7);
t126 = mrSges(7,2) + mrSges(6,3);
t20 = -mrSges(7,1) * t94 - t60 * mrSges(7,3);
t21 = -mrSges(6,1) * t94 + t60 * mrSges(6,2);
t125 = -t20 - t21;
t124 = -m(7) * pkin(5) - mrSges(7,1);
t123 = -t84 / 0.2e1;
t122 = t86 / 0.2e1;
t90 = cos(qJ(3));
t121 = pkin(3) * t90;
t87 = cos(pkin(9));
t76 = -pkin(1) * t87 - pkin(2);
t57 = -qJ(4) * t89 - t121 + t76;
t48 = t86 * t57;
t12 = -pkin(8) * t115 + t48 + (-t75 * t84 - pkin(4)) * t90;
t117 = t75 * t90;
t19 = t86 * t117 + t84 * t57;
t16 = -pkin(8) * t116 + t19;
t4 = t88 * t12 + t120 * t16;
t119 = Ifges(5,4) * t84;
t118 = Ifges(5,4) * t86;
t114 = t89 * mrSges(5,3);
t69 = t89 * t75;
t112 = Ifges(7,2) + Ifges(6,3);
t111 = pkin(8) + qJ(4);
t30 = -mrSges(7,2) * t44 - mrSges(7,3) * t90;
t31 = mrSges(6,2) * t90 - mrSges(6,3) * t44;
t110 = t30 + t31;
t32 = -mrSges(6,1) * t90 - mrSges(6,3) * t46;
t33 = t90 * mrSges(7,1) + t46 * mrSges(7,2);
t109 = -t32 + t33;
t64 = -t86 * mrSges(5,1) + t84 * mrSges(5,2);
t108 = t64 - mrSges(4,1);
t49 = pkin(4) * t116 + t69;
t107 = t84 ^ 2 + t86 ^ 2;
t82 = t89 ^ 2;
t83 = t90 ^ 2;
t106 = t82 + t83;
t105 = m(5) + t127;
t102 = t111 * t84;
t65 = t111 * t86;
t27 = t102 * t120 + t65 * t88;
t29 = -t102 * t88 + t120 * t65;
t104 = t27 ^ 2 + t29 ^ 2;
t77 = -pkin(4) * t86 - pkin(3);
t103 = t27 * t44 + t29 * t46;
t100 = qJ(4) * t107;
t99 = (-Ifges(7,4) - Ifges(6,5)) * t46 + (Ifges(6,6) - Ifges(7,6)) * t44;
t98 = -pkin(5) * t44 + qJ(6) * t46;
t18 = -t117 * t84 + t48;
t97 = -t18 * t84 + t19 * t86;
t61 = mrSges(5,2) * t90 - t114 * t84;
t62 = -mrSges(5,1) * t90 - t114 * t86;
t95 = t86 * t61 - t84 * t62;
t3 = t12 * t120 - t88 * t16;
t73 = t75 ^ 2;
t68 = t82 * t73;
t67 = Ifges(5,1) * t84 + t118;
t66 = Ifges(5,2) * t86 + t119;
t56 = Ifges(7,4) * t60;
t55 = Ifges(6,5) * t60;
t54 = Ifges(6,6) * t94;
t53 = Ifges(7,6) * t94;
t42 = -Ifges(5,5) * t90 + (Ifges(5,1) * t86 - t119) * t89;
t41 = -Ifges(5,6) * t90 + (-Ifges(5,2) * t84 + t118) * t89;
t25 = Ifges(6,1) * t60 + Ifges(6,4) * t94;
t24 = Ifges(7,1) * t60 - Ifges(7,5) * t94;
t23 = Ifges(6,4) * t60 + Ifges(6,2) * t94;
t22 = Ifges(7,5) * t60 - Ifges(7,3) * t94;
t17 = -pkin(5) * t94 - qJ(6) * t60 + t77;
t11 = Ifges(6,1) * t46 - Ifges(6,4) * t44 - Ifges(6,5) * t90;
t10 = Ifges(7,1) * t46 - Ifges(7,4) * t90 + Ifges(7,5) * t44;
t9 = Ifges(6,4) * t46 - Ifges(6,2) * t44 - Ifges(6,6) * t90;
t8 = Ifges(7,5) * t46 - Ifges(7,6) * t90 + Ifges(7,3) * t44;
t5 = -t98 + t49;
t2 = t90 * pkin(5) - t3;
t1 = -qJ(6) * t90 + t4;
t6 = [0.2e1 * t1 * t30 + 0.2e1 * t5 * t14 + 0.2e1 * t49 * t15 + 0.2e1 * t18 * t62 + 0.2e1 * t19 * t61 + 0.2e1 * t2 * t33 + 0.2e1 * t3 * t32 + 0.2e1 * t4 * t31 + Ifges(2,3) + Ifges(3,3) + (t10 + t11) * t46 + (t8 - t9) * t44 + (0.2e1 * t76 * mrSges(4,2) + Ifges(4,1) * t89 + t50 * t128 - t84 * t41 + t86 * t42) * t89 + (-0.2e1 * t76 * mrSges(4,1) + (Ifges(5,3) + Ifges(4,2) + t112) * t90 + (-Ifges(5,5) * t86 + Ifges(5,6) * t84 + (2 * Ifges(4,4))) * t89 + t99) * t90 + m(4) * (t73 * t83 + t76 ^ 2 + t68) + m(5) * (t18 ^ 2 + t19 ^ 2 + t68) + m(6) * (t3 ^ 2 + t4 ^ 2 + t49 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(3) * (t85 ^ 2 + t87 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t87 - mrSges(3,2) * t85) * pkin(1) + t106 * mrSges(4,3) * t128; t110 * t46 + t109 * t44 + t129 * t90 + m(6) * (-t3 * t44 + t4 * t46 - t49 * t90) + m(7) * (t1 * t46 + t2 * t44 - t5 * t90) + (m(5) * (t97 - t117) + t95) * t89; m(3) + m(5) * (t107 * t82 + t83) + m(4) * t106 + t127 * (t44 ^ 2 + t46 ^ 2 + t83); t84 * t42 / 0.2e1 + t41 * t122 + t77 * t15 + t49 * t21 - pkin(3) * t50 + t5 * t20 + t17 * t14 + (Ifges(5,5) * t123 - Ifges(5,6) * t86 / 0.2e1 - t55 / 0.2e1 - t54 / 0.2e1 - t56 / 0.2e1 + t53 / 0.2e1 + Ifges(4,6) - t75 * mrSges(4,2)) * t90 + (t24 / 0.2e1 + t25 / 0.2e1) * t46 + (t22 / 0.2e1 - t23 / 0.2e1) * t44 + t110 * t29 + t109 * t27 + t95 * qJ(4) + t97 * mrSges(5,3) + (t108 * t75 + t67 * t122 + t123 * t66 + Ifges(4,5)) * t89 + m(5) * (-pkin(3) * t69 + qJ(4) * t97) + m(6) * (-t27 * t3 + t29 * t4 + t49 * t77) + m(7) * (t1 * t29 + t17 * t5 + t2 * t27) + (t10 / 0.2e1 + t11 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3)) * t60 - (t8 / 0.2e1 - t9 / 0.2e1 - t4 * mrSges(6,3) - t1 * mrSges(7,2)) * t94; (mrSges(5,3) * t107 - mrSges(4,2)) * t89 + (-t108 + t125) * t90 + m(6) * (-t77 * t90 + t103) + m(7) * (-t17 * t90 + t103) + m(5) * (t100 * t89 + t121) + t126 * (t44 * t60 + t46 * t94); -0.2e1 * pkin(3) * t64 + 0.2e1 * t17 * t20 + 0.2e1 * t77 * t21 + t86 * t66 + t84 * t67 + Ifges(4,3) + m(7) * (t17 ^ 2 + t104) + m(6) * (t77 ^ 2 + t104) + m(5) * (qJ(4) ^ 2 * t107 + pkin(3) ^ 2) + (t24 + t25) * t60 - (t22 - t23) * t94 + 0.2e1 * mrSges(5,3) * t100 + 0.2e1 * (t27 * t60 + t29 * t94) * t126; m(5) * t69 + m(6) * t49 + m(7) * t5 - t129; -t105 * t90; -m(5) * pkin(3) + m(6) * t77 + m(7) * t17 - t125 + t64; t105; m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t30 + t1 * mrSges(7,3) - pkin(5) * t33 - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t112 * t90 - t99; m(7) * t98 + t93; t56 - t53 + t55 + t54 + (-pkin(5) * t60 + qJ(6) * t94) * mrSges(7,2) + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t29 + (-mrSges(6,1) + t124) * t27; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t112; m(7) * t2 + t33; m(7) * t44; m(7) * t27 + t60 * mrSges(7,2); 0; t124; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
