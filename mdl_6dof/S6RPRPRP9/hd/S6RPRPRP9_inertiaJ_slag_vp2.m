% Calculate joint inertia matrix for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:45
% EndTime: 2019-03-09 03:27:47
% DurationCPUTime: 0.90s
% Computational Cost: add. (1017->242), mult. (1982->320), div. (0->0), fcn. (1845->6), ass. (0->98)
t130 = Ifges(7,2) + Ifges(6,3);
t90 = (-pkin(1) - pkin(7));
t129 = -2 * t90;
t128 = m(6) + m(7);
t127 = mrSges(7,2) + mrSges(6,3);
t118 = cos(qJ(5));
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t87 = sin(qJ(5));
t61 = t118 * t85 + t87 * t86;
t93 = t118 * t86 - t87 * t85;
t18 = -mrSges(7,1) * t93 - t61 * mrSges(7,3);
t19 = -mrSges(6,1) * t93 + t61 * mrSges(6,2);
t126 = -t18 - t19;
t89 = cos(qJ(3));
t46 = t61 * t89;
t48 = t93 * t89;
t11 = t46 * mrSges(7,1) - t48 * mrSges(7,3);
t12 = t46 * mrSges(6,1) + t48 * mrSges(6,2);
t114 = t86 * t89;
t115 = t85 * t89;
t49 = mrSges(5,1) * t115 + mrSges(5,2) * t114;
t125 = -t11 - t12 - t49;
t124 = -m(7) * pkin(5) - mrSges(7,1);
t123 = mrSges(6,1) - t124;
t122 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t121 = 2 * qJ(2);
t120 = t85 / 0.2e1;
t119 = t86 / 0.2e1;
t117 = Ifges(5,4) * t85;
t116 = Ifges(5,4) * t86;
t88 = sin(qJ(3));
t113 = t88 * t90;
t112 = t89 * t90;
t110 = pkin(8) + qJ(4);
t64 = pkin(3) * t88 - qJ(4) * t89 + qJ(2);
t57 = t86 * t64;
t15 = -pkin(8) * t114 + t57 + (-t85 * t90 + pkin(4)) * t88;
t33 = t86 * t113 + t85 * t64;
t17 = -pkin(8) * t115 + t33;
t4 = t118 * t17 + t87 * t15;
t28 = -mrSges(7,2) * t46 + mrSges(7,3) * t88;
t29 = -mrSges(6,2) * t88 - mrSges(6,3) * t46;
t109 = t28 + t29;
t30 = mrSges(6,1) * t88 - mrSges(6,3) * t48;
t31 = -t88 * mrSges(7,1) + t48 * mrSges(7,2);
t108 = -t30 + t31;
t65 = -t86 * mrSges(5,1) + t85 * mrSges(5,2);
t107 = -t65 + mrSges(4,1);
t106 = t85 ^ 2 + t86 ^ 2;
t82 = t88 ^ 2;
t83 = t89 ^ 2;
t105 = t83 + t82;
t104 = m(5) + t128;
t101 = t110 * t85;
t66 = t110 * t86;
t25 = t101 * t118 + t66 * t87;
t27 = -t101 * t87 + t118 * t66;
t103 = t25 ^ 2 + t27 ^ 2;
t74 = -pkin(4) * t86 - pkin(3);
t44 = t61 * t88;
t47 = t93 * t88;
t102 = t25 * t44 + t27 * t47;
t100 = t105 * mrSges(4,3);
t58 = pkin(4) * t115 - t112;
t98 = qJ(4) * t106;
t32 = -t113 * t85 + t57;
t97 = -t32 * t85 + t33 * t86;
t62 = -mrSges(5,2) * t88 - mrSges(5,3) * t115;
t63 = mrSges(5,1) * t88 - mrSges(5,3) * t114;
t95 = t86 * t62 - t85 * t63;
t94 = t130 * t88 + (Ifges(7,4) + Ifges(6,5)) * t48 + (-Ifges(6,6) + Ifges(7,6)) * t46;
t3 = t118 * t15 - t87 * t17;
t92 = qJ(2) ^ 2;
t84 = t90 ^ 2;
t76 = t83 * t90;
t75 = t83 * t84;
t68 = Ifges(5,1) * t85 + t116;
t67 = Ifges(5,2) * t86 + t117;
t55 = Ifges(7,4) * t61;
t54 = Ifges(6,5) * t61;
t53 = Ifges(6,6) * t93;
t52 = Ifges(7,6) * t93;
t42 = Ifges(5,5) * t88 + (Ifges(5,1) * t86 - t117) * t89;
t41 = Ifges(5,6) * t88 + (-Ifges(5,2) * t85 + t116) * t89;
t23 = Ifges(6,1) * t61 + Ifges(6,4) * t93;
t22 = Ifges(7,1) * t61 - Ifges(7,5) * t93;
t21 = Ifges(6,4) * t61 + Ifges(6,2) * t93;
t20 = Ifges(7,5) * t61 - Ifges(7,3) * t93;
t14 = -pkin(5) * t93 - qJ(6) * t61 + t74;
t10 = Ifges(6,1) * t48 - Ifges(6,4) * t46 + Ifges(6,5) * t88;
t9 = Ifges(7,1) * t48 + Ifges(7,4) * t88 + Ifges(7,5) * t46;
t8 = Ifges(6,4) * t48 - Ifges(6,2) * t46 + Ifges(6,6) * t88;
t7 = Ifges(7,5) * t48 + Ifges(7,6) * t88 + Ifges(7,3) * t46;
t5 = pkin(5) * t46 - qJ(6) * t48 + t58;
t2 = -t88 * pkin(5) - t3;
t1 = qJ(6) * t88 + t4;
t6 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t121) + 0.2e1 * t1 * t28 + 0.2e1 * t5 * t11 + 0.2e1 * t58 * t12 + 0.2e1 * t2 * t31 + 0.2e1 * t4 * t29 + 0.2e1 * t3 * t30 + 0.2e1 * t32 * t63 + 0.2e1 * t33 * t62 + Ifges(3,1) + Ifges(2,3) + (t9 + t10) * t48 + (t7 - t8) * t46 + t100 * t129 + ((mrSges(4,2) * t121) + Ifges(4,1) * t89 + t49 * t129 - t85 * t41 + t86 * t42) * t89 + (mrSges(4,1) * t121 + (Ifges(4,2) + Ifges(5,3)) * t88 + (Ifges(5,5) * t86 - Ifges(5,6) * t85 - (2 * Ifges(4,4))) * t89 + t94) * t88 + m(4) * (t82 * t84 + t75 + t92) + (m(3) * (pkin(1) ^ 2 + t92)) + m(5) * (t32 ^ 2 + t33 ^ 2 + t75) + m(6) * (t3 ^ 2 + t4 ^ 2 + t58 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); -(m(3) * pkin(1)) + mrSges(3,2) + t95 * t88 + t109 * t47 + t108 * t44 - t100 + t125 * t89 + m(7) * (t1 * t47 + t2 * t44 - t5 * t89) + m(6) * (-t3 * t44 + t4 * t47 - t58 * t89) + m(5) * (t88 * t97 + t76) + m(4) * (t82 * t90 + t76); m(3) + m(5) * (t106 * t82 + t83) + m(4) * t105 + t128 * (t44 ^ 2 + t47 ^ 2 + t83); t41 * t119 + t42 * t120 + t74 * t12 + t58 * t19 - pkin(3) * t49 + t14 * t11 + t5 * t18 + (-(t90 * mrSges(4,2)) + Ifges(5,5) * t120 + Ifges(5,6) * t119 + t54 / 0.2e1 + t53 / 0.2e1 + t55 / 0.2e1 - t52 / 0.2e1 - Ifges(4,6)) * t88 + (t22 / 0.2e1 + t23 / 0.2e1) * t48 + (t20 / 0.2e1 - t21 / 0.2e1) * t46 + t109 * t27 + t108 * t25 + t95 * qJ(4) + t97 * mrSges(5,3) + (t68 * t119 - t85 * t67 / 0.2e1 + Ifges(4,5) + t107 * t90) * t89 + m(5) * (pkin(3) * t112 + qJ(4) * t97) + m(6) * (-t25 * t3 + t27 * t4 + t58 * t74) + m(7) * (t1 * t27 + t14 * t5 + t2 * t25) + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t10 / 0.2e1 + t9 / 0.2e1) * t61 - (-t1 * mrSges(7,2) - t4 * mrSges(6,3) + t7 / 0.2e1 - t8 / 0.2e1) * t93; (mrSges(5,3) * t106 - mrSges(4,2)) * t88 + (t107 + t126) * t89 + m(6) * (-t74 * t89 + t102) + m(7) * (-t14 * t89 + t102) + m(5) * (pkin(3) * t89 + t88 * t98) + t127 * (t44 * t61 + t47 * t93); -0.2e1 * pkin(3) * t65 + 0.2e1 * t14 * t18 + 0.2e1 * t74 * t19 + t86 * t67 + t85 * t68 + Ifges(4,3) + m(7) * (t14 ^ 2 + t103) + m(6) * (t74 ^ 2 + t103) + m(5) * (qJ(4) ^ 2 * t106 + pkin(3) ^ 2) + (t22 + t23) * t61 - (t20 - t21) * t93 + 0.2e1 * mrSges(5,3) * t98 + 0.2e1 * (t25 * t61 + t27 * t93) * t127; -m(5) * t112 + m(6) * t58 + m(7) * t5 - t125; -t104 * t89; -m(5) * pkin(3) + m(6) * t74 + m(7) * t14 - t126 + t65; t104; -pkin(5) * t31 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t28 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t94; t122 * t47 - t123 * t44; -t52 + t54 + t53 + t55 + (-pkin(5) * t61 + qJ(6) * t93) * mrSges(7,2) + t122 * t27 - t123 * t25; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t130; m(7) * t2 + t31; m(7) * t44; m(7) * t25 + t61 * mrSges(7,2); 0; t124; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
