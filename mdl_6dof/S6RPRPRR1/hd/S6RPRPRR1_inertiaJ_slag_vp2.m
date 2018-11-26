% Calculate joint inertia matrix for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 16:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:06
% EndTime: 2018-11-23 16:03:06
% DurationCPUTime: 0.77s
% Computational Cost: add. (1566->188), mult. (2850->275), div. (0->0), fcn. (3116->10), ass. (0->80)
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t97 = t73 ^ 2 + t76 ^ 2;
t123 = mrSges(7,3) * t97;
t77 = cos(qJ(3));
t122 = t77 ^ 2;
t121 = m(5) * pkin(3);
t103 = t73 * mrSges(7,3);
t112 = cos(qJ(5));
t69 = sin(pkin(11));
t71 = cos(pkin(11));
t75 = sin(qJ(3));
t49 = -t69 * t75 + t71 * t77;
t50 = t69 * t77 + t71 * t75;
t74 = sin(qJ(5));
t36 = -t112 * t49 + t50 * t74;
t38 = t112 * t50 + t74 * t49;
t15 = -mrSges(7,2) * t36 - t38 * t103;
t106 = t38 * t76;
t16 = mrSges(7,1) * t36 - mrSges(7,3) * t106;
t85 = t76 * t15 - t73 * t16;
t54 = -mrSges(7,1) * t76 + mrSges(7,2) * t73;
t120 = t123 * t38 + t36 * t54;
t70 = sin(pkin(10));
t60 = pkin(1) * t70 + pkin(7);
t95 = qJ(4) + t60;
t89 = t95 * t77;
t90 = t95 * t75;
t29 = -t69 * t90 + t71 * t89;
t18 = pkin(8) * t49 + t29;
t28 = -t69 * t89 - t71 * t90;
t82 = -t50 * pkin(8) + t28;
t10 = -t112 * t82 + t18 * t74;
t119 = t10 ^ 2;
t118 = t36 ^ 2;
t117 = t49 ^ 2;
t116 = 0.2e1 * t49;
t115 = pkin(3) * t69;
t114 = pkin(5) * t36;
t12 = t112 * t18 + t74 * t82;
t72 = cos(pkin(10));
t62 = -pkin(1) * t72 - pkin(2);
t53 = -pkin(3) * t77 + t62;
t39 = -pkin(4) * t49 + t53;
t13 = -pkin(9) * t38 + t114 + t39;
t3 = t12 * t76 + t13 * t73;
t113 = t3 * t76;
t110 = Ifges(7,4) * t73;
t109 = Ifges(7,4) * t76;
t108 = t10 * t36;
t107 = t38 * t73;
t61 = pkin(3) * t71 + pkin(4);
t45 = t112 * t61 - t74 * t115;
t105 = t45 * mrSges(6,1);
t46 = t112 * t115 + t74 * t61;
t104 = t46 * mrSges(6,2);
t100 = Ifges(7,5) * t106 + Ifges(7,3) * t36;
t99 = t36 * mrSges(6,1) + t38 * mrSges(6,2);
t98 = Ifges(7,5) * t73 + Ifges(7,6) * t76;
t96 = t75 ^ 2 + t122;
t55 = Ifges(7,2) * t76 + t110;
t56 = Ifges(7,1) * t73 + t109;
t94 = t76 * t55 + t73 * t56 + Ifges(6,3);
t93 = t97 * pkin(9);
t44 = pkin(9) + t46;
t92 = t97 * t44;
t91 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t2 = -t12 * t73 + t13 * t76;
t88 = -t2 * t73 + t113;
t87 = -t77 * mrSges(4,1) + t75 * mrSges(4,2);
t86 = -mrSges(7,1) * t73 - mrSges(7,2) * t76;
t84 = 0.2e1 * t123;
t83 = t91 + t99;
t8 = Ifges(7,6) * t36 + (-Ifges(7,2) * t73 + t109) * t38;
t9 = Ifges(7,5) * t36 + (Ifges(7,1) * t76 - t110) * t38;
t81 = -t12 * mrSges(6,2) + mrSges(7,3) * t113 - t2 * t103 - t55 * t107 / 0.2e1 + t56 * t106 / 0.2e1 + Ifges(6,5) * t38 + t73 * t9 / 0.2e1 + t76 * t8 / 0.2e1 + (t98 / 0.2e1 - Ifges(6,6)) * t36 + (t54 - mrSges(6,1)) * t10;
t43 = -pkin(5) - t45;
t35 = t38 ^ 2;
t14 = -mrSges(7,1) * t107 - mrSges(7,2) * t106;
t1 = [0.2e1 * t62 * t87 + Ifges(4,2) * t122 + Ifges(5,2) * t117 + 0.2e1 * t53 * t91 + 0.2e1 * t39 * t99 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + t29 * mrSges(5,3) * t116 - 0.2e1 * t10 * t14 + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * mrSges(5,3) * t28 + Ifges(5,1) * t50 + Ifges(5,4) * t116) * t50 + (-0.2e1 * mrSges(6,3) * t12 + Ifges(6,2) * t36 + t100) * t36 + (0.2e1 * mrSges(6,3) * t10 + Ifges(6,1) * t38 - t73 * t8 + t76 * t9 + (-Ifges(7,6) * t73 - (2 * Ifges(6,4))) * t36) * t38 + m(4) * (t96 * t60 ^ 2 + t62 ^ 2) + m(5) * (t28 ^ 2 + t29 ^ 2 + t53 ^ 2) + m(6) * (t12 ^ 2 + t39 ^ 2 + t119) + m(7) * (t2 ^ 2 + t3 ^ 2 + t119) + m(3) * (t70 ^ 2 + t72 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t75 + 0.2e1 * Ifges(4,4) * t77) * t75 + 0.2e1 * (mrSges(3,1) * t72 - mrSges(3,2) * t70) * pkin(1) + 0.2e1 * t96 * t60 * mrSges(4,3); -t36 * t14 + t85 * t38 + m(7) * (t88 * t38 + t108) + m(6) * (t12 * t38 + t108) + m(5) * (t28 * t49 + t29 * t50); m(3) + m(7) * (t97 * t35 + t118) + m(6) * (t35 + t118) + m(4) * t96 + m(5) * (t50 ^ 2 + t117); t85 * t44 + m(6) * (-t10 * t45 + t12 * t46) + (-mrSges(4,1) * t75 - mrSges(4,2) * t77) * t60 + (-t36 * t46 - t38 * t45) * mrSges(6,3) + Ifges(4,5) * t75 + Ifges(4,6) * t77 - t43 * t14 + Ifges(5,6) * t49 + Ifges(5,5) * t50 + t28 * mrSges(5,1) - t29 * mrSges(5,2) + t81 + m(7) * (t10 * t43 + t88 * t44) + (m(5) * (t28 * t71 + t29 * t69) + (t49 * t69 - t50 * t71) * mrSges(5,3)) * pkin(3); m(7) * (t36 * t43 + t38 * t92) + m(6) * (-t36 * t45 + t38 * t46) + (t49 * t71 + t50 * t69) * t121 - t83 - t87 + t120; 0.2e1 * t105 - 0.2e1 * t104 + 0.2e1 * t43 * t54 + Ifges(4,3) + Ifges(5,3) + t44 * t84 + m(7) * (t97 * t44 ^ 2 + t43 ^ 2) + m(6) * (t45 ^ 2 + t46 ^ 2) + t94 + (0.2e1 * mrSges(5,1) * t71 - 0.2e1 * mrSges(5,2) * t69 + (t69 ^ 2 + t71 ^ 2) * t121) * pkin(3); t73 * t15 + t76 * t16 + m(7) * (t2 * t76 + t3 * t73) + m(6) * t39 + m(5) * t53 + t83; 0; 0; m(7) * t97 + m(5) + m(6); t81 + (-m(7) * t10 + t14) * pkin(5) + (m(7) * t88 + t85) * pkin(9); m(7) * (t38 * t93 - t114) - t99 + t120; m(7) * (-pkin(5) * t43 + pkin(9) * t92) - t104 + t105 + (t43 - pkin(5)) * t54 + (t92 + t93) * mrSges(7,3) + t94; 0; -0.2e1 * pkin(5) * t54 + m(7) * (t97 * pkin(9) ^ 2 + pkin(5) ^ 2) + pkin(9) * t84 + t94; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t107 + t100; t14; t86 * t44 + t98; -t54; t86 * pkin(9) + t98; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
