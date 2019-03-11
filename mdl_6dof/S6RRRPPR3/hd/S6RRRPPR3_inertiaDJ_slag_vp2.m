% Calculate time derivative of joint inertia matrix for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:36
% EndTime: 2019-03-09 15:28:41
% DurationCPUTime: 2.42s
% Computational Cost: add. (2917->279), mult. (6282->384), div. (0->0), fcn. (5467->6), ass. (0->125)
t172 = cos(qJ(3));
t99 = sin(qJ(2));
t130 = t172 * t99;
t173 = cos(qJ(2));
t98 = sin(qJ(3));
t72 = t173 * t98 + t130;
t199 = 0.2e1 * t72;
t198 = -mrSges(5,1) - mrSges(4,1);
t100 = cos(qJ(6));
t97 = sin(qJ(6));
t171 = mrSges(7,2) * t97;
t77 = mrSges(7,1) * t100 - t171;
t197 = -mrSges(6,1) - t77;
t155 = Ifges(7,4) * t100;
t114 = -Ifges(7,2) * t97 + t155;
t79 = -Ifges(7,1) * t97 - t155;
t196 = (t114 - t79) * qJD(6);
t145 = qJD(6) * t97;
t190 = qJD(2) + qJD(3);
t54 = t190 * t72;
t151 = t100 * t54;
t117 = t172 * t173;
t157 = t98 * t99;
t71 = -t117 + t157;
t109 = t71 * t145 - t151;
t90 = -t173 * pkin(2) - pkin(1);
t121 = -t71 * pkin(3) + t72 * qJ(4) - t90;
t177 = -pkin(4) - pkin(9);
t23 = pkin(5) * t72 + t177 * t71 + t121;
t176 = -pkin(8) - pkin(7);
t73 = t176 * t130;
t80 = t176 * t173;
t60 = -t80 * t98 - t73;
t37 = -qJ(5) * t72 + t60;
t11 = t100 * t23 - t37 * t97;
t107 = qJD(2) * t80;
t123 = qJD(3) * t172;
t187 = t190 * t157;
t27 = -t172 * t107 - t80 * t123 + t176 * t187;
t53 = -t117 * t190 + t187;
t18 = t53 * qJ(5) - t72 * qJD(5) + t27;
t149 = qJD(2) * t99;
t138 = pkin(2) * t149;
t22 = t54 * pkin(3) + t53 * qJ(4) - t72 * qJD(4) + t138;
t4 = -pkin(5) * t53 + t177 * t54 - t22;
t2 = qJD(6) * t11 + t100 * t18 + t4 * t97;
t12 = t100 * t37 + t23 * t97;
t3 = -qJD(6) * t12 + t100 * t4 - t18 * t97;
t195 = t100 * t2 - t3 * t97;
t156 = t100 ^ 2 + t97 ^ 2;
t194 = (-mrSges(4,2) * t172 + (-t156 * mrSges(7,3) + mrSges(6,2) + t198) * t98) * pkin(2) * qJD(3);
t193 = -0.2e1 * t121;
t192 = Ifges(3,1) - Ifges(3,2);
t144 = qJD(6) * t100;
t164 = t54 * t97;
t108 = t71 * t144 + t164;
t163 = t71 * t97;
t42 = -mrSges(7,2) * t72 - mrSges(7,3) * t163;
t150 = t100 * t71;
t43 = mrSges(7,1) * t72 - mrSges(7,3) * t150;
t191 = t100 * t42 - t97 * t43;
t19 = -mrSges(7,1) * t53 + mrSges(7,3) * t109;
t20 = mrSges(7,2) * t53 - mrSges(7,3) * t108;
t103 = m(7) * ((-t100 * t11 - t12 * t97) * qJD(6) + t195) + t100 * t20 - t97 * t19;
t189 = -t43 * t144 - t42 * t145 + t103;
t148 = qJD(3) * t98;
t137 = pkin(2) * t148;
t120 = t72 * t137;
t118 = pkin(2) * t123;
t81 = t118 + qJD(4);
t175 = pkin(2) * t98;
t87 = qJ(4) + t175;
t188 = t87 * t54 + t81 * t71 - t120;
t184 = 2 * m(5);
t183 = 2 * m(6);
t182 = 0.2e1 * m(7);
t61 = t157 * t176 - t172 * t80;
t38 = t71 * qJ(5) + t61;
t181 = 0.2e1 * t38;
t180 = -0.2e1 * t71;
t96 = qJ(4) + pkin(5);
t179 = 0.2e1 * t96;
t178 = m(5) + m(6);
t169 = Ifges(7,4) * t97;
t26 = t98 * t107 + t148 * t80 + t190 * t73;
t17 = -t54 * qJ(5) - t71 * qJD(5) - t26;
t167 = t17 * t38;
t166 = t38 * t81;
t165 = t53 * mrSges(6,3);
t161 = t81 * t87;
t147 = qJD(4) * t38;
t139 = mrSges(5,3) - t197;
t136 = t172 * pkin(2);
t126 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t125 = t61 * t26 + t27 * t60;
t124 = qJD(2) * t173;
t119 = t139 * t81;
t89 = -t136 - pkin(3);
t116 = mrSges(7,1) * t97 + mrSges(7,2) * t100;
t115 = Ifges(7,1) * t100 - t169;
t113 = t100 * t12 - t11 * t97;
t86 = -pkin(4) + t89;
t112 = t156 * t137;
t111 = -qJ(4) * t54 - qJD(4) * t71;
t110 = qJ(4) * t81 + qJD(4) * t87;
t76 = t115 * qJD(6);
t78 = -Ifges(7,2) * t100 - t169;
t106 = t196 * t100 + t78 * t145 + t97 * t76;
t104 = -t109 * Ifges(7,5) - Ifges(7,6) * t108 - Ifges(7,3) * t53;
t31 = Ifges(7,6) * t72 + t114 * t71;
t32 = Ifges(7,5) * t72 + t115 * t71;
t7 = -Ifges(7,4) * t109 - Ifges(7,2) * t108 - Ifges(7,6) * t53;
t74 = t116 * qJD(6);
t8 = -Ifges(7,1) * t109 - Ifges(7,4) * t108 - Ifges(7,5) * t53;
t91 = Ifges(7,6) * t145;
t102 = t72 * (-Ifges(7,5) * t144 + t91) / 0.2e1 + t18 * mrSges(6,2) + t31 * t145 / 0.2e1 - t76 * t150 / 0.2e1 + t79 * t151 / 0.2e1 - t78 * t164 / 0.2e1 - t38 * t74 - t97 * t8 / 0.2e1 - t100 * t7 / 0.2e1 + t198 * t27 + (mrSges(5,3) - mrSges(4,2)) * t26 + t197 * t17 + t196 * t163 / 0.2e1 - (t71 * t78 + t32) * t144 / 0.2e1 + (Ifges(5,6) - Ifges(6,5) - Ifges(4,6)) * t54 + (-Ifges(6,6) - Ifges(5,4) - Ifges(4,5) + Ifges(7,5) * t97 / 0.2e1 + Ifges(7,6) * t100 / 0.2e1) * t53 + (t11 * t144 + t12 * t145 - t195) * mrSges(7,3);
t101 = -pkin(3) - pkin(4);
t93 = -pkin(3) + t177;
t85 = pkin(5) + t87;
t82 = -pkin(9) + t86;
t39 = t116 * t71;
t33 = -pkin(4) * t71 + t121;
t16 = mrSges(7,1) * t108 - mrSges(7,2) * t109;
t13 = -pkin(4) * t54 - t22;
t1 = [0.2e1 * (mrSges(4,1) * t71 + mrSges(4,2) * t72) * t138 + 0.2e1 * m(4) * (t138 * t90 + t125) + (t17 * t180 - 0.2e1 * t18 * t72) * mrSges(6,3) + (mrSges(4,3) + mrSges(5,2)) * (t26 * t180 + t27 * t199 - 0.2e1 * t53 * t60 - 0.2e1 * t54 * t61) + 0.2e1 * t33 * t126 + t16 * t181 + (t11 * t3 + t12 * t2 - t167) * t182 + (t13 * t33 + t18 * t37 - t167) * t183 + t8 * t150 + (-0.2e1 * t90 * mrSges(4,2) + mrSges(5,3) * t193 + (-(2 * Ifges(4,1)) - (2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t72 + (-Ifges(7,5) * t100 + Ifges(7,6) * t97 + (2 * Ifges(4,4)) + (2 * Ifges(6,4)) - (2 * Ifges(5,5))) * t71) * t53 + (0.2e1 * mrSges(4,1) * t90 + mrSges(5,1) * t193 + mrSges(6,3) * t181 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t199 + 0.2e1 * (Ifges(6,1) + Ifges(5,3) + Ifges(4,2)) * t71) * t54 - t109 * t32 - t108 * t31 + t72 * t104 + (0.2e1 * Ifges(3,4) * t173 + t192 * t99) * t124 + (-0.2e1 * Ifges(3,4) * t99 + t173 * t192) * t149 - 0.2e1 * pkin(1) * (mrSges(3,1) * t99 + mrSges(3,2) * t173) * qJD(2) + 0.2e1 * t11 * t19 + 0.2e1 * t12 * t20 - 0.2e1 * t17 * t39 + 0.2e1 * t2 * t42 + 0.2e1 * t3 * t43 - t7 * t163 + 0.2e1 * t37 * t165 + 0.2e1 * t22 * (mrSges(5,1) * t71 - mrSges(5,3) * t72) + 0.2e1 * t13 * (mrSges(6,1) * t72 + mrSges(6,2) * t71) + (-t121 * t22 + t125) * t184; m(5) * (t26 * t87 + t27 * t89 + t61 * t81) + t86 * t165 + Ifges(3,5) * t124 + t102 + m(4) * (-t172 * t27 + t26 * t98 + (t172 * t61 + t60 * t98) * qJD(3)) * pkin(2) - Ifges(3,6) * t149 + m(6) * (-t17 * t87 + t18 * t86 + t166) + m(7) * (-t17 * t85 + t166) + t81 * t39 + t85 * t16 + (m(5) * t60 + m(6) * t37 + m(7) * t113 + t191) * t137 + (-mrSges(3,1) * t124 + mrSges(3,2) * t149) * pkin(7) + t188 * mrSges(6,3) + t189 * t82 + (-t118 * t71 + t136 * t53 - t175 * t54 + t120) * mrSges(4,3) + (-t53 * t89 - t188) * mrSges(5,2); -0.2e1 * t85 * t74 + 0.2e1 * t119 + 0.2e1 * t194 + (t137 * t86 + t161) * t183 + (t112 * t82 + t81 * t85) * t182 + (t137 * t89 + t161) * t184 + t106; t189 * t93 + m(7) * (-t17 * t96 + t147) + m(5) * (-pkin(3) * t27 + qJ(4) * t26 + qJD(4) * t61) + m(6) * (-qJ(4) * t17 + t101 * t18 + t147) + (t101 * t53 - t111) * mrSges(6,3) + (pkin(3) * t53 + t111) * mrSges(5,2) + t102 + qJD(4) * t39 + t96 * t16; -(t85 + t96) * t74 + t119 + t139 * qJD(4) + t194 + m(6) * (t101 * t137 + t110) + m(7) * (qJD(4) * t85 + t112 * t93 + t81 * t96) + m(5) * (-pkin(3) * t137 + t110) + t106; -t74 * t179 + (m(7) * t179 + 0.2e1 * qJ(4) * t178 + 0.2e1 * mrSges(6,1) + 0.2e1 * mrSges(5,3) + 0.2e1 * t77) * qJD(4) + t106; (-mrSges(5,2) + mrSges(6,3)) * t53 + (-t100 * t43 - t97 * t42) * qJD(6) + m(6) * t18 + m(5) * t27 + t103; (m(7) * t156 + t178) * t137; 0; 0; t100 * t19 + t97 * t20 + t191 * qJD(6) + m(7) * (qJD(6) * t113 + t100 * t3 + t2 * t97) + m(6) * t13 + t126; 0; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t104; t91 - t116 * t137 + (t82 * t171 + (-mrSges(7,1) * t82 - Ifges(7,5)) * t100) * qJD(6); t91 + (t93 * t171 + (-mrSges(7,1) * t93 - Ifges(7,5)) * t100) * qJD(6); -t77 * qJD(6); -t74; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
