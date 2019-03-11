% Calculate time derivative of joint inertia matrix for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:36
% EndTime: 2019-03-09 03:27:41
% DurationCPUTime: 2.57s
% Computational Cost: add. (2374->357), mult. (5455->520), div. (0->0), fcn. (4542->6), ass. (0->150)
t204 = Ifges(7,4) + Ifges(6,5);
t202 = -Ifges(6,6) + Ifges(7,6);
t203 = Ifges(7,2) + Ifges(6,3);
t130 = (-pkin(1) - pkin(7));
t201 = 2 * t130;
t200 = m(6) + m(7);
t199 = mrSges(6,3) + mrSges(7,2);
t124 = sin(pkin(9));
t125 = cos(pkin(9));
t126 = sin(qJ(5));
t128 = cos(qJ(5));
t103 = t124 * t126 - t128 * t125;
t129 = cos(qJ(3));
t127 = sin(qJ(3));
t158 = qJD(3) * t127;
t163 = t125 * t126;
t104 = t124 * t128 + t163;
t192 = t104 * qJD(5);
t43 = t103 * t158 - t129 * t192;
t151 = t124 * t158;
t154 = qJD(5) * t128;
t155 = qJD(5) * t126;
t196 = -t124 * t155 + t125 * t154;
t45 = -t128 * t151 + t196 * t129 - t158 * t163;
t10 = t45 * mrSges(7,1) - t43 * mrSges(7,3);
t11 = t45 * mrSges(6,1) + t43 * mrSges(6,2);
t198 = -t10 - t11;
t197 = m(7) * qJD(6);
t182 = pkin(3) * t127;
t107 = -qJ(4) * t129 + qJ(2) + t182;
t100 = t125 * t107;
t160 = t127 * t130;
t77 = -t124 * t160 + t100;
t78 = t124 * t107 + t125 * t160;
t195 = -t124 * t77 + t125 * t78;
t156 = qJD(3) * t130;
t149 = t129 * t156;
t85 = -qJD(4) * t129 + qJD(2) + (pkin(3) * t129 + qJ(4) * t127) * qJD(3);
t80 = t125 * t85;
t60 = -t124 * t149 + t80;
t61 = t124 * t85 + t125 * t149;
t194 = -t124 * t60 + t125 * t61;
t94 = t103 * qJD(5);
t193 = t202 * t192 - t204 * t94;
t191 = m(7) * qJ(6) + mrSges(7,3);
t144 = -t124 * t130 + pkin(4);
t161 = t125 * t129;
t58 = -pkin(8) * t161 + t127 * t144 + t100;
t164 = t124 * t129;
t62 = -pkin(8) * t164 + t78;
t177 = t126 * t58 + t128 * t62;
t162 = t125 * t127;
t27 = t80 + (pkin(8) * t162 + t129 * t144) * qJD(3);
t46 = pkin(8) * t151 + t61;
t4 = -qJD(5) * t177 - t126 * t46 + t128 * t27;
t190 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t189 = -mrSges(6,2) + t191;
t188 = 2 * m(5);
t187 = 2 * m(6);
t186 = 0.2e1 * m(7);
t185 = 2 * mrSges(4,1);
t123 = t125 ^ 2;
t184 = t124 / 0.2e1;
t181 = pkin(8) + qJ(4);
t157 = qJD(3) * t129;
t23 = -mrSges(7,2) * t45 + mrSges(7,3) * t157;
t26 = -mrSges(6,2) * t157 - mrSges(6,3) * t45;
t180 = t23 + t26;
t24 = mrSges(6,1) * t157 - mrSges(6,3) * t43;
t25 = -mrSges(7,1) * t157 + t43 * mrSges(7,2);
t179 = -t24 + t25;
t51 = mrSges(7,1) * t192 + t94 * mrSges(7,3);
t52 = mrSges(6,1) * t192 - t94 * mrSges(6,2);
t178 = -t51 - t52;
t82 = t104 * t129;
t71 = -mrSges(7,2) * t82 + mrSges(7,3) * t127;
t72 = -mrSges(6,2) * t127 - mrSges(6,3) * t82;
t176 = t71 + t72;
t84 = t103 * t129;
t73 = mrSges(6,1) * t127 + mrSges(6,3) * t84;
t74 = -mrSges(7,1) * t127 - mrSges(7,2) * t84;
t175 = -t73 + t74;
t172 = Ifges(5,4) * t124;
t171 = Ifges(5,4) * t125;
t170 = t124 * Ifges(5,2);
t165 = -mrSges(5,1) * t125 + mrSges(5,2) * t124 - mrSges(4,1);
t159 = t124 ^ 2 + t123;
t153 = t126 * qJD(4);
t152 = t128 * qJD(4);
t119 = -pkin(4) * t125 - pkin(3);
t150 = t127 * t157;
t113 = t181 * t125;
t143 = qJD(5) * t181;
t28 = t125 * t152 - t113 * t155 + (-t128 * t143 - t153) * t124;
t29 = t125 * t153 + t113 * t154 + (-t126 * t143 + t152) * t124;
t145 = t181 * t124;
t69 = t126 * t113 + t128 * t145;
t70 = t128 * t113 - t126 * t145;
t146 = t70 * t28 + t29 * t69;
t142 = t159 * mrSges(5,3);
t141 = t159 * qJ(4);
t101 = pkin(4) * t164 - t129 * t130;
t140 = mrSges(5,1) * t124 + mrSges(5,2) * t125;
t139 = -Ifges(5,5) * t125 + Ifges(5,6) * t124;
t15 = -t126 * t62 + t128 * t58;
t136 = t203 * t157 + t202 * t45 + t204 * t43;
t93 = -pkin(4) * t151 + t127 * t156;
t3 = t126 * t27 + t128 * t46 + t58 * t154 - t62 * t155;
t42 = -qJD(3) * t84 - t127 * t192;
t44 = qJD(3) * t82 + t196 * t127;
t81 = t104 * t127;
t83 = t103 * t127;
t135 = -t28 * t83 + t29 * t81 + t70 * t42 + t44 * t69;
t106 = mrSges(5,1) * t127 - mrSges(5,3) * t161;
t105 = -mrSges(5,2) * t127 - mrSges(5,3) * t164;
t98 = (mrSges(5,1) * t129 + mrSges(5,3) * t162) * qJD(3);
t97 = (mrSges(5,3) * t124 * t127 - mrSges(5,2) * t129) * qJD(3);
t96 = t140 * t129;
t86 = t140 * t158;
t76 = (t129 * Ifges(5,5) + (-t125 * Ifges(5,1) + t172) * t127) * qJD(3);
t75 = (t129 * Ifges(5,6) + (t170 - t171) * t127) * qJD(3);
t68 = Ifges(6,1) * t104 - Ifges(6,4) * t103;
t67 = Ifges(7,1) * t104 + Ifges(7,5) * t103;
t66 = Ifges(6,4) * t104 - Ifges(6,2) * t103;
t65 = Ifges(7,5) * t104 + Ifges(7,3) * t103;
t64 = mrSges(6,1) * t103 + mrSges(6,2) * t104;
t63 = mrSges(7,1) * t103 - mrSges(7,3) * t104;
t57 = pkin(5) * t103 - qJ(6) * t104 + t119;
t56 = -Ifges(6,1) * t94 - Ifges(6,4) * t192;
t55 = -Ifges(7,1) * t94 + Ifges(7,5) * t192;
t54 = -Ifges(6,4) * t94 - Ifges(6,2) * t192;
t53 = -Ifges(7,5) * t94 + Ifges(7,3) * t192;
t48 = mrSges(6,1) * t82 - mrSges(6,2) * t84;
t47 = mrSges(7,1) * t82 + mrSges(7,3) * t84;
t33 = -Ifges(6,1) * t84 - Ifges(6,4) * t82 + Ifges(6,5) * t127;
t32 = -Ifges(7,1) * t84 + Ifges(7,4) * t127 + Ifges(7,5) * t82;
t31 = -Ifges(6,4) * t84 - Ifges(6,2) * t82 + Ifges(6,6) * t127;
t30 = -Ifges(7,5) * t84 + Ifges(7,6) * t127 + Ifges(7,3) * t82;
t21 = pkin(5) * t192 + qJ(6) * t94 - qJD(6) * t104;
t20 = t82 * pkin(5) + t84 * qJ(6) + t101;
t13 = -pkin(5) * t127 - t15;
t12 = qJ(6) * t127 + t177;
t9 = Ifges(6,1) * t43 - Ifges(6,4) * t45 + Ifges(6,5) * t157;
t8 = Ifges(7,1) * t43 + Ifges(7,4) * t157 + Ifges(7,5) * t45;
t7 = Ifges(6,4) * t43 - Ifges(6,2) * t45 + Ifges(6,6) * t157;
t6 = Ifges(7,5) * t43 + Ifges(7,6) * t157 + Ifges(7,3) * t45;
t5 = pkin(5) * t45 - qJ(6) * t43 + qJD(6) * t84 + t93;
t2 = -pkin(5) * t157 - t4;
t1 = qJ(6) * t157 + qJD(6) * t127 + t3;
t14 = [(t127 * t185 + 0.2e1 * t129 * mrSges(4,2) + (2 * mrSges(3,3)) + 0.2e1 * (m(3) + m(4)) * qJ(2)) * qJD(2) + 0.2e1 * t177 * t26 + (t101 * t93 + t15 * t4 + t177 * t3) * t187 + ((qJ(2) * t185 + (-(2 * Ifges(4,4)) - t139) * t129 - t204 * t84 + t202 * t82) * t129 + (t96 * t201 - 0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * (Ifges(4,4) + t139) * t127 + (-Ifges(5,1) * t123 + (2 * Ifges(4,2)) - (2 * Ifges(4,1)) + (2 * Ifges(5,3)) - (2 * m(5) * t130 ^ 2) + (-t170 + 0.2e1 * t171) * t124 + t203) * t129) * t127) * qJD(3) - (t8 + t9) * t84 + (t6 - t7) * t82 + 0.2e1 * t15 * t24 + 0.2e1 * t13 * t25 + 0.2e1 * t20 * t10 + 0.2e1 * t12 * t23 + (t1 * t12 + t13 * t2 + t20 * t5) * t186 + (t30 - t31) * t45 + t136 * t127 + (t32 + t33) * t43 + (-t124 * t75 + t125 * t76 + t86 * t201) * t129 + (t77 * t60 + t78 * t61) * t188 + 0.2e1 * t5 * t47 + 0.2e1 * t1 * t71 + 0.2e1 * t3 * t72 + 0.2e1 * t4 * t73 + 0.2e1 * t2 * t74 + 0.2e1 * t93 * t48 + 0.2e1 * t78 * t97 + 0.2e1 * t77 * t98 + 0.2e1 * t101 * t11 + 0.2e1 * t61 * t105 + 0.2e1 * t60 * t106; -t180 * t83 + t179 * t81 + t175 * t44 + t176 * t42 + (-t124 * t98 + t125 * t97) * t127 + (t86 + t198) * t129 + ((t105 * t125 - t106 * t124) * t129 + (t47 + t48 + t96) * t127) * qJD(3) + m(7) * (-t1 * t83 + t12 * t42 - t129 * t5 + t13 * t44 + t20 * t158 + t2 * t81) + m(6) * (t101 * t158 - t129 * t93 - t15 * t44 + t177 * t42 - t3 * t83 - t4 * t81) + m(5) * (t194 * t127 + (-0.2e1 * t160 + t195) * t157); 0.2e1 * m(5) * (-0.1e1 + t159) * t150 + 0.2e1 * t200 * (-t83 * t42 + t44 * t81 - t150); ((-t130 * mrSges(4,2) + Ifges(5,5) * t184 + Ifges(5,6) * t125 / 0.2e1 - Ifges(4,6) + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t104 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t103) * t129 + ((Ifges(5,2) * t125 + t172) * t184 - t125 * (Ifges(5,1) * t124 + t171) / 0.2e1 - Ifges(4,5) + (-m(5) * pkin(3) + t165) * t130) * t127) * qJD(3) + t179 * t69 + t180 * t70 + t175 * t29 + t176 * t28 + m(6) * (t119 * t93 - t15 * t29 + t177 * t28 + t3 * t70 - t4 * t69) + m(7) * (t1 * t70 + t12 * t28 + t13 * t29 + t2 * t69 + t20 * t21 + t5 * t57) - (t55 / 0.2e1 + t56 / 0.2e1) * t84 + (t53 / 0.2e1 - t54 / 0.2e1) * t82 + (t65 / 0.2e1 - t66 / 0.2e1) * t45 + (t67 / 0.2e1 + t68 / 0.2e1) * t43 + (t61 * mrSges(5,3) + qJ(4) * t97 + qJD(4) * t105 + t75 / 0.2e1) * t125 + (-t60 * mrSges(5,3) - qJ(4) * t98 - qJD(4) * t106 + t76 / 0.2e1) * t124 + (-t103 * t3 - t104 * t4 + t15 * t94 - t177 * t192) * mrSges(6,3) + (-t1 * t103 + t104 * t2 - t12 * t192 - t13 * t94) * mrSges(7,2) + (t30 / 0.2e1 - t31 / 0.2e1) * t192 - (t32 / 0.2e1 + t33 / 0.2e1) * t94 + (t8 / 0.2e1 + t9 / 0.2e1) * t104 + (t6 / 0.2e1 - t7 / 0.2e1) * t103 + t21 * t47 + t20 * t51 + t57 * t10 + t5 * t63 + pkin(3) * t86 + t93 * t64 + t101 * t52 + t119 * t11 + t193 * t127 / 0.2e1 + m(5) * (t194 * qJ(4) + t195 * qJD(4)); t178 * t129 + ((-mrSges(4,2) + t142) * t129 + (t63 + t64 + t165) * t127) * qJD(3) + m(6) * (t119 * t158 + t135) + m(7) * (-t129 * t21 + t57 * t158 + t135) + m(5) * (t159 * t127 * qJD(4) + (t129 * t141 - t182) * qJD(3)) + t199 * (-t103 * t42 + t104 * t44 + t192 * t83 - t81 * t94); 0.2e1 * t119 * t52 + 0.2e1 * t21 * t63 + 0.2e1 * t57 * t51 + (t65 - t66) * t192 - (t68 + t67) * t94 + (t56 + t55) * t104 + (t53 - t54) * t103 + (t21 * t57 + t146) * t186 + t146 * t187 + (t141 * t188 + 0.2e1 * t142) * qJD(4) + 0.2e1 * t199 * (-t103 * t28 + t104 * t29 - t192 * t70 - t69 * t94); m(6) * t93 + m(7) * t5 + ((m(5) * t130) - t140) * t158 - t198; (m(5) + t200) * t158; m(7) * t21 - t178; 0; -pkin(5) * t25 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + qJD(6) * t71 + qJ(6) * t23 + t1 * mrSges(7,3) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) + t136; t189 * t42 + t190 * t44 - t83 * t197; t70 * t197 + (pkin(5) * t94 - qJ(6) * t192 - qJD(6) * t103) * mrSges(7,2) + t190 * t29 + t189 * t28 + t193; 0; 0.2e1 * t191 * qJD(6); m(7) * t2 + t25; m(7) * t44; m(7) * t29 - t94 * mrSges(7,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
