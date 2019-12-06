% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:03
% EndTime: 2019-12-05 16:12:09
% DurationCPUTime: 2.06s
% Computational Cost: add. (1816->273), mult. (4819->404), div. (0->0), fcn. (4001->6), ass. (0->149)
t125 = cos(pkin(8));
t124 = sin(pkin(8));
t127 = sin(qJ(2));
t128 = cos(qJ(3));
t129 = cos(qJ(2));
t171 = t128 * t129;
t78 = t127 * t124 + t125 * t171;
t185 = t78 * t125;
t174 = t127 * t125;
t77 = t124 * t171 - t174;
t186 = t77 * t124;
t140 = (t185 + t186) * qJ(4);
t126 = sin(qJ(3));
t175 = t126 * t129;
t210 = -m(6) / 0.2e1;
t212 = -m(5) / 0.2e1;
t221 = mrSges(5,3) + mrSges(6,2);
t95 = -pkin(4) * t125 - qJ(5) * t124 - pkin(3);
t224 = -t221 * (t185 / 0.2e1 + t186 / 0.2e1) + (-pkin(3) * t175 + t140) * t212 + (t175 * t95 + t140) * t210;
t222 = m(5) + m(6);
t220 = Ifges(6,4) + Ifges(5,5);
t219 = Ifges(5,6) - Ifges(6,6);
t180 = t124 * t126;
t87 = t128 * mrSges(5,2) - mrSges(5,3) * t180;
t94 = -mrSges(6,2) * t180 - t128 * mrSges(6,3);
t218 = t87 + t94;
t178 = t125 * t126;
t89 = -t128 * mrSges(5,1) - mrSges(5,3) * t178;
t90 = t128 * mrSges(6,1) + mrSges(6,2) * t178;
t217 = -t89 + t90;
t98 = -mrSges(6,1) * t125 - mrSges(6,3) * t124;
t160 = m(6) * t95 + t98;
t121 = t125 ^ 2;
t170 = t124 ^ 2 + t121;
t183 = qJ(4) * t128;
t100 = pkin(3) * t126 - t183;
t182 = t100 * t125;
t58 = pkin(6) * t180 + t182;
t86 = t124 * t100;
t59 = -pkin(6) * t178 + t86;
t215 = -t124 * t58 + t125 * t59;
t41 = t86 + (-pkin(6) * t125 + qJ(5)) * t126;
t159 = pkin(6) * t124 + pkin(4);
t42 = -t126 * t159 - t182;
t214 = t124 * t42 + t125 * t41;
t213 = t210 + t212;
t209 = m(6) / 0.2e1;
t211 = m(5) / 0.2e1;
t166 = t211 + t209;
t208 = -t98 / 0.2e1;
t188 = t128 * mrSges(4,2);
t101 = t126 * mrSges(4,1) + t188;
t206 = -t101 / 0.2e1;
t205 = -t124 / 0.2e1;
t204 = t124 / 0.2e1;
t203 = t125 / 0.2e1;
t202 = t126 / 0.2e1;
t201 = m(6) * t124;
t200 = pkin(6) * t129;
t199 = Ifges(5,4) * t124;
t198 = Ifges(5,4) * t125;
t197 = Ifges(6,5) * t124;
t196 = Ifges(6,5) * t125;
t97 = -pkin(3) * t128 - qJ(4) * t126 - pkin(2);
t191 = t125 * t97;
t190 = t126 * mrSges(6,1);
t189 = t126 * mrSges(4,2);
t173 = t127 * t128;
t76 = -t124 * t129 + t125 * t173;
t187 = t128 * t76;
t122 = t126 ^ 2;
t172 = t127 * t129;
t106 = t122 * t172;
t123 = t128 ^ 2;
t75 = t124 * t173 + t125 * t129;
t8 = m(4) * (t106 + (t123 - 0.1e1) * t172) + t222 * (t75 * t77 + t76 * t78 + t106);
t184 = t8 * qJD(1);
t177 = t125 * t128;
t57 = pkin(6) * t177 + t124 * t97;
t148 = t124 * t75 + t125 * t76;
t176 = t126 * t127;
t11 = t222 * (t128 * t127 ^ 2 * t126 - t148 * t176);
t181 = t11 * qJD(1);
t179 = t124 * t128;
t169 = qJ(4) * qJD(3);
t168 = qJD(3) * t124;
t167 = qJD(3) * t128;
t165 = m(6) * t178;
t164 = qJD(5) * t201;
t163 = t122 * t174;
t99 = -mrSges(5,1) * t125 + mrSges(5,2) * t124;
t157 = -m(5) * pkin(3) - mrSges(4,1) + t99;
t156 = t166 * t129;
t154 = t124 * mrSges(5,1) + t125 * mrSges(5,2);
t153 = t124 * mrSges(6,1) - t125 * mrSges(6,3);
t39 = -qJ(5) * t128 + t57;
t40 = t128 * t159 - t191;
t56 = -pkin(6) * t179 + t191;
t62 = t126 * Ifges(6,6) + (t124 * Ifges(6,3) + t196) * t128;
t63 = t126 * Ifges(5,6) + (-t124 * Ifges(5,2) + t198) * t128;
t64 = t126 * Ifges(6,4) + (t125 * Ifges(6,1) + t197) * t128;
t65 = t126 * Ifges(5,5) + (t125 * Ifges(5,1) - t199) * t128;
t143 = pkin(4) * t124 - qJ(5) * t125 + pkin(6);
t66 = t143 * t126;
t67 = t143 * t128;
t79 = t153 * t126;
t80 = t154 * t126;
t81 = t153 * t128;
t82 = t154 * t128;
t88 = -mrSges(5,2) * t126 - mrSges(5,3) * t179;
t91 = mrSges(5,1) * t126 - mrSges(5,3) * t177;
t92 = mrSges(6,2) * t177 - t190;
t93 = -mrSges(6,2) * t179 + mrSges(6,3) * t126;
t1 = -pkin(2) * t101 + t39 * t93 + t40 * t92 + t41 * t94 + t42 * t90 + t56 * t91 + t57 * t88 + t58 * t89 + t59 * t87 + t66 * t81 + t67 * t79 + m(5) * (t56 * t58 + t57 * t59) + m(6) * (t39 * t41 + t40 * t42 + t66 * t67) + (-Ifges(4,4) * t126 + pkin(6) * t82 + (t64 / 0.2e1 + t65 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t126) * t125 + (t62 / 0.2e1 - t63 / 0.2e1 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t126) * t124) * t126 + (pkin(6) * t80 + (m(5) * pkin(6) ^ 2 - Ifges(6,2) - Ifges(5,3) + Ifges(4,1) - Ifges(4,2) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t121 + ((Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t124 + (-Ifges(5,4) + Ifges(6,5)) * t125) * t124) * t126 + (t219 * t124 - t220 * t125 + Ifges(4,4)) * t128) * t128;
t136 = (-t87 / 0.2e1 - t94 / 0.2e1) * t125 + (t89 / 0.2e1 - t90 / 0.2e1) * t124;
t150 = -t124 * t56 + t125 * t57;
t151 = t124 * t40 + t125 * t39;
t130 = (t88 / 0.2e1 + t93 / 0.2e1) * t76 + (-t91 / 0.2e1 + t92 / 0.2e1) * t75 + ((t66 * t209 + t80 / 0.2e1 + t79 / 0.2e1) * t128 + (t82 / 0.2e1 + t81 / 0.2e1 + (0.2e1 * pkin(6) * t128 - t150) * t211 + (-t151 + t67) * t209 + t136) * t126) * t127 + (-t58 * t75 + t59 * t76) * t211 + (t41 * t76 + t42 * t75) * t209;
t3 = (t188 / 0.2e1 + t206 + (mrSges(4,1) / 0.2e1 - t99 / 0.2e1 + t208) * t126) * t129 + t130 + t224;
t152 = t3 * qJD(1) + t1 * qJD(2);
t149 = -t124 * t76 + t125 * t75;
t13 = (t213 * t149 + t156) * t126;
t7 = (t217 * t125 - t218 * t124 + m(5) * (-t124 * t57 - t125 * t56) + m(6) * (-t124 * t39 + t125 * t40)) * t126;
t147 = -qJD(1) * t13 + qJD(2) * t7;
t146 = t166 * t173;
t17 = t128 * t94 - m(6) * (-t128 * t39 - t178 * t66) + t79 * t178;
t23 = 0.2e1 * (t77 / 0.4e1 + t163 / 0.4e1 + t187 / 0.4e1) * m(6);
t145 = -qJD(1) * t23 - qJD(2) * t17;
t133 = (-t124 * t66 + (-t126 * t95 - t183) * t125) * t209 + t79 * t205;
t142 = t42 * t209 - t190 / 0.2e1;
t10 = (t128 * mrSges(6,2) + t202 * t98) * t125 - t133 + t142;
t37 = t160 * t124;
t144 = -qJD(2) * t10 - qJD(3) * t37;
t134 = t166 * t148;
t15 = t146 - t134;
t28 = (0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * qJ(4) + t221) * t170;
t132 = (pkin(6) * t211 + (-mrSges(6,3) / 0.2e1 + mrSges(5,2) / 0.2e1) * t125 + (mrSges(6,1) / 0.2e1 + mrSges(5,1) / 0.2e1) * t124) * t128 + t67 * t209;
t135 = t150 * t212 + t151 * t210;
t4 = t132 + t135 + t136;
t138 = -qJD(1) * t15 - qJD(2) * t4 + qJD(3) * t28;
t115 = t122 * t200;
t83 = (qJD(2) * t178 + t168) * m(6);
t24 = (-t163 - t187 + t77) * t209;
t16 = t146 + t134;
t14 = t222 * t149 * t202 + t126 * t156;
t12 = t178 * t208 + t133 + t142;
t5 = t218 * t203 + t204 * t90 + t205 * t89 + t132 - t135;
t2 = -mrSges(4,1) * t175 / 0.2e1 - mrSges(4,2) * t171 / 0.2e1 + t129 * t206 + t130 + (t99 + t98) * t175 / 0.2e1 - t224;
t6 = [qJD(2) * t8 + qJD(3) * t11, t2 * qJD(3) + t14 * qJD(4) + t24 * qJD(5) + t184 + (t218 * t78 + t217 * t77 + (-t128 * mrSges(4,1) - mrSges(3,1) + t189) * t127 + (-mrSges(3,2) + (t79 + t80) * t126 + (t122 + t123) * mrSges(4,3)) * t129 + m(4) * (-t127 * pkin(2) + t123 * t200 + t115) + 0.2e1 * (-t56 * t77 + t57 * t78 + t115) * t211 + 0.2e1 * (t175 * t66 + t39 * t78 + t40 * t77) * t209) * qJD(2), t181 + t2 * qJD(2) + t16 * qJD(4) + ((t157 + t160) * t167 + (mrSges(4,2) * qJD(3) - t164 + (-t221 * qJD(3) + 0.2e1 * t213 * t169) * t170) * t126) * t127, qJD(2) * t14 + qJD(3) * t16, -m(6) * t168 * t176 + qJD(2) * t24; qJD(3) * t3 - qJD(4) * t13 - qJD(5) * t23 - t184, qJD(3) * t1 + qJD(4) * t7 - qJD(5) * t17, (-Ifges(4,6) * t126 - t125 * t62 / 0.2e1 + t63 * t203 + t95 * t81 - pkin(3) * t82 + pkin(6) * t189 + t160 * t67 + (t64 + t65) * t204 + (t220 * t124 + t219 * t125) * t202 + t215 * mrSges(5,3) + t214 * mrSges(6,2)) * qJD(3) + t5 * qJD(4) + t12 * qJD(5) + ((t88 + t93) * t125 + (-t91 + t92) * t124 + m(5) * t215 + m(6) * t214) * t169 + ((-Ifges(6,3) * t125 + t197) * t204 + (Ifges(5,2) * t125 + t199) * t205 + Ifges(4,5) + t157 * pkin(6) + (-t196 + t198 + (Ifges(5,1) + Ifges(6,1)) * t124) * t203) * t167 + t152, qJD(3) * t5 + t147, qJD(3) * t12 + t145; -qJD(2) * t3 - qJD(4) * t15 - t181, -qJD(4) * t4 - qJD(5) * t10 - t152, qJD(4) * t28 - qJD(5) * t37, t138, t144; qJD(2) * t13 + qJD(3) * t15, qJD(3) * t4 - qJD(5) * t165 - t147, -t138 - t164, 0, -t83; t23 * qJD(2), qJD(3) * t10 + qJD(4) * t165 - t145, qJD(4) * t201 - t144, t83, 0;];
Cq = t6;
