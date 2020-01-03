% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:33
% EndTime: 2019-12-31 17:59:35
% DurationCPUTime: 1.21s
% Computational Cost: add. (1964->202), mult. (4092->300), div. (0->0), fcn. (3100->6), ass. (0->123)
t125 = sin(qJ(5));
t122 = t125 ^ 2;
t127 = cos(qJ(5));
t124 = t127 ^ 2;
t215 = t122 + t124;
t221 = t215 - 0.1e1;
t220 = mrSges(6,3) * t215;
t120 = Ifges(6,5) * t127;
t216 = -t125 * Ifges(6,6) + t120;
t219 = -Ifges(5,4) + t216;
t218 = qJD(4) * (-mrSges(5,2) + t220);
t189 = t127 * mrSges(6,1);
t191 = t125 * mrSges(6,2);
t158 = t189 - t191;
t217 = -mrSges(5,1) - t158;
t126 = sin(qJ(4));
t128 = cos(qJ(4));
t175 = t126 * t128;
t121 = Ifges(6,4) * t127;
t214 = Ifges(6,2) * t125 - t121;
t123 = t126 ^ 2;
t212 = t128 ^ 2;
t211 = -mrSges(6,1) / 0.2e1;
t32 = (-t123 + t212) * t221;
t210 = -t32 / 0.2e1;
t188 = t127 * mrSges(6,2);
t192 = t125 * mrSges(6,1);
t99 = t188 + t192;
t209 = -t99 / 0.2e1;
t208 = pkin(7) * mrSges(6,3);
t113 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t207 = t113 / 0.2e1;
t206 = t125 / 0.2e1;
t205 = -t126 / 0.2e1;
t204 = t126 / 0.2e1;
t203 = -t127 / 0.2e1;
t202 = -t128 / 0.2e1;
t201 = t128 / 0.2e1;
t200 = m(6) * t113;
t199 = pkin(4) * t128;
t198 = pkin(7) * t128;
t197 = t126 * pkin(4);
t196 = t126 * pkin(7);
t195 = m(6) * qJD(4);
t194 = Ifges(6,4) * t125;
t114 = sin(pkin(8)) * pkin(1) + qJ(3);
t83 = t114 + t197 - t198;
t187 = t127 * t83;
t186 = t128 * mrSges(5,2);
t176 = t125 * t128;
t150 = t126 * mrSges(6,2) + mrSges(6,3) * t176;
t134 = t150 * t203;
t172 = t127 * t128;
t152 = t126 * mrSges(6,1) - mrSges(6,3) * t172;
t135 = t152 * t205;
t177 = t125 * t126;
t146 = t113 * t177 - t187;
t136 = t146 * t125;
t149 = -mrSges(6,2) * t128 + mrSges(6,3) * t177;
t139 = t127 * t149;
t173 = t127 * t126;
t151 = mrSges(6,1) * t128 + mrSges(6,3) * t173;
t141 = t125 * t151;
t182 = t83 * t125;
t31 = t113 * t173 + t182;
t6 = t126 * t134 + t139 * t202 + t125 * t135 + t141 * t201 + t212 * t209 + t123 * t99 / 0.2e1 - m(6) * (-t126 * t136 - t31 * t173 + (t221 * t212 + t123) * t113) / 0.2e1;
t184 = t6 * qJD(1);
t138 = t99 * t175;
t143 = t126 * t99;
t178 = t113 * t126;
t8 = m(6) * ((t31 - t182) * t127 + (t124 - 0.2e1 + 0.2e1 * t122) * t178) * t201 + t128 * t134 + t139 * t204 - t152 * t176 / 0.2e1 + t141 * t205 + t138 / 0.2e1 - t143 * t202;
t183 = t8 * qJD(1);
t140 = t125 * t150;
t84 = t158 * t128;
t129 = t127 * t135 - t140 * t205 + t84 * t202 - t175 * t220 / 0.2e1;
t147 = t191 / 0.2e1 - t189 / 0.2e1;
t9 = t129 + t147;
t181 = t9 * qJD(1);
t13 = -t147 * t175 + t84 * t205;
t180 = qJD(1) * t13;
t14 = -t140 + t127 * t152 + t186 + t126 * mrSges(5,1) + mrSges(4,3) + m(6) * (t31 * t125 - t146 * t127) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t114;
t179 = qJD(1) * t14;
t104 = t196 + t199;
t174 = t127 * t104;
t160 = t215 * t128;
t46 = t126 * t160 - t175;
t169 = t46 * t195;
t168 = -pkin(4) * t84 / 0.2e1;
t167 = t208 / 0.2e1;
t166 = Ifges(6,4) * t172;
t165 = t177 / 0.2e1;
t164 = -t173 / 0.2e1;
t161 = qJD(4) * t217;
t157 = Ifges(6,1) * t127 - t194;
t156 = Ifges(6,5) * t125 + Ifges(6,6) * t127;
t4 = -t146 * t150 + t31 * t152 + (t136 * mrSges(6,3) + t113 * t84 + (-Ifges(6,4) * t176 + Ifges(6,5) * t126) * t125 + (t166 + t31 * mrSges(6,3) + Ifges(6,6) * t126 + (Ifges(6,1) - Ifges(6,2)) * t176) * t127) * t128;
t155 = -t4 * qJD(1) - t13 * qJD(2);
t154 = qJD(2) * t210 - t46 * qJD(3);
t153 = t46 * qJD(2) + qJD(3) * t210;
t148 = t192 / 0.2e1 + t188 / 0.2e1;
t145 = t125 * t104 + t113 * t172;
t144 = t113 * t176 - t174;
t1 = -t113 * t138 + t145 * t150 - t31 * t149 + t144 * t152 + t146 * t151 - m(6) * (-t144 * t187 + t31 * t145) - t114 * (t128 * mrSges(5,1) - t126 * mrSges(5,2)) + (t219 * t126 + (t124 * Ifges(6,1) + t122 * Ifges(6,2) + Ifges(5,1) - Ifges(5,2) - Ifges(6,3) + (-(-0.1e1 + t122) * t200 - t99) * t113) * t128 + (t174 * t200 - 0.2e1 * t166) * t125) * t126 - t219 * t212;
t142 = -t1 * qJD(1) - t6 * qJD(2) + t8 * qJD(3);
t137 = t209 + t148;
t133 = t113 * t209 + (-t124 / 0.2e1 - t122 / 0.2e1) * t208;
t100 = Ifges(6,2) * t127 + t194;
t132 = (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.4e1 + t167) * t127 - t100 / 0.4e1;
t101 = Ifges(6,1) * t125 + t121;
t130 = -t101 / 0.4e1 - t121 / 0.4e1 + (Ifges(6,2) / 0.2e1 + t167 - Ifges(6,1) / 0.4e1) * t125;
t2 = t168 + t147 * t104 + (t147 * pkin(7) + t216) * t126 + (-Ifges(6,3) / 0.2e1 + (mrSges(6,2) * t207 + t132) * t127 + (-0.5e1 / 0.4e1 * t121 + mrSges(6,1) * t207 + t130) * t125 + t133) * t128;
t20 = pkin(4) * t99 - t125 * t157 / 0.2e1 + t100 * t206 + (t101 - t214) * t203;
t33 = t137 * t126;
t34 = t137 * t128;
t131 = t2 * qJD(1) - t33 * qJD(2) + t34 * qJD(3) - t20 * qJD(4);
t36 = -t148 * t128 + t99 * t202;
t35 = t148 * t126 + t99 * t204;
t30 = t32 * t195 / 0.2e1;
t10 = t129 - t147;
t7 = t8 * qJD(4);
t5 = -qJD(4) * t6 - qJD(5) * t13;
t3 = t126 * t120 / 0.4e1 + t168 + Ifges(6,5) * t164 + Ifges(6,6) * t165 + Ifges(6,3) * t201 - t145 * mrSges(6,2) / 0.2e1 + t144 * t211 + t133 * t128 + ((-Ifges(6,6) / 0.2e1 + pkin(7) * mrSges(6,2) / 0.2e1) * t126 + t130 * t128) * t125 + ((pkin(7) * t211 + Ifges(6,5) / 0.4e1) * t126 + (-0.5e1 / 0.4e1 * t194 + t132) * t128) * t127;
t11 = [qJD(3) * t14 - qJD(4) * t1 - qJD(5) * t4, t5, qJD(5) * t10 + t179 + t7, (t101 * t164 + t100 * t165 - t113 * t186 + t127 * (Ifges(6,6) * t128 + t214 * t126) / 0.2e1 + t156 * t201 + (t198 * t215 - t197) * t200 + (Ifges(6,5) * t128 - t157 * t126) * t206 + pkin(4) * t143 - Ifges(5,5) * t126 - Ifges(5,6) * t128 + t217 * t178 + (t139 - t141) * pkin(7) + (t144 * t125 + t145 * t127) * mrSges(6,3)) * qJD(4) + t3 * qJD(5) + t142, t10 * qJD(3) + t3 * qJD(4) + (-t31 * mrSges(6,1) + t146 * mrSges(6,2) - t128 * t156) * qJD(5) + t155; t5, -t169, t30, -t184 + t35 * qJD(5) + t128 * t161 - t126 * t218 + ((-t196 * t215 - t199) * qJD(4) - t153) * m(6), qJD(4) * t35 - qJD(5) * t84 - t180; qJD(5) * t9 - t179 + t7, t30, t169, t183 + t36 * qJD(5) + t126 * t161 + t128 * t218 + ((pkin(7) * t160 - t197) * qJD(4) - t154) * m(6), -qJD(5) * t126 * t158 + t36 * qJD(4) + t181; qJD(5) * t2 - t142, t153 * m(6) - t33 * qJD(5) + t184, t154 * m(6) + t34 * qJD(5) - t183, -t20 * qJD(5), (-pkin(7) * t158 + t216) * qJD(5) + t131; -qJD(3) * t9 - qJD(4) * t2 - t155, qJD(4) * t33 + t180, -qJD(4) * t34 - t181, -t131, 0;];
Cq = t11;
