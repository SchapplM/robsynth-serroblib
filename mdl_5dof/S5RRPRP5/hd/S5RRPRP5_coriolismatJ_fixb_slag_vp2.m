% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:04
% EndTime: 2019-12-31 19:54:10
% DurationCPUTime: 2.72s
% Computational Cost: add. (5942->220), mult. (11561->286), div. (0->0), fcn. (12656->6), ass. (0->114)
t124 = sin(pkin(8));
t166 = cos(pkin(8));
t194 = sin(qJ(2));
t196 = cos(qJ(2));
t114 = -t124 * t194 + t166 * t196;
t125 = sin(qJ(4));
t137 = t124 * t196 + t166 * t194;
t195 = cos(qJ(4));
t212 = t125 * t114 + t195 * t137;
t120 = -t196 * pkin(2) - pkin(1);
t104 = -t114 * pkin(3) + t120;
t173 = qJ(5) * t212;
t93 = -t195 * t114 + t125 * t137;
t35 = t93 * pkin(4) + t104 - t173;
t221 = m(6) * t35 + mrSges(6,1) * t93 - mrSges(6,3) * t212;
t247 = Ifges(5,1) + Ifges(6,1);
t90 = Ifges(6,5) * t212;
t250 = t90 + (Ifges(6,3) - t247) * t93;
t160 = t194 * pkin(6);
t138 = -t194 * qJ(3) - t160;
t161 = t196 * pkin(6);
t139 = t196 * qJ(3) + t161;
t213 = -t124 * t139 + t166 * t138;
t220 = -pkin(7) * t137 + t213;
t99 = t124 * t138 + t166 * t139;
t63 = t114 * pkin(7) + t99;
t38 = t125 * t63 - t195 * t220;
t249 = t125 * t220 + t195 * t63;
t222 = t212 * mrSges(6,1);
t152 = t93 * mrSges(6,3) + t222;
t238 = t93 * mrSges(5,2);
t153 = t212 * mrSges(5,1) - t238;
t190 = Ifges(6,5) * t93;
t193 = Ifges(5,4) * t212;
t203 = t93 / 0.2e1;
t230 = -t212 / 0.2e1;
t241 = -t93 / 0.2e1;
t248 = (-Ifges(5,2) * t93 + t193) * t230 + (Ifges(6,3) * t212 - t190) * t203 + t104 * t153 + t35 * t152 + (t190 - 0.2e1 * Ifges(5,4) * t93 + (-Ifges(5,2) + t247) * t212) * t241;
t167 = qJ(5) * t93;
t198 = pkin(4) * t212;
t46 = t167 + t198;
t123 = t194 * pkin(2);
t239 = m(4) * t123;
t237 = t93 * mrSges(6,2);
t229 = t212 / 0.2e1;
t231 = t230 + t229;
t217 = mrSges(6,1) + mrSges(5,1);
t215 = Ifges(5,5) + Ifges(6,4);
t214 = -Ifges(5,6) + Ifges(6,6);
t121 = m(6) * qJ(5) + mrSges(6,3);
t211 = 0.2e1 * m(6);
t210 = m(5) / 0.2e1;
t209 = m(6) / 0.2e1;
t208 = m(4) * pkin(2);
t207 = t249 / 0.2e1;
t199 = m(6) * t249;
t187 = pkin(2) * t124;
t148 = t195 * t187;
t154 = t166 * pkin(2);
t119 = t154 + pkin(3);
t163 = t125 * t119;
t109 = t148 + t163;
t197 = t109 / 0.2e1;
t192 = Ifges(6,4) * t93;
t191 = Ifges(5,5) * t93;
t189 = Ifges(5,6) * t212;
t188 = Ifges(6,6) * t212;
t185 = t249 * mrSges(5,1);
t184 = t249 * mrSges(6,1);
t183 = t38 * mrSges(5,2);
t182 = t38 * mrSges(6,3);
t128 = t137 ^ 2;
t5 = (t114 ^ 2 + t128) * mrSges(4,3) + m(4) * (t99 * t114 - t137 * t213) + (m(6) + m(5)) * (t212 * t38 - t249 * t93) + (mrSges(5,3) + mrSges(6,2)) * (t212 ^ 2 + t93 ^ 2);
t172 = qJD(1) * t5;
t105 = pkin(3) * t137 + t123;
t131 = mrSges(4,1) * t137 + t114 * mrSges(4,2);
t36 = t105 + t46;
t1 = -Ifges(4,4) * t128 - pkin(1) * (t194 * mrSges(3,1) + t196 * mrSges(3,2)) + (-mrSges(4,1) * t123 + Ifges(4,4) * t114) * t114 + (-Ifges(3,2) + Ifges(3,1)) * t196 * t194 + (mrSges(4,2) * t123 + (Ifges(4,1) - Ifges(4,2)) * t114) * t137 + (-t194 ^ 2 + t196 ^ 2) * Ifges(3,4) + ((-Ifges(5,4) + Ifges(6,5)) * t212 + t250) * t229 + (t131 + t239) * t120 + t248 + t221 * t36 + (m(5) * t104 + mrSges(5,1) * t93 + mrSges(5,2) * t212) * t105;
t171 = t1 * qJD(1);
t2 = (-t193 + t90 + t250) * t229 + t221 * t46 + t248;
t170 = t2 * qJD(1);
t106 = qJ(5) + t109;
t108 = t195 * t119 - t125 * t187;
t107 = -pkin(4) - t108;
t126 = (-t108 * t212 - t109 * t93) * t210 + (-t106 * t93 + t107 * t212) * t209 + (t124 * t114 - t166 * t137) * t208 / 0.2e1;
t132 = t105 * t210 + t36 * t209 + t239 / 0.2e1;
t6 = t126 - t131 - t132 - t152 - t153;
t169 = t6 * qJD(1);
t7 = 0.2e1 * t241 * mrSges(6,3) + 0.2e1 * t230 * mrSges(5,1) + (-t198 / 0.4e1 - t167 / 0.4e1 - t46 / 0.4e1) * t211 - t222 + t238;
t168 = t7 * qJD(1);
t14 = t221 * t212;
t165 = qJD(1) * t14;
t41 = t211 * t229;
t162 = t41 * qJD(1);
t156 = t241 + t203;
t135 = -(-mrSges(6,3) + mrSges(5,2)) * t108 - t217 * t109;
t15 = m(6) * (t106 * t108 + t107 * t109) + t135;
t129 = ((-t106 + t109) * t38 + (t107 + t108) * t249) * t209;
t133 = (-t106 / 0.2e1 + t197) * t212 + (-t107 / 0.2e1 - t108 / 0.2e1) * t93;
t140 = m(6) * (-pkin(4) * t249 - qJ(5) * t38);
t3 = -t140 / 0.2e1 + (t173 / 0.2e1 + pkin(4) * t241 + t133) * mrSges(6,2) + t129 + t217 * (-t249 / 0.2e1 + t207) + t215 * t156 + t214 * t231;
t144 = t3 * qJD(1) + t15 * qJD(2);
t101 = m(6) * t106 + mrSges(6,3);
t13 = t156 * mrSges(6,2);
t143 = qJD(1) * t13 + qJD(2) * t101;
t59 = mrSges(6,3) + (t148 / 0.4e1 + t163 / 0.4e1 + qJ(5) / 0.2e1 - t109 / 0.4e1) * t211;
t142 = qJD(2) * t59 + qJD(4) * t121;
t58 = (0.2e1 * qJ(5) + t109) * t209 + mrSges(6,3) + m(6) * t197;
t42 = t231 * m(6);
t19 = -t237 + t199;
t12 = -t237 + t199 / 0.2e1 + m(6) * t207;
t11 = t126 + t132;
t4 = -t192 / 0.2e1 - t191 / 0.2e1 - t189 / 0.2e1 + t188 / 0.2e1 + t140 / 0.2e1 - t184 / 0.2e1 - t185 / 0.2e1 - t182 / 0.2e1 + t183 / 0.2e1 + pkin(4) * t237 / 0.2e1 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t212 + (-Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t93 + (-mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1) * t249 + (mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1) * t38 + t129 + (-t173 / 0.2e1 + t133) * mrSges(6,2);
t8 = [qJD(2) * t1 + qJD(3) * t5 + qJD(4) * t2 - qJD(5) * t14, t171 + (-Ifges(4,6) * t137 + Ifges(4,5) * t114 - t99 * mrSges(4,1) - t213 * mrSges(4,2) - t191 - t192 + t188 - t189 - t185 - t184 - t182 + t183 + mrSges(3,2) * t160 - t107 * t237 + m(6) * (-t106 * t38 + t107 * t249) + m(5) * (-t108 * t249 - t109 * t38) + (t124 * t213 - t166 * t99) * t208 - mrSges(3,1) * t161 - Ifges(3,6) * t194 + Ifges(3,5) * t196 - t106 * t212 * mrSges(6,2) + (t108 * t93 - t109 * t212) * mrSges(5,3) + (-t114 * t154 - t137 * t187) * mrSges(4,3)) * qJD(2) + t11 * qJD(3) + t4 * qJD(4) + t12 * qJD(5), qJD(2) * t11 + qJD(5) * t42 + t172, t4 * qJD(2) + t19 * qJD(5) + t170 + ((-m(6) * pkin(4) - t217) * t249 + (mrSges(5,2) - t121) * t38 + (-qJ(5) * mrSges(6,2) + t214) * t212 + (pkin(4) * mrSges(6,2) - t215) * t93) * qJD(4), qJD(2) * t12 + qJD(3) * t42 + qJD(4) * t19 - t165; qJD(3) * t6 + qJD(4) * t3 + qJD(5) * t13 - t171, qJD(4) * t15 + qJD(5) * t101, t169, (m(6) * (-pkin(4) * t109 + qJ(5) * t108) + t135) * qJD(4) + t58 * qJD(5) + t144, qJD(4) * t58 + t143; -qJD(2) * t6 - qJD(4) * t7 - qJD(5) * t41 - t172, -t169, 0, -t168, -t162; -qJD(2) * t3 + qJD(3) * t7 - t170, qJD(5) * t59 - t144, t168, t121 * qJD(5), t142; -qJD(2) * t13 + qJD(3) * t41 + t165, -qJD(4) * t59 - t143, t162, -t142, 0;];
Cq = t8;
