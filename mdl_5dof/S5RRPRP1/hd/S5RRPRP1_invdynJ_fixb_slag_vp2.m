% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:58:59
% DurationCPUTime: 2.62s
% Computational Cost: add. (2026->337), mult. (3250->414), div. (0->0), fcn. (1649->12), ass. (0->152)
t232 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t151 = cos(qJ(4));
t230 = -mrSges(6,1) - mrSges(5,1);
t231 = t151 * t230 - mrSges(4,1);
t229 = Ifges(5,1) + Ifges(6,1);
t227 = Ifges(6,5) + Ifges(5,5);
t226 = Ifges(5,2) + Ifges(6,2);
t225 = Ifges(6,6) + Ifges(5,6);
t142 = qJDD(1) + qJDD(2);
t145 = sin(pkin(8));
t146 = cos(pkin(8));
t149 = sin(qJ(2));
t200 = qJD(1) * pkin(1);
t181 = t149 * t200;
t152 = cos(qJ(2));
t214 = pkin(1) * t152;
t85 = -qJD(2) * t181 + qJDD(1) * t214;
t63 = pkin(2) * t142 + t85;
t187 = qJD(2) * t152;
t86 = (qJD(1) * t187 + qJDD(1) * t149) * pkin(1);
t26 = t145 * t63 + t146 * t86;
t18 = pkin(7) * t142 + t26;
t224 = qJD(3) * qJD(4) + t18;
t223 = (Ifges(5,4) + Ifges(6,4)) * t151;
t222 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t207 = mrSges(5,2) + mrSges(6,2);
t136 = t151 * qJD(3);
t148 = sin(qJ(4));
t143 = qJD(1) + qJD(2);
t180 = t152 * t200;
t95 = pkin(2) * t143 + t180;
t45 = t145 * t95 + t146 * t181;
t33 = pkin(7) * t143 + t45;
t173 = qJ(5) * t143 + t33;
t15 = -t173 * t148 + t136;
t13 = qJD(4) * pkin(4) + t15;
t28 = -t148 * t33 + t136;
t220 = -t28 * mrSges(5,3) - t13 * mrSges(6,3);
t144 = qJ(1) + qJ(2);
t137 = sin(t144);
t125 = pkin(2) * t137;
t138 = cos(t144);
t126 = pkin(2) * t138;
t213 = pkin(2) * t145;
t212 = pkin(2) * t146;
t211 = pkin(4) * t151;
t185 = qJD(4) * t148;
t5 = t148 * qJDD(3) + t224 * t151 - t33 * t185;
t209 = t151 * t5;
t193 = t143 * t148;
t91 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t193;
t92 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t193;
t206 = t91 + t92;
t192 = t143 * t151;
t93 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t192;
t94 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t192;
t205 = t93 + t94;
t204 = Ifges(5,4) * t148;
t202 = Ifges(6,4) * t148;
t128 = pkin(2) + t214;
t190 = t146 * t149;
t72 = pkin(1) * t190 + t145 * t128;
t65 = pkin(7) + t72;
t198 = -qJ(5) - t65;
t186 = qJD(3) * t148;
t29 = t151 * t33 + t186;
t197 = qJD(4) * t29;
t196 = qJD(4) * t65;
t191 = t145 * t149;
t122 = pkin(7) + t213;
t189 = -qJ(5) - t122;
t150 = sin(qJ(1));
t140 = t150 * pkin(1);
t188 = t125 + t140;
t184 = qJD(4) * t151;
t135 = t151 * qJD(5);
t179 = pkin(4) * t185;
t134 = pkin(8) + t144;
t120 = sin(t134);
t121 = cos(t134);
t127 = pkin(3) + t211;
t147 = -qJ(5) - pkin(7);
t178 = t120 * t127 + t121 * t147 + t125;
t177 = t121 * pkin(3) + t120 * pkin(7) + t126;
t176 = t143 * t185;
t77 = t142 * t151 - t176;
t78 = t142 * t148 + t143 * t184;
t30 = -t77 * mrSges(6,1) + t78 * mrSges(6,2);
t25 = -t145 * t86 + t146 * t63;
t104 = t145 * t181;
t44 = t146 * t95 - t104;
t172 = qJD(4) * t198;
t71 = -pkin(1) * t191 + t128 * t146;
t171 = qJD(4) * t189;
t64 = -pkin(3) - t71;
t167 = -t120 * t147 + t121 * t127 + t126;
t16 = t173 * t151 + t186;
t165 = t29 * mrSges(5,3) + t16 * mrSges(6,3);
t164 = g(2) * t121 + g(3) * t120;
t102 = -mrSges(5,1) * t151 + mrSges(5,2) * t148;
t163 = -mrSges(6,1) * t151 + mrSges(6,2) * t148;
t162 = Ifges(5,2) * t151 + t204;
t161 = Ifges(6,2) * t151 + t202;
t160 = -t13 * t148 + t151 * t16;
t159 = -t148 * t28 + t151 * t29;
t17 = -pkin(3) * t142 - t25;
t158 = pkin(1) * (t145 * t152 + t190);
t67 = qJD(2) * t158;
t157 = -t138 * mrSges(3,1) + t137 * mrSges(3,2) + t222 * t120 + t231 * t121;
t156 = -t137 * mrSges(3,1) - t138 * mrSges(3,2) + (m(5) * pkin(7) - t222) * t121 + t231 * t120;
t155 = t207 * t164;
t27 = -t127 * t143 + qJD(5) - t44;
t3 = qJ(5) * t77 + t143 * t135 + t5;
t32 = -pkin(3) * t143 - t44;
t59 = Ifges(6,6) * qJD(4) + t161 * t143;
t60 = Ifges(5,6) * qJD(4) + t162 * t143;
t105 = Ifges(6,4) * t192;
t61 = Ifges(6,1) * t193 + Ifges(6,5) * qJD(4) + t105;
t106 = Ifges(5,4) * t192;
t62 = Ifges(5,1) * t193 + Ifges(5,5) * qJD(4) + t106;
t8 = -pkin(4) * t77 + qJDD(5) + t17;
t154 = t3 * t151 * mrSges(6,3) + t85 * mrSges(3,1) + t25 * mrSges(4,1) - t86 * mrSges(3,2) + mrSges(5,3) * t209 + t17 * t102 + t8 * t163 + (-t225 * t148 + t227 * t151) * qJD(4) ^ 2 / 0.2e1 - (t60 + t59) * t185 / 0.2e1 + (t229 * t151 - t202 - t204) * t176 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t142 + (t32 * (mrSges(5,1) * t148 + mrSges(5,2) * t151) + t27 * (mrSges(6,1) * t148 + mrSges(6,2) * t151)) * qJD(4) + (t162 / 0.2e1 + t161 / 0.2e1 + t148 * t232 + t226 * t151 / 0.2e1) * t77 + (t229 * t148 + t223 / 0.2e1 + t151 * t232) * t78 + (t62 + t61 + (-t226 * t148 + t223) * t143) * t184 / 0.2e1 + (t227 * t148 + t225 * t151) * qJDD(4);
t153 = cos(qJ(1));
t141 = t153 * pkin(1);
t139 = t151 * qJ(5);
t133 = t151 * qJDD(3);
t123 = -pkin(3) - t212;
t113 = t120 * pkin(3);
t100 = -t127 - t212;
t84 = t122 * t151 + t139;
t83 = t189 * t148;
t76 = t102 * t143;
t75 = t163 * t143;
t69 = (t146 * t152 - t191) * qJD(2) * pkin(1);
t68 = t146 * t180 - t104;
t66 = qJD(1) * t158;
t54 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t78;
t53 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t78;
t52 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t77;
t51 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t77;
t50 = -qJD(5) * t148 + t151 * t171;
t49 = t148 * t171 + t135;
t48 = t64 - t211;
t46 = t67 + t179;
t39 = t151 * t65 + t139;
t38 = t198 * t148;
t31 = -mrSges(5,1) * t77 + mrSges(5,2) * t78;
t11 = (-qJD(5) - t69) * t148 + t151 * t172;
t10 = t148 * t172 + t151 * t69 + t135;
t6 = -t148 * t18 + t133 - t197;
t2 = -t33 * t184 + qJDD(4) * pkin(4) - qJ(5) * t78 + t133 + (-qJD(5) * t143 - t224) * t148;
t1 = [(-t150 * mrSges(2,1) - mrSges(2,2) * t153 - m(4) * t188 - m(6) * (t140 + t178) - m(5) * (t113 + t188) + t156) * g(3) + (-t142 * t72 - t143 * t69 - t26) * mrSges(4,2) + m(5) * (t17 * t64 + t32 * t67) + t11 * t91 + t10 * t93 + t46 * t75 + t67 * t76 + t39 * t51 + t38 * t53 + t64 * t31 + (-mrSges(2,1) * t153 + t150 * mrSges(2,2) - m(4) * (t126 + t141) - m(6) * (t141 + t167) - m(5) * (t141 + t177) + t157) * g(2) + t154 + t48 * t30 + (t142 * t71 - t143 * t67) * mrSges(4,1) + Ifges(2,3) * qJDD(1) + (m(5) * (-t29 * t196 - t28 * t69 - t6 * t65) - t65 * t54 - t69 * t92 - t94 * t196 + (-qJD(4) * t16 - t2) * mrSges(6,3) + (-t6 - t197) * mrSges(5,3) + t155) * t148 + (m(5) * (-t28 * t196 + t29 * t69 + t5 * t65) + t65 * t52 + t69 * t94 - t92 * t196 + t220 * qJD(4)) * t151 + ((-t142 * t149 - t143 * t187) * mrSges(3,2) + (-qJD(2) * t143 * t149 + t142 * t152) * mrSges(3,1) + (-g(2) * t153 - g(3) * t150 + t149 * t86 + t152 * t85) * m(3)) * pkin(1) + m(6) * (t10 * t16 + t11 * t13 + t2 * t38 + t27 * t46 + t3 * t39 + t48 * t8) + m(4) * (t25 * t71 + t26 * t72 - t44 * t67 + t45 * t69); (-t6 * mrSges(5,3) - t2 * mrSges(6,3) - t122 * t54 + t206 * t68 + (pkin(4) * t75 - t122 * t94 - t165) * qJD(4) + t155) * t148 + (t122 * t52 - t205 * t68 + (-t122 * t92 + t220) * qJD(4)) * t151 + t156 * g(3) + (t66 * mrSges(4,1) + t68 * mrSges(4,2) + (mrSges(3,1) * t149 + mrSges(3,2) * t152) * t200) * t143 + (-t75 - t76) * t66 + (-t142 * t213 - t26) * mrSges(4,2) + t157 * g(2) + t123 * t31 + t100 * t30 + t83 * t53 + t84 * t51 + t50 * t91 + t49 * t93 + t154 + t142 * mrSges(4,1) * t212 + (-t167 * g(2) - t178 * g(3) + t100 * t8 + t13 * t50 + t16 * t49 - t160 * t68 + t2 * t83 + t3 * t84 + (t179 - t66) * t27) * m(6) + ((t145 * t26 + t146 * t25) * pkin(2) - t125 * g(3) + t44 * t66 - t45 * t68 - t126 * g(2)) * m(4) + (t123 * t17 + (-t6 * t148 + t209 + (-t148 * t29 - t151 * t28) * qJD(4)) * t122 + (-t113 - t125) * g(3) - t177 * g(2) - t159 * t68 - t32 * t66) * m(5); m(4) * qJDD(3) + (t53 + t54) * t151 + (t51 + t52) * t148 + (-t206 * t148 + t205 * t151) * qJD(4) + m(5) * (t159 * qJD(4) + t148 * t5 + t151 * t6) + m(6) * (t160 * qJD(4) + t148 * t3 + t151 * t2) + (-m(4) - m(5) - m(6)) * g(1); t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t15 * t93 - t28 * t94 + t29 * t92 + t227 * t78 + t225 * t77 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t102 + t163) * g(1) + (t53 + (-g(1) * t151 + t2) * m(6)) * pkin(4) + (-m(6) * (-t13 + t15) + t91) * t16 + ((-t105 / 0.2e1 - t106 / 0.2e1 - t61 / 0.2e1 - t62 / 0.2e1 - t27 * mrSges(6,2) - t32 * mrSges(5,2) + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4) - t220) * t151 + (t59 / 0.2e1 + t60 / 0.2e1 - t27 * mrSges(6,1) - t32 * mrSges(5,1) + t232 * t193 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t27 - t75) * pkin(4) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t192 + t165) * t148) * t143 + (g(2) * t120 - g(3) * t121) * ((m(6) * pkin(4) - t230) * t148 + t207 * t151); (t148 * t91 - t151 * t93) * t143 + (-t160 * t143 + t164 + t8) * m(6) + t30;];
tau = t1;
