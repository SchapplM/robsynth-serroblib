% Calculate vector of inverse dynamics joint torques for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:46
% DurationCPUTime: 4.53s
% Computational Cost: add. (1701->312), mult. (3913->410), div. (0->0), fcn. (2365->8), ass. (0->144)
t229 = mrSges(4,1) + mrSges(5,1);
t221 = Ifges(4,1) + Ifges(5,1);
t219 = Ifges(5,4) + Ifges(4,5);
t116 = qJ(2) + qJ(3);
t112 = sin(t116);
t113 = cos(t116);
t212 = -t229 * t113 + (mrSges(4,2) - mrSges(5,3)) * t112;
t118 = sin(qJ(2));
t120 = cos(qJ(2));
t88 = -mrSges(3,1) * t120 + mrSges(3,2) * t118;
t228 = t212 + t88;
t220 = -Ifges(4,4) + Ifges(5,5);
t218 = -Ifges(4,6) + Ifges(5,6);
t115 = qJD(2) + qJD(3);
t117 = sin(qJ(3));
t192 = cos(qJ(3));
t149 = t192 * t120;
t164 = t118 * qJD(1);
t70 = -qJD(1) * t149 + t117 * t164;
t181 = t70 * Ifges(5,5);
t67 = Ifges(4,4) * t70;
t79 = t117 * t120 + t118 * t192;
t71 = t79 * qJD(1);
t227 = t219 * t115 + t221 * t71 + t181 - t67;
t165 = qJD(3) * t117;
t122 = -pkin(6) - pkin(5);
t91 = t122 * t120;
t82 = qJD(1) * t91;
t153 = t192 * t82;
t90 = t122 * t118;
t81 = qJD(1) * t90;
t51 = t117 * t81 - t153;
t226 = -pkin(2) * t165 + t51;
t162 = qJD(1) * qJD(2);
t84 = qJDD(1) * t120 - t118 * t162;
t225 = t113 * (-m(5) * qJ(4) - mrSges(5,3));
t197 = t79 / 0.2e1;
t224 = t84 / 0.2e1;
t223 = -m(5) - m(4);
t193 = t120 / 0.2e1;
t60 = -mrSges(4,2) * t115 - mrSges(4,3) * t70;
t182 = t70 * mrSges(5,2);
t63 = mrSges(5,3) * t115 - t182;
t217 = -t60 - t63;
t180 = t71 * mrSges(4,3);
t216 = -mrSges(5,2) * t71 + t229 * t115 - t180;
t215 = t120 * Ifges(3,2);
t214 = t113 * pkin(3) + t112 * qJ(4);
t163 = t120 * qJD(1);
t187 = pkin(5) * t120;
t188 = pkin(5) * t118;
t211 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t164) * t187 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t163) * t188;
t76 = t84 * pkin(5);
t85 = qJDD(1) * t118 + t120 * t162;
t77 = t85 * pkin(5);
t210 = t118 * t77 + t120 * t76;
t119 = sin(qJ(1));
t121 = cos(qJ(1));
t209 = g(1) * t121 + g(2) * t119;
t174 = t117 * t82;
t75 = qJD(2) * pkin(2) + t81;
t48 = t192 * t75 + t174;
t208 = -t48 + qJD(4);
t206 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t49 = t117 * t75 - t153;
t56 = qJDD(2) * pkin(2) - pkin(6) * t85 - t77;
t57 = pkin(6) * t84 + t76;
t7 = -qJD(3) * t49 - t117 * t57 + t192 * t56;
t205 = m(3) * pkin(1) + mrSges(2,1) - t228;
t201 = -t70 / 0.2e1;
t200 = t70 / 0.2e1;
t198 = t71 / 0.2e1;
t194 = t115 / 0.2e1;
t191 = pkin(2) * t117;
t190 = pkin(2) * t118;
t189 = pkin(2) * t120;
t41 = t115 * qJ(4) + t49;
t184 = t41 * mrSges(5,2);
t183 = t48 * mrSges(4,3);
t179 = t71 * Ifges(4,4);
t178 = mrSges(4,2) * t113;
t177 = Ifges(3,4) * t118;
t176 = Ifges(3,4) * t120;
t170 = qJDD(1) * pkin(1);
t166 = qJD(2) * t118;
t158 = t192 * pkin(2);
t157 = pkin(2) * t164;
t156 = pkin(2) * t166;
t152 = t119 * t225;
t151 = t121 * t225;
t107 = pkin(1) + t189;
t150 = qJD(2) * t122;
t148 = qJD(3) * t192;
t114 = qJDD(2) + qJDD(3);
t131 = -t117 * t118 + t149;
t126 = t131 * qJD(3);
t30 = qJD(1) * t126 + t117 * t84 + t192 * t85;
t17 = -t114 * mrSges(5,1) + t30 * mrSges(5,2);
t145 = t162 / 0.2e1;
t142 = pkin(2) * t148;
t141 = t120 * t150;
t42 = pkin(3) * t71 + qJ(4) * t70;
t89 = t107 * qJD(1);
t65 = -pkin(2) * t84 - t170;
t140 = mrSges(3,1) * t118 + mrSges(3,2) * t120;
t137 = t177 + t215;
t136 = Ifges(3,5) * t120 - Ifges(3,6) * t118;
t133 = t117 * t91 + t192 * t90;
t59 = t117 * t90 - t192 * t91;
t132 = pkin(1) * t140;
t6 = t117 * t56 + t75 * t148 + t165 * t82 + t192 * t57;
t130 = t118 * (Ifges(3,1) * t120 - t177);
t127 = t79 * qJD(3);
t125 = t178 + (m(5) * pkin(3) + t229) * t112;
t124 = m(5) * (-pkin(3) * t112 - t190) - t112 * mrSges(5,1);
t3 = qJ(4) * t114 + qJD(4) * t115 + t6;
t31 = qJD(1) * t127 + t117 * t85 - t192 * t84;
t32 = pkin(3) * t70 - qJ(4) * t71 - t89;
t33 = -t115 * pkin(3) + t208;
t66 = Ifges(5,5) * t71;
t37 = t115 * Ifges(5,6) + t70 * Ifges(5,3) + t66;
t38 = -t70 * Ifges(4,2) + t115 * Ifges(4,6) + t179;
t4 = -t114 * pkin(3) + qJDD(4) - t7;
t123 = t182 * t33 - t32 * (mrSges(5,1) * t71 + mrSges(5,3) * t70) + t89 * (mrSges(4,1) * t71 - mrSges(4,2) * t70) - t6 * mrSges(4,2) + t7 * mrSges(4,1) + t3 * mrSges(5,3) - t4 * mrSges(5,1) - t70 * t183 + t71 * t184 + t38 * t198 + (Ifges(5,3) * t71 - t181) * t201 + t218 * t31 + t219 * t30 - (t218 * t71 - t219 * t70) * t115 / 0.2e1 + (Ifges(5,2) + Ifges(4,3)) * t114 + (-Ifges(4,2) * t71 + t227 - t67) * t200 - (-t221 * t70 - t179 + t37 + t66) * t71 / 0.2e1;
t109 = Ifges(3,4) * t163;
t106 = -t158 - pkin(3);
t103 = qJ(4) + t191;
t97 = t142 + qJD(4);
t83 = t118 * t150;
t69 = Ifges(3,1) * t164 + Ifges(3,5) * qJD(2) + t109;
t68 = Ifges(3,6) * qJD(2) + qJD(1) * t137;
t54 = qJD(2) * t79 + t127;
t53 = qJD(2) * t131 + t126;
t52 = t192 * t81 + t174;
t47 = -pkin(3) * t131 - qJ(4) * t79 - t107;
t46 = mrSges(4,1) * t70 + mrSges(4,2) * t71;
t45 = mrSges(5,1) * t70 - mrSges(5,3) * t71;
t34 = t42 + t157;
t19 = -mrSges(5,2) * t31 + mrSges(5,3) * t114;
t18 = -mrSges(4,2) * t114 - mrSges(4,3) * t31;
t16 = mrSges(4,1) * t114 - mrSges(4,3) * t30;
t8 = pkin(3) * t54 - qJ(4) * t53 - qJD(4) * t79 + t156;
t1 = pkin(3) * t31 - qJ(4) * t30 - qJD(4) * t71 + t65;
t2 = [-(-m(4) * t7 + m(5) * t4 - t16 + t17) * t133 + m(4) * (-t107 * t65 - t156 * t89) + (m(4) * t6 + m(5) * t3 + t18 + t19) * t59 + (-t107 * mrSges(4,1) + t47 * mrSges(5,1) - (Ifges(4,2) + Ifges(5,3)) * t131 + 0.2e1 * t220 * t197) * t31 + (-m(4) * t48 + m(5) * t33 - t216) * (qJD(3) * t59 + t117 * t83 - t141 * t192) + (m(4) * t49 + m(5) * t41 - t217) * (qJD(3) * t133 + t117 * t141 + t192 * t83) + (t187 * t84 + t188 * t85 + t210) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t210) + (t69 * t193 + t136 * qJD(2) / 0.2e1 - t211) * qJD(2) - t88 * t170 - t132 * t162 - t68 * t166 / 0.2e1 + (-t184 - t89 * mrSges(4,1) - t38 / 0.2e1 + t32 * mrSges(5,1) + t37 / 0.2e1 - t49 * mrSges(4,3) + Ifges(5,3) * t200 - Ifges(4,2) * t201 + t220 * t198 + t218 * t194) * t54 + (t114 * t219 + t221 * t30) * t197 + (Ifges(3,4) * t85 + Ifges(3,2) * t84) * t193 + (Ifges(3,1) * t85 + Ifges(3,4) * t224 + Ifges(3,5) * qJDD(2) - t145 * t215) * t118 + (-mrSges(3,1) * t188 - mrSges(3,2) * t187 + 0.2e1 * Ifges(3,6) * t193) * qJDD(2) + (t120 * t176 + t130) * t145 + m(5) * (t1 * t47 + t32 * t8) - pkin(1) * (-mrSges(3,1) * t84 + mrSges(3,2) * t85) + t8 * t45 + t46 * t156 + (t131 * t3 + t4 * t79) * mrSges(5,2) + t1 * (-mrSges(5,1) * t131 - mrSges(5,3) * t79) + t65 * (-mrSges(4,1) * t131 + mrSges(4,2) * t79) + (t131 * t6 - t7 * t79) * mrSges(4,3) + (-t131 * t218 + t219 * t79) * t114 / 0.2e1 + (-t131 * t220 + t221 * t79) * t30 / 0.2e1 + t131 * (Ifges(4,4) * t30 + Ifges(4,6) * t114) / 0.2e1 - t131 * (Ifges(5,5) * t30 + Ifges(5,6) * t114) / 0.2e1 + (t223 * (t121 * t107 - t119 * t122) + (-m(5) * t214 - t205) * t121 + t206 * t119) * g(2) + ((-t122 * t223 + t206) * t121 + (-m(5) * (-t107 - t214) + m(4) * t107 + t205) * t119) * g(1) - t47 * mrSges(5,3) * t30 - t107 * mrSges(4,2) * t30 + t85 * t176 / 0.2e1 + (t227 / 0.2e1 - mrSges(4,2) * t89 + t33 * mrSges(5,2) - mrSges(5,3) * t32 + Ifges(4,4) * t201 + Ifges(5,5) * t200 + t194 * t219 + t198 * t221 - t183) * t53 + Ifges(2,3) * qJDD(1) + t137 * t224; -(-Ifges(3,2) * t164 + t109 + t69) * t163 / 0.2e1 + ((t192 * t7 + t117 * t6 + (-t117 * t48 + t192 * t49) * qJD(3)) * pkin(2) + t89 * t157 + t48 * t51 - t49 * t52) * m(4) + (m(4) * t190 + mrSges(4,1) * t112 + t140 + t178) * t209 + t217 * t52 + t226 * t216 + (t103 * t3 + t106 * t4 - t32 * t34 + (t97 - t52) * t41 - t226 * t33) * m(5) + t68 * t164 / 0.2e1 - t136 * t162 / 0.2e1 - t46 * t157 - g(1) * (t121 * t124 - t151) - g(2) * (t119 * t124 - t152) + t123 + (t211 + (-t130 / 0.2e1 + t132) * qJD(1)) * qJD(1) + t97 * t63 + t103 * t19 + t106 * t17 + Ifges(3,6) * t84 + Ifges(3,5) * t85 - t76 * mrSges(3,2) - t77 * mrSges(3,1) - t34 * t45 + t60 * t142 + (-m(4) * t189 - m(5) * (t214 + t189) + t228) * g(3) + t16 * t158 + t49 * t180 + t18 * t191 + Ifges(3,3) * qJDD(2); (t121 * t125 + t151) * g(1) + (t119 * t125 + t152) * g(2) + (t180 + t216) * t49 + t217 * t48 + t123 + t212 * g(3) + qJD(4) * t63 - t42 * t45 - pkin(3) * t17 + qJ(4) * t19 + (-pkin(3) * t4 - g(3) * t214 + qJ(4) * t3 + t208 * t41 - t32 * t42 - t33 * t49) * m(5); -t115 * t63 + t71 * t45 + (g(3) * t113 - t112 * t209 - t41 * t115 + t32 * t71 + t4) * m(5) + t17;];
tau = t2;
