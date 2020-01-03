% Calculate vector of inverse dynamics joint torques for
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:07:07
% DurationCPUTime: 5.62s
% Computational Cost: add. (3058->381), mult. (6312->513), div. (0->0), fcn. (4128->10), ass. (0->172)
t115 = sin(pkin(8));
t116 = cos(pkin(8));
t174 = t115 ^ 2 + t116 ^ 2;
t151 = t174 * mrSges(4,3);
t118 = -pkin(1) - qJ(3);
t239 = -qJD(1) * qJD(3) + qJDD(1) * t118;
t119 = sin(qJ(5));
t122 = cos(qJ(5));
t120 = sin(qJ(4));
t210 = cos(qJ(4));
t159 = t210 * t116;
t128 = -t120 * t115 + t159;
t129 = -t115 * t210 - t120 * t116;
t71 = t129 * qJD(1);
t49 = qJD(4) * t71 + qJDD(1) * t128;
t173 = qJD(1) * t115;
t72 = qJD(1) * t159 - t120 * t173;
t56 = qJD(4) * t122 - t119 * t72;
t26 = qJD(5) * t56 + qJDD(4) * t119 + t122 * t49;
t220 = t26 / 0.2e1;
t57 = qJD(4) * t119 + t122 * t72;
t27 = -qJD(5) * t57 + qJDD(4) * t122 - t119 * t49;
t219 = t27 / 0.2e1;
t50 = -qJD(1) * qJD(4) * t128 + qJDD(1) * t129;
t46 = qJDD(5) - t50;
t218 = t46 / 0.2e1;
t121 = sin(qJ(1));
t123 = cos(qJ(1));
t230 = -g(1) * t121 + g(2) * t123;
t238 = Ifges(6,1) * t220 + Ifges(6,4) * t219 + Ifges(6,5) * t218;
t237 = -m(5) - m(6);
t7 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t236 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t49 + t7;
t85 = qJD(1) * t118 + qJD(2);
t152 = -pkin(6) * qJD(1) + t85;
t63 = t152 * t115;
t64 = t152 * t116;
t37 = t120 * t64 + t210 * t63;
t32 = qJD(4) * pkin(7) + t37;
t101 = qJD(1) * qJ(2) + qJD(3);
t83 = pkin(3) * t173 + t101;
t33 = -pkin(4) * t71 - pkin(7) * t72 + t83;
t12 = -t119 * t32 + t122 * t33;
t235 = t12 * mrSges(6,1);
t13 = t119 * t33 + t122 * t32;
t234 = t13 * mrSges(6,2);
t195 = t72 * mrSges(5,3);
t191 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t56 - mrSges(6,2) * t57 - t195;
t233 = qJD(5) - t71;
t113 = qJD(1) * qJD(2);
t86 = -qJDD(1) * qJ(2) - t113;
t8 = mrSges(6,1) * t46 - mrSges(6,3) * t26;
t9 = -mrSges(6,2) * t46 + mrSges(6,3) * t27;
t232 = -t119 * t8 + t122 * t9;
t155 = qJD(4) * t210;
t172 = qJD(4) * t120;
t80 = qJDD(2) + t239;
t148 = -pkin(6) * qJDD(1) + t80;
t59 = t148 * t115;
t60 = t148 * t116;
t14 = t120 * t60 + t64 * t155 - t172 * t63 + t210 * t59;
t10 = qJDD(4) * pkin(7) + t14;
t166 = qJDD(1) * t115;
t84 = qJDD(3) - t86;
t75 = pkin(3) * t166 + t84;
t18 = -pkin(4) * t50 - pkin(7) * t49 + t75;
t1 = qJD(5) * t12 + t10 * t122 + t119 * t18;
t2 = -qJD(5) * t13 - t10 * t119 + t122 * t18;
t231 = t1 * t122 - t119 * t2;
t227 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t15 = -qJD(4) * t37 - t120 * t59 + t210 * t60;
t36 = -t120 * t63 + t210 * t64;
t31 = -qJD(4) * pkin(4) - t36;
t226 = -m(6) * t31 + t191;
t73 = -t115 * t155 - t116 * t172;
t74 = -t115 * t172 + t116 * t155;
t225 = -t128 * t15 + t129 * t14 - t36 * t73 - t37 * t74;
t224 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t111 = pkin(8) + qJ(4);
t103 = sin(t111);
t104 = cos(t111);
t145 = mrSges(5,1) * t103 + mrSges(5,2) * t104;
t147 = mrSges(4,1) * t115 + mrSges(4,2) * t116;
t223 = -t147 - mrSges(3,3) - t145 - m(6) * (pkin(4) * t103 - pkin(7) * t104) + t104 * mrSges(6,3) + mrSges(2,2);
t221 = m(4) * t101 + m(5) * t83 - mrSges(5,1) * t71 + mrSges(5,2) * t72 + t147 * qJD(1);
t217 = -t56 / 0.2e1;
t216 = -t57 / 0.2e1;
t215 = t57 / 0.2e1;
t214 = -t233 / 0.2e1;
t211 = t72 / 0.2e1;
t209 = mrSges(5,3) * t71;
t208 = Ifges(5,4) * t72;
t207 = Ifges(6,4) * t57;
t11 = -qJDD(4) * pkin(4) - t15;
t203 = t11 * t128;
t105 = t115 * pkin(3);
t192 = -pkin(6) + t118;
t165 = qJDD(1) * t116;
t190 = mrSges(4,1) * t166 + mrSges(4,2) * t165;
t189 = Ifges(6,4) * t119;
t188 = Ifges(6,4) * t122;
t187 = t119 * t71;
t186 = t119 * t128;
t185 = t122 * t71;
t184 = t122 * t73;
t183 = t122 * t128;
t182 = qJDD(1) * pkin(1);
t181 = t119 * t121;
t180 = t119 * t123;
t178 = t121 * t122;
t177 = t122 * t123;
t94 = qJ(2) + t105;
t175 = t123 * pkin(1) + t121 * qJ(2);
t171 = qJD(5) * t119;
t170 = qJD(5) * t122;
t169 = -m(4) + t237;
t164 = Ifges(6,5) * t26 + Ifges(6,6) * t27 + Ifges(6,3) * t46;
t162 = m(6) * pkin(7) + mrSges(6,3);
t22 = Ifges(6,2) * t56 + Ifges(6,6) * t233 + t207;
t161 = -t119 * t22 / 0.2e1;
t156 = -t50 * mrSges(5,1) + t49 * mrSges(5,2);
t154 = t170 / 0.2e1;
t107 = t123 * qJ(2);
t153 = -pkin(1) * t121 + t107;
t150 = t174 * t85;
t149 = t174 * t80;
t144 = -mrSges(6,1) * t122 + mrSges(6,2) * t119;
t143 = mrSges(6,1) * t119 + mrSges(6,2) * t122;
t142 = Ifges(6,1) * t122 - t189;
t141 = -Ifges(6,2) * t119 + t188;
t140 = Ifges(6,5) * t122 - Ifges(6,6) * t119;
t139 = t13 * t119 + t12 * t122;
t138 = -t12 * t119 + t13 * t122;
t34 = -mrSges(6,2) * t233 + mrSges(6,3) * t56;
t35 = mrSges(6,1) * t233 - mrSges(6,3) * t57;
t137 = -t119 * t34 - t122 * t35;
t47 = -pkin(4) * t129 - pkin(7) * t128 + t94;
t81 = t192 * t115;
t82 = t192 * t116;
t52 = t120 * t82 + t210 * t81;
t19 = -t119 * t52 + t122 * t47;
t20 = t119 * t47 + t122 * t52;
t133 = -t120 * t81 + t210 * t82;
t132 = -t119 * t73 - t128 * t170;
t131 = -t128 * t171 + t184;
t130 = t31 * t143;
t127 = m(6) * pkin(4) - t144;
t125 = -qJD(5) * t139 + t231;
t124 = qJD(1) ^ 2;
t117 = -pkin(6) - qJ(3);
t102 = qJDD(2) - t182;
t70 = t103 * t177 - t181;
t69 = t103 * t180 + t178;
t68 = t103 * t178 + t180;
t67 = -t103 * t181 + t177;
t65 = Ifges(5,4) * t71;
t61 = -qJD(4) * mrSges(5,2) + t209;
t53 = Ifges(6,4) * t56;
t48 = pkin(4) * t72 - pkin(7) * t71;
t43 = pkin(4) * t74 - pkin(7) * t73 + qJD(2);
t41 = t72 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t65;
t40 = t71 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t208;
t39 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t50;
t29 = qJD(3) * t129 + qJD(4) * t133;
t23 = Ifges(6,1) * t57 + Ifges(6,5) * t233 + t53;
t21 = t57 * Ifges(6,5) + t56 * Ifges(6,6) + Ifges(6,3) * t233;
t17 = t119 * t48 + t122 * t36;
t16 = -t119 * t36 + t122 * t48;
t6 = -qJD(5) * t20 - t119 * t29 + t122 * t43;
t5 = qJD(5) * t19 + t119 * t43 + t122 * t29;
t3 = t26 * Ifges(6,4) + t27 * Ifges(6,2) + t46 * Ifges(6,6);
t4 = [t225 * mrSges(5,3) + t221 * qJD(2) + t183 * t238 + t74 * t235 + t233 * (Ifges(6,5) * t131 + Ifges(6,6) * t132 + Ifges(6,3) * t74) / 0.2e1 - (-t75 * mrSges(5,2) - Ifges(5,1) * t49 - Ifges(5,4) * t50 - Ifges(5,5) * qJDD(4) - t140 * t218 - t141 * t219 - t142 * t220 + t154 * t22) * t128 - (Ifges(6,3) * t218 + Ifges(6,6) * t219 + Ifges(6,5) * t220 + t164 / 0.2e1 - Ifges(5,6) * qJDD(4) - Ifges(5,4) * t49 - Ifges(5,2) * t50 + t75 * mrSges(5,1) + t224) * t129 + m(5) * (t14 * t52 + t29 * t37 + t75 * t94) + m(6) * (t1 * t20 + t12 * t6 + t13 * t5 + t19 * t2) + t143 * t203 + qJ(2) * t190 + t23 * t184 / 0.2e1 - (-m(5) * t15 + m(6) * t11 + t236) * t133 + (-m(3) * t153 - m(4) * t107 - t70 * mrSges(6,1) + t69 * mrSges(6,2) + t237 * (t123 * t105 + t121 * t117 + t153) + (-m(4) * t118 - t227) * t121 + t223 * t123) * g(1) + (-t68 * mrSges(6,1) - t67 * mrSges(6,2) + (-m(4) - m(3)) * t175 + t237 * (t121 * t105 - t117 * t123 + t175) + (-m(4) * qJ(3) + t227) * t123 + t223 * t121) * g(2) - t74 * t234 + (Ifges(5,1) * t73 - Ifges(5,4) * t74) * t211 + (Ifges(6,1) * t131 + Ifges(6,4) * t132 + Ifges(6,5) * t74) * t215 + (Ifges(4,1) * t116 - Ifges(4,4) * t115) * t165 + (-m(5) * t36 - t226) * (qJD(3) * t128 + qJD(4) * t52) + m(3) * (-pkin(1) * t102 + (-t86 + t113) * qJ(2)) - (Ifges(4,4) * t116 - Ifges(4,2) * t115) * t166 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t84 * t147 + t56 * (Ifges(6,4) * t131 + Ifges(6,2) * t132 + Ifges(6,6) * t74) / 0.2e1 + t31 * (-mrSges(6,1) * t132 + mrSges(6,2) * t131) + (-t1 * t186 - t12 * t131 + t13 * t132 - t183 * t2) * mrSges(6,3) - (qJD(5) * t23 + t3) * t186 / 0.2e1 + m(4) * (qJ(2) * t84 - qJD(3) * t150 + t118 * t149) + (-t182 + t102) * mrSges(3,2) + t94 * t156 + t83 * (mrSges(5,1) * t74 + mrSges(5,2) * t73) - 0.2e1 * t86 * mrSges(3,3) + t73 * t41 / 0.2e1 + t74 * t21 / 0.2e1 + t71 * (Ifges(5,4) * t73 - Ifges(5,2) * t74) / 0.2e1 + qJD(4) * (Ifges(5,5) * t73 - Ifges(5,6) * t74) / 0.2e1 - t74 * t40 / 0.2e1 + t29 * t61 + t52 * t39 + t5 * t34 + t6 * t35 + t19 * t8 + t20 * t9 + (-t80 - t239) * t151 + t73 * t161; -t236 * t128 + t191 * t73 + (-m(3) * qJ(2) - mrSges(3,3)) * t124 + (-t119 * t35 + t122 * t34 + t61) * t74 + (mrSges(3,2) - t151) * qJDD(1) - (qJD(5) * t137 + t232 + t39) * t129 + m(6) * (-t125 * t129 + t138 * t74 - t31 * t73 - t203) + m(3) * t102 - m(5) * t225 + m(4) * t149 + (-m(6) * t139 + t137 - t221) * qJD(1) + t230 * (m(3) - t169); -t71 * t61 + t191 * t72 - t124 * t151 + (t233 * t34 + t8) * t122 + (-t233 * t35 + t9) * t119 + t156 + t190 + (g(1) * t123 + g(2) * t121) * t169 + (t1 * t119 + t122 * t2 + t138 * t233 - t31 * t72) * m(6) + (t36 * t72 - t37 * t71 + t75) * m(5) + (qJD(1) * t150 + t84) * m(4); (t103 * t127 - t104 * t162 + t145) * g(3) + ((-t171 + t187) * t13 + (-t170 + t185) * t12 + t231) * mrSges(6,3) + (m(6) * t125 - t170 * t35 - t171 * t34 + t232) * pkin(7) + (t195 + t226) * t37 + t119 * t238 + t72 * t234 + (t140 * t233 + t141 * t56 + t142 * t57) * qJD(5) / 0.2e1 + t22 * t187 / 0.2e1 + ((mrSges(5,1) + t127) * t104 + (-mrSges(5,2) + t162) * t103) * t230 - t72 * t235 + (Ifges(6,6) * t72 + t141 * t71) * t217 + (Ifges(6,5) * t119 + Ifges(6,6) * t122) * t218 + (Ifges(6,2) * t122 + t189) * t219 + (Ifges(6,1) * t119 + t188) * t220 + t40 * t211 + (Ifges(6,3) * t72 + t140 * t71) * t214 + (Ifges(6,5) * t72 + t142 * t71) * t216 + (t209 - t61) * t36 + (-pkin(4) * t11 - t12 * t16 - t13 * t17) * m(6) + (-t185 / 0.2e1 + t154) * t23 + t11 * t144 - t71 * t130 - (Ifges(5,1) * t71 - t208 + t21) * t72 / 0.2e1 - (-Ifges(5,2) * t72 + t41 + t65) * t71 / 0.2e1 + t122 * t3 / 0.2e1 - t83 * (mrSges(5,1) * t72 + mrSges(5,2) * t71) - qJD(4) * (Ifges(5,5) * t71 - Ifges(5,6) * t72) / 0.2e1 + Ifges(5,6) * t50 + Ifges(5,5) * t49 - t17 * t34 - t16 * t35 - t14 * mrSges(5,2) + t15 * mrSges(5,1) - pkin(4) * t7 + (t130 + t161) * qJD(5) + Ifges(5,3) * qJDD(4); -t31 * (mrSges(6,1) * t57 + mrSges(6,2) * t56) + (Ifges(6,1) * t56 - t207) * t216 + t22 * t215 + (Ifges(6,5) * t56 - Ifges(6,6) * t57) * t214 - t12 * t34 + t13 * t35 - g(1) * (mrSges(6,1) * t67 - mrSges(6,2) * t68) - g(2) * (mrSges(6,1) * t69 + mrSges(6,2) * t70) + g(3) * t143 * t104 + (t12 * t56 + t13 * t57) * mrSges(6,3) + t164 + (-Ifges(6,2) * t57 + t23 + t53) * t217 + t224;];
tau = t4;
