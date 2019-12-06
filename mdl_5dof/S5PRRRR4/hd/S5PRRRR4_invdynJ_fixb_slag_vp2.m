% Calculate vector of inverse dynamics joint torques for
% S5PRRRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:46
% DurationCPUTime: 2.01s
% Computational Cost: add. (2764->302), mult. (4181->412), div. (0->0), fcn. (2432->12), ass. (0->152)
t154 = sin(qJ(4));
t206 = Ifges(5,4) * t154;
t230 = t206 / 0.2e1;
t159 = -pkin(8) - pkin(7);
t179 = qJD(4) * t159;
t103 = t154 * t179;
t157 = cos(qJ(4));
t104 = t157 * t179;
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t158 = cos(qJ(3));
t204 = pkin(2) * qJD(2);
t183 = t158 * t204;
t114 = t159 * t154;
t146 = t157 * pkin(8);
t115 = pkin(7) * t157 + t146;
t67 = t114 * t153 + t115 * t156;
t99 = t153 * t157 + t154 * t156;
t229 = -qJD(5) * t67 - t103 * t153 + t104 * t156 + t99 * t183;
t66 = t114 * t156 - t115 * t153;
t98 = -t153 * t154 + t156 * t157;
t228 = qJD(5) * t66 + t103 * t156 + t104 * t153 - t98 * t183;
t205 = Ifges(5,4) * t157;
t227 = t157 * Ifges(5,2);
t148 = qJDD(2) + qJDD(3);
t155 = sin(qJ(3));
t203 = pkin(2) * qJD(3);
t180 = qJD(2) * t203;
t193 = pkin(2) * qJDD(2);
t95 = -t155 * t180 + t158 * t193;
t86 = -pkin(3) * t148 - t95;
t151 = qJD(2) + qJD(3);
t186 = qJD(4) * t154;
t178 = t151 * t186;
t89 = t148 * t157 - t178;
t185 = qJD(4) * t157;
t90 = t148 * t154 + t151 * t185;
t226 = m(5) * t86 - mrSges(5,1) * t89 + mrSges(5,2) * t90;
t113 = -t157 * mrSges(5,1) + t154 * mrSges(5,2);
t152 = qJ(4) + qJ(5);
t144 = sin(t152);
t145 = cos(t152);
t174 = t145 * mrSges(6,1) - t144 * mrSges(6,2);
t225 = t113 - t174;
t108 = pkin(7) * t151 + t155 * t204;
t143 = t157 * qJD(1);
t96 = t155 * t193 + t158 * t180;
t87 = pkin(7) * t148 + t96;
t40 = qJD(4) * t143 + t154 * qJDD(1) - t108 * t186 + t157 * t87;
t187 = qJD(1) * t154;
t80 = t108 * t157 + t187;
t41 = -qJD(4) * t80 + t157 * qJDD(1) - t154 * t87;
t224 = -t41 * t154 + t157 * t40;
t149 = pkin(9) + qJ(2);
t142 = qJ(3) + t149;
t129 = sin(t142);
t130 = cos(t142);
t220 = g(1) * t130 + g(2) * t129;
t223 = mrSges(4,2) - mrSges(6,3) - mrSges(5,3);
t222 = -mrSges(4,1) + t225;
t83 = t98 * t151;
t84 = t99 * t151;
t48 = -mrSges(6,1) * t83 + mrSges(6,2) * t84;
t133 = pkin(4) * t157 + pkin(3);
t85 = -t133 * t151 - t183;
t221 = -m(6) * t85 - t48;
t218 = t84 / 0.2e1;
t217 = m(5) + m(6);
t214 = mrSges(6,3) * t83;
t213 = Ifges(6,4) * t84;
t212 = pkin(2) * t158;
t138 = sin(t149);
t210 = g(1) * t138;
t208 = t84 * mrSges(6,3);
t132 = pkin(2) * t155 + pkin(7);
t207 = -pkin(8) - t132;
t176 = pkin(8) * t151 + t108;
t59 = t157 * t176 + t187;
t201 = t153 * t59;
t199 = t156 * t59;
t195 = Ifges(5,5) * qJD(4);
t194 = Ifges(5,6) * qJD(4);
t192 = t151 * t154;
t191 = t151 * t157;
t190 = t154 * t158;
t189 = t157 * t158;
t188 = t130 * pkin(3) + t129 * pkin(7);
t184 = m(2) + m(3) + m(4);
t182 = t158 * t203;
t181 = pkin(4) * t186;
t175 = qJD(4) * t207;
t173 = -t129 * t159 + t130 * t133;
t172 = mrSges(5,1) * t154 + mrSges(5,2) * t157;
t171 = mrSges(6,1) * t144 + mrSges(6,2) * t145;
t170 = t206 + t227;
t58 = -t154 * t176 + t143;
t54 = qJD(4) * pkin(4) + t58;
t27 = t156 * t54 - t201;
t28 = t153 * t54 + t199;
t93 = t207 * t154;
t94 = t132 * t157 + t146;
t51 = -t153 * t94 + t156 * t93;
t52 = t153 * t93 + t156 * t94;
t106 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t192;
t107 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t191;
t169 = -t154 * t106 + t157 * t107;
t168 = t99 * qJD(5);
t167 = t98 * qJD(5);
t109 = -pkin(3) * t151 - t183;
t79 = -t108 * t154 + t143;
t165 = m(5) * (t109 * t155 + t189 * t80 - t190 * t79);
t147 = qJDD(4) + qJDD(5);
t150 = qJD(4) + qJD(5);
t19 = qJDD(4) * pkin(4) - pkin(8) * t90 + t41;
t25 = pkin(8) * t89 + t40;
t3 = qJD(5) * t27 + t153 * t19 + t156 * t25;
t35 = t151 * t167 + t153 * t89 + t156 * t90;
t36 = -t151 * t168 - t153 * t90 + t156 * t89;
t4 = -qJD(5) * t28 - t153 * t25 + t156 * t19;
t43 = Ifges(6,2) * t83 + Ifges(6,6) * t150 + t213;
t76 = Ifges(6,4) * t83;
t44 = Ifges(6,1) * t84 + Ifges(6,5) * t150 + t76;
t164 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t27 * t214 + t43 * t218 - t85 * (mrSges(6,1) * t84 + mrSges(6,2) * t83) + Ifges(6,3) * t147 - t84 * (Ifges(6,1) * t83 - t213) / 0.2e1 + Ifges(6,6) * t36 + Ifges(6,5) * t35 - t150 * (Ifges(6,5) * t83 - Ifges(6,6) * t84) / 0.2e1 - (-Ifges(6,2) * t84 + t44 + t76) * t83 / 0.2e1;
t163 = t223 * t129 + t222 * t130;
t162 = (m(5) * pkin(3) + m(6) * t133 - t222) * t129 + (-m(5) * pkin(7) + m(6) * t159 + t223) * t130;
t74 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t89;
t75 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t90;
t161 = m(5) * (-t185 * t79 - t186 * t80 + t224) + t157 * t74 - t154 * t75 - t106 * t185 - t107 * t186;
t49 = -pkin(4) * t89 + t86;
t55 = qJD(4) * t98 + t167;
t56 = -qJD(4) * t99 - t168;
t81 = t151 * t170 + t194;
t116 = Ifges(5,4) * t191;
t82 = Ifges(5,1) * t192 + t116 + t195;
t160 = qJD(4) ^ 2 * (Ifges(5,5) * t157 - Ifges(5,6) * t154) / 0.2e1 + t109 * t172 * qJD(4) + Ifges(4,3) * t148 + t150 * (Ifges(6,5) * t55 + Ifges(6,6) * t56) / 0.2e1 + t86 * t113 + t95 * mrSges(4,1) - t96 * mrSges(4,2) + t85 * (-mrSges(6,1) * t56 + mrSges(6,2) * t55) + t83 * (Ifges(6,4) * t55 + Ifges(6,2) * t56) / 0.2e1 + t55 * t44 / 0.2e1 + t56 * t43 / 0.2e1 + (Ifges(5,1) * t157 - t206) * t178 / 0.2e1 - t81 * t186 / 0.2e1 + (Ifges(6,1) * t55 + Ifges(6,4) * t56) * t218 + (t227 / 0.2e1 + t230 + t170 / 0.2e1) * t89 + t90 * (t154 * Ifges(5,1) + t205) + (t151 * (-Ifges(5,2) * t154 + t205) + t82) * t185 / 0.2e1 + (-t27 * t55 + t28 * t56) * mrSges(6,3) + ((-t154 * t80 - t157 * t79) * qJD(4) + t224) * mrSges(5,3) + qJDD(4) * (Ifges(5,5) * t154 + Ifges(5,6) * t157) + (t49 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t35 + Ifges(6,4) * t36 + Ifges(6,5) * t147) * t99 + (-t49 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t35 + Ifges(6,2) * t36 + Ifges(6,6) * t147) * t98;
t139 = cos(t149);
t127 = pkin(2) * t139;
t112 = -t133 - t212;
t105 = t155 * t203 + t181;
t88 = t113 * t151;
t69 = -t154 * t182 + t157 * t175;
t68 = t154 * t175 + t157 * t182;
t63 = mrSges(6,1) * t150 - t208;
t62 = -mrSges(6,2) * t150 + t214;
t32 = t156 * t58 - t201;
t31 = -t153 * t58 - t199;
t24 = -mrSges(6,2) * t147 + mrSges(6,3) * t36;
t23 = mrSges(6,1) * t147 - mrSges(6,3) * t35;
t10 = -qJD(5) * t52 - t153 * t68 + t156 * t69;
t9 = qJD(5) * t51 + t153 * t69 + t156 * t68;
t7 = -mrSges(6,1) * t36 + mrSges(6,2) * t35;
t1 = [t154 * t74 + t157 * t75 + t98 * t23 + t99 * t24 + t55 * t62 + t56 * t63 + t169 * qJD(4) + m(5) * (t154 * t40 + t157 * t41 + (-t154 * t79 + t157 * t80) * qJD(4)) + m(6) * (t27 * t56 + t28 * t55 + t3 * t99 + t4 * t98) + t184 * qJDD(1) + (-t184 - t217) * g(3); t161 * t132 + ((t158 * mrSges(4,1) - t155 * mrSges(4,2)) * t148 + t217 * t210 + (-g(2) * t139 + t155 * t96 + t158 * t95 + t210) * m(4) + (t165 + t155 * t88 - t106 * t190 + t107 * t189 + (-t155 * mrSges(4,1) - t158 * mrSges(4,2)) * t151) * qJD(3)) * pkin(2) + t105 * t48 + t112 * t7 + t9 * t62 + t10 * t63 + t51 * t23 + t52 * t24 + (mrSges(3,1) * t138 + mrSges(3,2) * t139 + t162) * g(1) + (-m(6) * (t127 + t173) - m(5) * (t127 + t188) - mrSges(3,1) * t139 + mrSges(3,2) * t138 + t163) * g(2) + t160 + m(6) * (t10 * t27 + t105 * t85 + t112 * t49 + t28 * t9 + t3 * t52 + t4 * t51) + Ifges(3,3) * qJDD(2) + t226 * (-pkin(3) - t212); t162 * g(1) + (-m(5) * t188 + t163) * g(2) + t229 * t63 + t228 * t62 + t161 * pkin(7) - t133 * t7 + t66 * t23 + t67 * t24 + ((t151 * mrSges(4,2) - t169) * t158 - t165 + (t151 * mrSges(4,1) + t221 - t88) * t155) * t204 + t48 * t181 + t160 - t226 * pkin(3) + (-t173 * g(2) - t133 * t49 + t85 * t181 + t228 * t28 + t229 * t27 + t3 * t67 + t4 * t66) * m(6); t225 * g(3) - m(6) * (t27 * t31 + t28 * t32) + ((-t116 / 0.2e1 - t82 / 0.2e1 - t109 * mrSges(5,2) - t195 / 0.2e1 + t79 * mrSges(5,3)) * t157 + (t81 / 0.2e1 - t109 * mrSges(5,1) + t194 / 0.2e1 + t80 * mrSges(5,3) + (t230 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t157) * t151 + t221 * pkin(4)) * t154) * t151 + (t153 * t24 + t156 * t23 + (-g(3) * t157 + t153 * t3 + t220 * t154 + t156 * t4) * m(6) + (-t153 * t63 + t156 * t62 + (-t153 * t27 + t156 * t28) * m(6)) * qJD(5)) * pkin(4) + t164 + t80 * t106 - t79 * t107 + Ifges(5,5) * t90 + Ifges(5,6) * t89 - t32 * t62 - t31 * t63 - t40 * mrSges(5,2) + t41 * mrSges(5,1) + t28 * t208 + Ifges(5,3) * qJDD(4) + t220 * (t171 + t172); (t63 + t208) * t28 + t164 - g(3) * t174 - t27 * t62 + t220 * t171;];
tau = t1;
