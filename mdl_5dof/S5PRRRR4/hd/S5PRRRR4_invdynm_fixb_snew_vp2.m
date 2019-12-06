% Calculate vector of cutting torques with Newton-Euler for
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:47
% EndTime: 2019-12-05 17:07:51
% DurationCPUTime: 3.78s
% Computational Cost: add. (65287->214), mult. (87051->277), div. (0->0), fcn. (55165->10), ass. (0->96)
t195 = sin(pkin(9));
t196 = cos(pkin(9));
t180 = t195 * g(1) - t196 * g(2);
t181 = -t196 * g(1) - t195 * g(2);
t200 = sin(qJ(2));
t204 = cos(qJ(2));
t161 = t204 * t180 - t200 * t181;
t156 = qJDD(2) * pkin(2) + t161;
t162 = t200 * t180 + t204 * t181;
t205 = qJD(2) ^ 2;
t157 = -t205 * pkin(2) + t162;
t199 = sin(qJ(3));
t203 = cos(qJ(3));
t142 = t199 * t156 + t203 * t157;
t191 = qJD(2) + qJD(3);
t187 = t191 ^ 2;
t189 = qJDD(2) + qJDD(3);
t138 = -t187 * pkin(3) + t189 * pkin(7) + t142;
t194 = -g(3) + qJDD(1);
t198 = sin(qJ(4));
t202 = cos(qJ(4));
t134 = -t198 * t138 + t202 * t194;
t221 = qJD(4) * t191;
t219 = t202 * t221;
t171 = t198 * t189 + t219;
t131 = (-t171 + t219) * pkin(8) + (t187 * t198 * t202 + qJDD(4)) * pkin(4) + t134;
t135 = t202 * t138 + t198 * t194;
t172 = t202 * t189 - t198 * t221;
t223 = t191 * t198;
t179 = qJD(4) * pkin(4) - pkin(8) * t223;
t193 = t202 ^ 2;
t132 = -t193 * t187 * pkin(4) + t172 * pkin(8) - qJD(4) * t179 + t135;
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t129 = t201 * t131 - t197 * t132;
t166 = (-t197 * t198 + t201 * t202) * t191;
t147 = t166 * qJD(5) + t201 * t171 + t197 * t172;
t167 = (t197 * t202 + t198 * t201) * t191;
t152 = -t166 * mrSges(6,1) + t167 * mrSges(6,2);
t190 = qJD(4) + qJD(5);
t159 = -t190 * mrSges(6,2) + t166 * mrSges(6,3);
t188 = qJDD(4) + qJDD(5);
t126 = m(6) * t129 + t188 * mrSges(6,1) - t147 * mrSges(6,3) - t167 * t152 + t190 * t159;
t130 = t197 * t131 + t201 * t132;
t146 = -t167 * qJD(5) - t197 * t171 + t201 * t172;
t160 = t190 * mrSges(6,1) - t167 * mrSges(6,3);
t127 = m(6) * t130 - t188 * mrSges(6,2) + t146 * mrSges(6,3) + t166 * t152 - t190 * t160;
t117 = t201 * t126 + t197 * t127;
t164 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t198 + Ifges(5,2) * t202) * t191;
t165 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t198 + Ifges(5,4) * t202) * t191;
t149 = Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * t190;
t150 = Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * t190;
t210 = -mrSges(6,1) * t129 + mrSges(6,2) * t130 - Ifges(6,5) * t147 - Ifges(6,6) * t146 - Ifges(6,3) * t188 - t167 * t149 + t166 * t150;
t224 = mrSges(5,1) * t134 - mrSges(5,2) * t135 + Ifges(5,5) * t171 + Ifges(5,6) * t172 + Ifges(5,3) * qJDD(4) + pkin(4) * t117 + (t198 * t164 - t202 * t165) * t191 - t210;
t170 = (-mrSges(5,1) * t202 + mrSges(5,2) * t198) * t191;
t222 = t191 * t202;
t177 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t222;
t115 = m(5) * t134 + qJDD(4) * mrSges(5,1) - t171 * mrSges(5,3) + qJD(4) * t177 - t170 * t223 + t117;
t176 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t223;
t215 = -t197 * t126 + t201 * t127;
t116 = m(5) * t135 - qJDD(4) * mrSges(5,2) + t172 * mrSges(5,3) - qJD(4) * t176 + t170 * t222 + t215;
t216 = -t198 * t115 + t202 * t116;
t109 = m(4) * t142 - t187 * mrSges(4,1) - t189 * mrSges(4,2) + t216;
t141 = t203 * t156 - t199 * t157;
t213 = -t189 * pkin(3) - t141;
t137 = -t187 * pkin(7) + t213;
t133 = t179 * t223 - t172 * pkin(4) + (-pkin(8) * t193 - pkin(7)) * t187 + t213;
t211 = m(6) * t133 - t146 * mrSges(6,1) + t147 * mrSges(6,2) - t166 * t159 + t167 * t160;
t207 = -m(5) * t137 + t172 * mrSges(5,1) - t171 * mrSges(5,2) - t176 * t223 + t177 * t222 - t211;
t121 = m(4) * t141 + t189 * mrSges(4,1) - t187 * mrSges(4,2) + t207;
t104 = t199 * t109 + t203 * t121;
t101 = m(3) * t161 + qJDD(2) * mrSges(3,1) - t205 * mrSges(3,2) + t104;
t217 = t203 * t109 - t199 * t121;
t102 = m(3) * t162 - t205 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t217;
t94 = t204 * t101 + t200 * t102;
t111 = t202 * t115 + t198 * t116;
t220 = m(4) * t194 + t111;
t218 = -t200 * t101 + t204 * t102;
t148 = Ifges(6,5) * t167 + Ifges(6,6) * t166 + Ifges(6,3) * t190;
t118 = -mrSges(6,1) * t133 + mrSges(6,3) * t130 + Ifges(6,4) * t147 + Ifges(6,2) * t146 + Ifges(6,6) * t188 - t167 * t148 + t190 * t150;
t119 = mrSges(6,2) * t133 - mrSges(6,3) * t129 + Ifges(6,1) * t147 + Ifges(6,4) * t146 + Ifges(6,5) * t188 + t166 * t148 - t190 * t149;
t163 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t198 + Ifges(5,6) * t202) * t191;
t106 = mrSges(5,2) * t137 - mrSges(5,3) * t134 + Ifges(5,1) * t171 + Ifges(5,4) * t172 + Ifges(5,5) * qJDD(4) - pkin(8) * t117 - qJD(4) * t164 - t197 * t118 + t201 * t119 + t163 * t222;
t97 = -mrSges(5,1) * t137 + mrSges(5,3) * t135 + Ifges(5,4) * t171 + Ifges(5,2) * t172 + Ifges(5,6) * qJDD(4) - pkin(4) * t211 + pkin(8) * t215 + qJD(4) * t165 + t201 * t118 + t197 * t119 - t163 * t223;
t212 = mrSges(4,1) * t141 - mrSges(4,2) * t142 + Ifges(4,3) * t189 + pkin(3) * t207 + pkin(7) * t216 + t198 * t106 + t202 * t97;
t209 = mrSges(3,1) * t161 - mrSges(3,2) * t162 + Ifges(3,3) * qJDD(2) + pkin(2) * t104 + t212;
t208 = mrSges(2,1) * t180 - mrSges(2,2) * t181 + pkin(1) * t94 + t209;
t95 = -mrSges(4,1) * t194 + mrSges(4,3) * t142 + t187 * Ifges(4,5) + Ifges(4,6) * t189 - pkin(3) * t111 - t224;
t92 = m(2) * t181 + t218;
t91 = m(2) * t180 + t94;
t90 = mrSges(4,2) * t194 - mrSges(4,3) * t141 + Ifges(4,5) * t189 - t187 * Ifges(4,6) - pkin(7) * t111 + t202 * t106 - t198 * t97;
t89 = mrSges(3,2) * t194 - mrSges(3,3) * t161 + Ifges(3,5) * qJDD(2) - t205 * Ifges(3,6) - pkin(6) * t104 - t199 * t95 + t203 * t90;
t88 = -mrSges(3,1) * t194 + mrSges(3,3) * t162 + t205 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t220 + pkin(6) * t217 + t199 * t90 + t203 * t95;
t87 = mrSges(2,2) * t194 - mrSges(2,3) * t180 - pkin(5) * t94 - t200 * t88 + t204 * t89;
t86 = -mrSges(2,1) * t194 + mrSges(2,3) * t181 + t200 * t89 + t204 * t88 - pkin(1) * (m(3) * t194 + t220) + pkin(5) * t218;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t196 * t87 - t195 * t86 - qJ(1) * (t195 * t92 + t196 * t91), t87, t89, t90, t106, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t195 * t87 + t196 * t86 + qJ(1) * (-t195 * t91 + t196 * t92), t86, t88, t95, t97, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t208, t208, t209, t212, t224, -t210;];
m_new = t1;
