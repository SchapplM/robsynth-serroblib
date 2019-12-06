% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR6
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:40
% EndTime: 2019-12-05 17:09:46
% DurationCPUTime: 3.56s
% Computational Cost: add. (59876->216), mult. (76028->277), div. (0->0), fcn. (46790->10), ass. (0->96)
t196 = sin(pkin(9));
t223 = cos(pkin(9));
t183 = -g(1) * t223 - g(2) * t196;
t195 = -g(3) + qJDD(1);
t200 = sin(qJ(2));
t204 = cos(qJ(2));
t163 = -t183 * t200 + t204 * t195;
t161 = qJDD(2) * pkin(2) + t163;
t164 = t204 * t183 + t200 * t195;
t205 = qJD(2) ^ 2;
t162 = -pkin(2) * t205 + t164;
t199 = sin(qJ(3));
t203 = cos(qJ(3));
t150 = t199 * t161 + t203 * t162;
t193 = qJD(2) + qJD(3);
t189 = t193 ^ 2;
t191 = qJDD(2) + qJDD(3);
t142 = -pkin(3) * t189 + pkin(7) * t191 + t150;
t182 = g(1) * t196 - g(2) * t223;
t198 = sin(qJ(4));
t202 = cos(qJ(4));
t137 = -t198 * t142 - t202 * t182;
t220 = qJD(4) * t193;
t219 = t202 * t220;
t173 = t191 * t198 + t219;
t134 = (-t173 + t219) * pkin(8) + (t189 * t198 * t202 + qJDD(4)) * pkin(4) + t137;
t138 = t202 * t142 - t198 * t182;
t174 = t191 * t202 - t198 * t220;
t222 = t193 * t198;
t181 = qJD(4) * pkin(4) - pkin(8) * t222;
t194 = t202 ^ 2;
t135 = -pkin(4) * t189 * t194 + pkin(8) * t174 - qJD(4) * t181 + t138;
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t132 = t134 * t201 - t135 * t197;
t168 = (-t197 * t198 + t201 * t202) * t193;
t147 = qJD(5) * t168 + t173 * t201 + t174 * t197;
t169 = (t197 * t202 + t198 * t201) * t193;
t155 = -mrSges(6,1) * t168 + mrSges(6,2) * t169;
t192 = qJD(4) + qJD(5);
t156 = -mrSges(6,2) * t192 + mrSges(6,3) * t168;
t190 = qJDD(4) + qJDD(5);
t129 = m(6) * t132 + mrSges(6,1) * t190 - t147 * mrSges(6,3) - t155 * t169 + t156 * t192;
t133 = t134 * t197 + t135 * t201;
t146 = -qJD(5) * t169 - t173 * t197 + t174 * t201;
t157 = mrSges(6,1) * t192 - mrSges(6,3) * t169;
t130 = m(6) * t133 - mrSges(6,2) * t190 + t146 * mrSges(6,3) + t155 * t168 - t157 * t192;
t119 = t201 * t129 + t197 * t130;
t166 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t198 + Ifges(5,2) * t202) * t193;
t167 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t198 + Ifges(5,4) * t202) * t193;
t152 = Ifges(6,4) * t169 + Ifges(6,2) * t168 + Ifges(6,6) * t192;
t153 = Ifges(6,1) * t169 + Ifges(6,4) * t168 + Ifges(6,5) * t192;
t209 = -mrSges(6,1) * t132 + mrSges(6,2) * t133 - Ifges(6,5) * t147 - Ifges(6,6) * t146 - Ifges(6,3) * t190 - t169 * t152 + t168 * t153;
t225 = mrSges(5,1) * t137 - mrSges(5,2) * t138 + Ifges(5,5) * t173 + Ifges(5,6) * t174 + Ifges(5,3) * qJDD(4) + pkin(4) * t119 + (t166 * t198 - t167 * t202) * t193 - t209;
t224 = m(3) + m(4);
t221 = t193 * t202;
t172 = (-mrSges(5,1) * t202 + mrSges(5,2) * t198) * t193;
t179 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t221;
t117 = m(5) * t137 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t173 + qJD(4) * t179 - t172 * t222 + t119;
t178 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t222;
t215 = -t129 * t197 + t201 * t130;
t118 = m(5) * t138 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t174 - qJD(4) * t178 + t172 * t221 + t215;
t216 = -t117 * t198 + t118 * t202;
t108 = m(4) * t150 - mrSges(4,1) * t189 - mrSges(4,2) * t191 + t216;
t149 = t203 * t161 - t199 * t162;
t213 = -t191 * pkin(3) - t149;
t141 = -pkin(7) * t189 + t213;
t136 = t181 * t222 - t174 * pkin(4) + (-pkin(8) * t194 - pkin(7)) * t189 + t213;
t210 = m(6) * t136 - t146 * mrSges(6,1) + t147 * mrSges(6,2) - t168 * t156 + t157 * t169;
t207 = -m(5) * t141 + t174 * mrSges(5,1) - mrSges(5,2) * t173 - t178 * t222 + t179 * t221 - t210;
t123 = m(4) * t149 + mrSges(4,1) * t191 - mrSges(4,2) * t189 + t207;
t103 = t108 * t199 + t123 * t203;
t112 = t117 * t202 + t118 * t198;
t101 = m(3) * t163 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t205 + t103;
t217 = t108 * t203 - t123 * t199;
t102 = m(3) * t164 - mrSges(3,1) * t205 - qJDD(2) * mrSges(3,2) + t217;
t218 = -t101 * t200 + t102 * t204;
t151 = Ifges(6,5) * t169 + Ifges(6,6) * t168 + Ifges(6,3) * t192;
t120 = -mrSges(6,1) * t136 + mrSges(6,3) * t133 + Ifges(6,4) * t147 + Ifges(6,2) * t146 + Ifges(6,6) * t190 - t151 * t169 + t153 * t192;
t121 = mrSges(6,2) * t136 - mrSges(6,3) * t132 + Ifges(6,1) * t147 + Ifges(6,4) * t146 + Ifges(6,5) * t190 + t151 * t168 - t152 * t192;
t165 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t198 + Ifges(5,6) * t202) * t193;
t105 = mrSges(5,2) * t141 - mrSges(5,3) * t137 + Ifges(5,1) * t173 + Ifges(5,4) * t174 + Ifges(5,5) * qJDD(4) - pkin(8) * t119 - qJD(4) * t166 - t120 * t197 + t121 * t201 + t165 * t221;
t99 = -mrSges(5,1) * t141 + mrSges(5,3) * t138 + Ifges(5,4) * t173 + Ifges(5,2) * t174 + Ifges(5,6) * qJDD(4) - pkin(4) * t210 + pkin(8) * t215 + qJD(4) * t167 + t201 * t120 + t197 * t121 - t165 * t222;
t93 = -mrSges(4,2) * t182 - mrSges(4,3) * t149 + Ifges(4,5) * t191 - Ifges(4,6) * t189 - pkin(7) * t112 + t105 * t202 - t198 * t99;
t97 = mrSges(4,1) * t182 + mrSges(4,3) * t150 + t189 * Ifges(4,5) + Ifges(4,6) * t191 - pkin(3) * t112 - t225;
t89 = Ifges(3,6) * qJDD(2) + t205 * Ifges(3,5) + mrSges(3,1) * t182 + mrSges(3,3) * t164 + t199 * t93 + t203 * t97 - pkin(2) * (-m(4) * t182 + t112) + pkin(6) * t217;
t92 = -mrSges(3,2) * t182 - mrSges(3,3) * t163 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t205 - pkin(6) * t103 - t199 * t97 + t203 * t93;
t212 = -mrSges(2,2) * t183 + pkin(5) * t218 + t200 * t92 + t204 * t89 + pkin(1) * (t182 * t224 - t112) + mrSges(2,1) * t182;
t211 = mrSges(4,1) * t149 - mrSges(4,2) * t150 + Ifges(4,3) * t191 + pkin(3) * t207 + pkin(7) * t216 + t105 * t198 + t202 * t99;
t208 = mrSges(3,1) * t163 - mrSges(3,2) * t164 + Ifges(3,3) * qJDD(2) + pkin(2) * t103 + t211;
t109 = (m(2) + t224) * t182 - t112;
t96 = t101 * t204 + t102 * t200;
t94 = m(2) * t183 + t218;
t90 = -mrSges(2,1) * t195 + mrSges(2,3) * t183 - pkin(1) * t96 - t208;
t87 = mrSges(2,2) * t195 - mrSges(2,3) * t182 - pkin(5) * t96 - t200 * t89 + t204 * t92;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t223 * t87 - t196 * t90 - qJ(1) * (t109 * t223 + t196 * t94), t87, t92, t93, t105, t121; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t196 * t87 + t223 * t90 + qJ(1) * (-t109 * t196 + t223 * t94), t90, t89, t97, t99, t120; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t212, t212, t208, t211, t225, -t209;];
m_new = t1;
