% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:54
% EndTime: 2020-01-03 11:53:59
% DurationCPUTime: 4.17s
% Computational Cost: add. (72311->225), mult. (97595->288), div. (0->0), fcn. (55165->10), ass. (0->98)
t205 = sin(qJ(1));
t209 = cos(qJ(1));
t184 = -t209 * g(2) - t205 * g(3);
t177 = qJDD(1) * pkin(1) + t184;
t183 = -t205 * g(2) + t209 * g(3);
t210 = qJD(1) ^ 2;
t178 = -t210 * pkin(1) + t183;
t200 = sin(pkin(9));
t201 = cos(pkin(9));
t160 = t201 * t177 - t200 * t178;
t157 = qJDD(1) * pkin(2) + t160;
t161 = t200 * t177 + t201 * t178;
t158 = -t210 * pkin(2) + t161;
t204 = sin(qJ(3));
t208 = cos(qJ(3));
t142 = t204 * t157 + t208 * t158;
t195 = qJD(1) + qJD(3);
t191 = t195 ^ 2;
t193 = qJDD(1) + qJDD(3);
t139 = -t191 * pkin(3) + t193 * pkin(7) + t142;
t199 = -g(1) + qJDD(2);
t203 = sin(qJ(4));
t207 = cos(qJ(4));
t135 = -t203 * t139 + t207 * t199;
t226 = qJD(4) * t195;
t224 = t207 * t226;
t172 = t203 * t193 + t224;
t132 = (-t172 + t224) * pkin(8) + (t191 * t203 * t207 + qJDD(4)) * pkin(4) + t135;
t136 = t207 * t139 + t203 * t199;
t173 = t207 * t193 - t203 * t226;
t228 = t195 * t203;
t181 = qJD(4) * pkin(4) - pkin(8) * t228;
t198 = t207 ^ 2;
t133 = -t198 * t191 * pkin(4) + t173 * pkin(8) - qJD(4) * t181 + t136;
t202 = sin(qJ(5));
t206 = cos(qJ(5));
t130 = t206 * t132 - t202 * t133;
t167 = (-t202 * t203 + t206 * t207) * t195;
t148 = t167 * qJD(5) + t206 * t172 + t202 * t173;
t168 = (t202 * t207 + t203 * t206) * t195;
t153 = -t167 * mrSges(6,1) + t168 * mrSges(6,2);
t194 = qJD(4) + qJD(5);
t162 = -t194 * mrSges(6,2) + t167 * mrSges(6,3);
t192 = qJDD(4) + qJDD(5);
t127 = m(6) * t130 + t192 * mrSges(6,1) - t148 * mrSges(6,3) - t168 * t153 + t194 * t162;
t131 = t202 * t132 + t206 * t133;
t147 = -t168 * qJD(5) - t202 * t172 + t206 * t173;
t163 = t194 * mrSges(6,1) - t168 * mrSges(6,3);
t128 = m(6) * t131 - t192 * mrSges(6,2) + t147 * mrSges(6,3) + t167 * t153 - t194 * t163;
t118 = t206 * t127 + t202 * t128;
t165 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t203 + Ifges(5,2) * t207) * t195;
t166 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t203 + Ifges(5,4) * t207) * t195;
t150 = Ifges(6,4) * t168 + Ifges(6,2) * t167 + Ifges(6,6) * t194;
t151 = Ifges(6,1) * t168 + Ifges(6,4) * t167 + Ifges(6,5) * t194;
t215 = -mrSges(6,1) * t130 + mrSges(6,2) * t131 - Ifges(6,5) * t148 - Ifges(6,6) * t147 - Ifges(6,3) * t192 - t168 * t150 + t167 * t151;
t229 = mrSges(5,1) * t135 - mrSges(5,2) * t136 + Ifges(5,5) * t172 + Ifges(5,6) * t173 + Ifges(5,3) * qJDD(4) + pkin(4) * t118 + (t203 * t165 - t207 * t166) * t195 - t215;
t171 = (-mrSges(5,1) * t207 + mrSges(5,2) * t203) * t195;
t227 = t195 * t207;
t180 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t227;
t116 = m(5) * t135 + qJDD(4) * mrSges(5,1) - t172 * mrSges(5,3) + qJD(4) * t180 - t171 * t228 + t118;
t179 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t228;
t220 = -t202 * t127 + t206 * t128;
t117 = m(5) * t136 - qJDD(4) * mrSges(5,2) + t173 * mrSges(5,3) - qJD(4) * t179 + t171 * t227 + t220;
t221 = -t203 * t116 + t207 * t117;
t110 = m(4) * t142 - t191 * mrSges(4,1) - t193 * mrSges(4,2) + t221;
t141 = t208 * t157 - t204 * t158;
t218 = -t193 * pkin(3) - t141;
t138 = -t191 * pkin(7) + t218;
t134 = t181 * t228 - t173 * pkin(4) + (-pkin(8) * t198 - pkin(7)) * t191 + t218;
t216 = m(6) * t134 - t147 * mrSges(6,1) + t148 * mrSges(6,2) - t167 * t162 + t168 * t163;
t212 = -m(5) * t138 + t173 * mrSges(5,1) - t172 * mrSges(5,2) - t179 * t228 + t180 * t227 - t216;
t122 = m(4) * t141 + t193 * mrSges(4,1) - t191 * mrSges(4,2) + t212;
t105 = t204 * t110 + t208 * t122;
t102 = m(3) * t160 + qJDD(1) * mrSges(3,1) - t210 * mrSges(3,2) + t105;
t222 = t208 * t110 - t204 * t122;
t103 = m(3) * t161 - t210 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t222;
t95 = t201 * t102 + t200 * t103;
t112 = t207 * t116 + t203 * t117;
t225 = m(4) * t199 + t112;
t223 = -t200 * t102 + t201 * t103;
t149 = Ifges(6,5) * t168 + Ifges(6,6) * t167 + Ifges(6,3) * t194;
t119 = -mrSges(6,1) * t134 + mrSges(6,3) * t131 + Ifges(6,4) * t148 + Ifges(6,2) * t147 + Ifges(6,6) * t192 - t168 * t149 + t194 * t151;
t120 = mrSges(6,2) * t134 - mrSges(6,3) * t130 + Ifges(6,1) * t148 + Ifges(6,4) * t147 + Ifges(6,5) * t192 + t167 * t149 - t194 * t150;
t164 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t203 + Ifges(5,6) * t207) * t195;
t107 = mrSges(5,2) * t138 - mrSges(5,3) * t135 + Ifges(5,1) * t172 + Ifges(5,4) * t173 + Ifges(5,5) * qJDD(4) - pkin(8) * t118 - qJD(4) * t165 - t202 * t119 + t206 * t120 + t164 * t227;
t98 = -mrSges(5,1) * t138 + mrSges(5,3) * t136 + Ifges(5,4) * t172 + Ifges(5,2) * t173 + Ifges(5,6) * qJDD(4) - pkin(4) * t216 + pkin(8) * t220 + qJD(4) * t166 + t206 * t119 + t202 * t120 - t164 * t228;
t217 = mrSges(4,1) * t141 - mrSges(4,2) * t142 + Ifges(4,3) * t193 + pkin(3) * t212 + pkin(7) * t221 + t203 * t107 + t207 * t98;
t214 = mrSges(3,1) * t160 - mrSges(3,2) * t161 + Ifges(3,3) * qJDD(1) + pkin(2) * t105 + t217;
t213 = mrSges(2,1) * t184 - mrSges(2,2) * t183 + Ifges(2,3) * qJDD(1) + pkin(1) * t95 + t214;
t96 = -mrSges(4,1) * t199 + mrSges(4,3) * t142 + t191 * Ifges(4,5) + Ifges(4,6) * t193 - pkin(3) * t112 - t229;
t93 = m(2) * t184 + qJDD(1) * mrSges(2,1) - t210 * mrSges(2,2) + t95;
t92 = m(2) * t183 - t210 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t223;
t91 = mrSges(4,2) * t199 - mrSges(4,3) * t141 + Ifges(4,5) * t193 - t191 * Ifges(4,6) - pkin(7) * t112 + t207 * t107 - t203 * t98;
t90 = mrSges(3,2) * t199 - mrSges(3,3) * t160 + Ifges(3,5) * qJDD(1) - t210 * Ifges(3,6) - pkin(6) * t105 - t204 * t96 + t208 * t91;
t89 = -mrSges(3,1) * t199 + mrSges(3,3) * t161 + t210 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t225 + pkin(6) * t222 + t204 * t91 + t208 * t96;
t88 = -mrSges(2,2) * g(1) - mrSges(2,3) * t184 + Ifges(2,5) * qJDD(1) - t210 * Ifges(2,6) - qJ(2) * t95 - t200 * t89 + t201 * t90;
t87 = Ifges(2,6) * qJDD(1) + t210 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t183 + t200 * t90 + t201 * t89 - pkin(1) * (m(3) * t199 + t225) + qJ(2) * t223;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t213, t88, t90, t91, t107, t120; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t205 * t88 + t209 * t87 - pkin(5) * (t205 * t93 - t209 * t92), t87, t89, t96, t98, t119; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t209 * t88 + t205 * t87 + pkin(5) * (t205 * t92 + t209 * t93), t213, t214, t217, t229, -t215;];
m_new = t1;
