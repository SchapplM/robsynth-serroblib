% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:56
% EndTime: 2019-12-05 15:50:07
% DurationCPUTime: 5.43s
% Computational Cost: add. (72112->223), mult. (130217->293), div. (0->0), fcn. (88944->12), ass. (0->105)
t201 = sin(pkin(9));
t204 = cos(pkin(9));
t190 = t201 * g(1) - t204 * g(2);
t191 = -t204 * g(1) - t201 * g(2);
t199 = -g(3) + qJDD(1);
t208 = sin(qJ(2));
t205 = cos(pkin(5));
t211 = cos(qJ(2));
t229 = t205 * t211;
t202 = sin(pkin(5));
t231 = t202 * t211;
t158 = t190 * t229 - t208 * t191 + t199 * t231;
t156 = qJDD(2) * pkin(2) + t158;
t230 = t205 * t208;
t232 = t202 * t208;
t159 = t190 * t230 + t211 * t191 + t199 * t232;
t213 = qJD(2) ^ 2;
t157 = -t213 * pkin(2) + t159;
t200 = sin(pkin(10));
t203 = cos(pkin(10));
t152 = t200 * t156 + t203 * t157;
t149 = -t213 * pkin(3) + qJDD(2) * pkin(7) + t152;
t173 = -t202 * t190 + t205 * t199;
t172 = qJDD(3) + t173;
t207 = sin(qJ(4));
t210 = cos(qJ(4));
t146 = t210 * t149 + t207 * t172;
t186 = (-pkin(4) * t210 - pkin(8) * t207) * qJD(2);
t212 = qJD(4) ^ 2;
t226 = t210 * qJD(2);
t143 = -t212 * pkin(4) + qJDD(4) * pkin(8) + t186 * t226 + t146;
t151 = t203 * t156 - t200 * t157;
t148 = -qJDD(2) * pkin(3) - t213 * pkin(7) - t151;
t225 = qJD(2) * qJD(4);
t222 = t210 * t225;
t187 = t207 * qJDD(2) + t222;
t223 = t207 * t225;
t188 = t210 * qJDD(2) - t223;
t144 = (-t187 - t222) * pkin(8) + (-t188 + t223) * pkin(4) + t148;
t206 = sin(qJ(5));
t209 = cos(qJ(5));
t140 = -t206 * t143 + t209 * t144;
t227 = qJD(2) * t207;
t183 = t209 * qJD(4) - t206 * t227;
t166 = t183 * qJD(5) + t206 * qJDD(4) + t209 * t187;
t184 = t206 * qJD(4) + t209 * t227;
t167 = -t183 * mrSges(6,1) + t184 * mrSges(6,2);
t196 = qJD(5) - t226;
t170 = -t196 * mrSges(6,2) + t183 * mrSges(6,3);
t181 = qJDD(5) - t188;
t137 = m(6) * t140 + t181 * mrSges(6,1) - t166 * mrSges(6,3) - t184 * t167 + t196 * t170;
t141 = t209 * t143 + t206 * t144;
t165 = -t184 * qJD(5) + t209 * qJDD(4) - t206 * t187;
t171 = t196 * mrSges(6,1) - t184 * mrSges(6,3);
t138 = m(6) * t141 - t181 * mrSges(6,2) + t165 * mrSges(6,3) + t183 * t167 - t196 * t171;
t131 = -t206 * t137 + t209 * t138;
t228 = t210 * t172;
t142 = -qJDD(4) * pkin(4) - t212 * pkin(8) - t228 + (qJD(2) * t186 + t149) * t207;
t160 = Ifges(6,5) * t184 + Ifges(6,6) * t183 + Ifges(6,3) * t196;
t162 = Ifges(6,1) * t184 + Ifges(6,4) * t183 + Ifges(6,5) * t196;
t132 = -mrSges(6,1) * t142 + mrSges(6,3) * t141 + Ifges(6,4) * t166 + Ifges(6,2) * t165 + Ifges(6,6) * t181 - t184 * t160 + t196 * t162;
t161 = Ifges(6,4) * t184 + Ifges(6,2) * t183 + Ifges(6,6) * t196;
t133 = mrSges(6,2) * t142 - mrSges(6,3) * t140 + Ifges(6,1) * t166 + Ifges(6,4) * t165 + Ifges(6,5) * t181 + t183 * t160 - t196 * t161;
t139 = -m(6) * t142 + t165 * mrSges(6,1) - t166 * mrSges(6,2) + t183 * t170 - t184 * t171;
t145 = -t207 * t149 + t228;
t177 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t207 + Ifges(5,2) * t210) * qJD(2);
t178 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t207 + Ifges(5,4) * t210) * qJD(2);
t234 = mrSges(5,1) * t145 - mrSges(5,2) * t146 + Ifges(5,5) * t187 + Ifges(5,6) * t188 + Ifges(5,3) * qJDD(4) + pkin(4) * t139 + pkin(8) * t131 + t209 * t132 + t206 * t133 + (t207 * t177 - t210 * t178) * qJD(2);
t185 = (-mrSges(5,1) * t210 + mrSges(5,2) * t207) * qJD(2);
t192 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t227;
t128 = m(5) * t146 - qJDD(4) * mrSges(5,2) + t188 * mrSges(5,3) - qJD(4) * t192 + t185 * t226 + t131;
t193 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t226;
t135 = m(5) * t145 + qJDD(4) * mrSges(5,1) - t187 * mrSges(5,3) + qJD(4) * t193 - t185 * t227 + t139;
t220 = t210 * t128 - t207 * t135;
t120 = m(4) * t152 - t213 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t220;
t130 = t209 * t137 + t206 * t138;
t216 = -m(5) * t148 + t188 * mrSges(5,1) - t187 * mrSges(5,2) - t192 * t227 + t193 * t226 - t130;
t125 = m(4) * t151 + qJDD(2) * mrSges(4,1) - t213 * mrSges(4,2) + t216;
t113 = t200 * t120 + t203 * t125;
t111 = m(3) * t158 + qJDD(2) * mrSges(3,1) - t213 * mrSges(3,2) + t113;
t221 = t203 * t120 - t200 * t125;
t112 = m(3) * t159 - t213 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t221;
t106 = -t208 * t111 + t211 * t112;
t233 = pkin(6) * t106;
t123 = t207 * t128 + t210 * t135;
t224 = m(4) * t172 + t123;
t121 = m(3) * t173 + t224;
t102 = t111 * t229 + t112 * t230 - t202 * t121;
t176 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t207 + Ifges(5,6) * t210) * qJD(2);
t115 = mrSges(5,2) * t148 - mrSges(5,3) * t145 + Ifges(5,1) * t187 + Ifges(5,4) * t188 + Ifges(5,5) * qJDD(4) - pkin(8) * t130 - qJD(4) * t177 - t206 * t132 + t209 * t133 + t176 * t226;
t215 = mrSges(6,1) * t140 - mrSges(6,2) * t141 + Ifges(6,5) * t166 + Ifges(6,6) * t165 + Ifges(6,3) * t181 + t184 * t161 - t183 * t162;
t117 = -mrSges(5,1) * t148 + mrSges(5,3) * t146 + Ifges(5,4) * t187 + Ifges(5,2) * t188 + Ifges(5,6) * qJDD(4) - pkin(4) * t130 + qJD(4) * t178 - t176 * t227 - t215;
t103 = mrSges(4,2) * t172 - mrSges(4,3) * t151 + Ifges(4,5) * qJDD(2) - t213 * Ifges(4,6) - pkin(7) * t123 + t210 * t115 - t207 * t117;
t107 = -mrSges(4,1) * t172 + mrSges(4,3) * t152 + t213 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t123 - t234;
t94 = -mrSges(3,1) * t173 + mrSges(3,3) * t159 + t213 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t224 + qJ(3) * t221 + t200 * t103 + t203 * t107;
t96 = mrSges(3,2) * t173 - mrSges(3,3) * t158 + Ifges(3,5) * qJDD(2) - t213 * Ifges(3,6) - qJ(3) * t113 + t203 * t103 - t200 * t107;
t217 = mrSges(4,1) * t151 - mrSges(4,2) * t152 + Ifges(4,3) * qJDD(2) + pkin(3) * t216 + pkin(7) * t220 + t207 * t115 + t210 * t117;
t98 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + Ifges(3,3) * qJDD(2) + pkin(2) * t113 + t217;
t218 = mrSges(2,1) * t190 - mrSges(2,2) * t191 + pkin(1) * t102 + t202 * t233 + t205 * t98 + t94 * t231 + t96 * t232;
t104 = m(2) * t191 + t106;
t101 = t205 * t121 + (t111 * t211 + t112 * t208) * t202;
t99 = m(2) * t190 + t102;
t92 = mrSges(2,2) * t199 - mrSges(2,3) * t190 - t208 * t94 + t211 * t96 + (-t101 * t202 - t102 * t205) * pkin(6);
t91 = -mrSges(2,1) * t199 + mrSges(2,3) * t191 - pkin(1) * t101 - t202 * t98 + (t208 * t96 + t211 * t94 + t233) * t205;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t204 * t92 - t201 * t91 - qJ(1) * (t201 * t104 + t204 * t99), t92, t96, t103, t115, t133; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t201 * t92 + t204 * t91 + qJ(1) * (t204 * t104 - t201 * t99), t91, t94, t107, t117, t132; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t218, t218, t98, t217, t234, t215;];
m_new = t1;
