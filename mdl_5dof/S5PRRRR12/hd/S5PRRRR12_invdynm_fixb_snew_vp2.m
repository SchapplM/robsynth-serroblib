% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:40
% EndTime: 2024-09-28 18:07:50
% DurationCPUTime: 5.66s
% Computational Cost: add. (204370->193), mult. (277052->267), div. (0->0), fcn. (211667->14), ass. (0->103)
t203 = sin(pkin(11));
t206 = cos(pkin(11));
t190 = t203 * g(1) - t206 * g(2);
t191 = -t206 * g(1) - t203 * g(2);
t202 = -g(3) + qJDD(1);
t212 = sin(qJ(2));
t208 = cos(pkin(5));
t216 = cos(qJ(2));
t229 = t208 * t216;
t205 = sin(pkin(5));
t231 = t205 * t216;
t171 = t190 * t229 - t212 * t191 + t202 * t231;
t169 = qJDD(2) * pkin(2) + t171;
t230 = t208 * t212;
t232 = t205 * t212;
t172 = t190 * t230 + t216 * t191 + t202 * t232;
t217 = qJD(2) ^ 2;
t170 = -t217 * pkin(2) + t172;
t211 = sin(qJ(3));
t215 = cos(qJ(3));
t164 = t215 * t169 - t211 * t170;
t200 = qJDD(2) + qJDD(3);
t161 = t200 * pkin(3) + t164;
t165 = t211 * t169 + t215 * t170;
t201 = qJD(2) + qJD(3);
t199 = t201 ^ 2;
t162 = -t199 * pkin(3) + t165;
t210 = sin(qJ(4));
t214 = cos(qJ(4));
t157 = t210 * t161 + t214 * t162;
t198 = qJD(4) + t201;
t196 = t198 ^ 2;
t197 = qJDD(4) + t200;
t204 = sin(pkin(6));
t236 = pkin(10) * t204;
t154 = -t196 * pkin(4) + t197 * t236 + t157;
t209 = sin(qJ(5));
t213 = cos(qJ(5));
t156 = t214 * t161 - t210 * t162;
t153 = t197 * pkin(4) + t196 * t236 + t156;
t183 = -t205 * t190 + t208 * t202;
t207 = cos(pkin(6));
t221 = t153 * t207 + t183 * t204;
t148 = -t209 * t154 + t221 * t213;
t189 = t207 * t198 + qJD(5);
t233 = t204 * t213;
t226 = t198 * t233;
t177 = -t189 * mrSges(6,2) + mrSges(6,3) * t226;
t235 = t198 * t204;
t178 = (-mrSges(6,1) * t213 + mrSges(6,2) * t209) * t235;
t228 = qJD(5) * t198;
t179 = (t197 * t209 + t213 * t228) * t204;
t188 = t207 * t197 + qJDD(5);
t234 = t204 * t209;
t227 = t198 * t234;
t146 = m(6) * t148 + t188 * mrSges(6,1) - t179 * mrSges(6,3) + t189 * t177 - t178 * t227;
t149 = t213 * t154 + t221 * t209;
t176 = t189 * mrSges(6,1) - mrSges(6,3) * t227;
t180 = (t197 * t213 - t209 * t228) * t204;
t147 = m(6) * t149 - t188 * mrSges(6,2) + t180 * mrSges(6,3) - t189 * t176 + t178 * t226;
t152 = -t204 * t153 + t207 * t183;
t151 = m(6) * t152 - t180 * mrSges(6,1) + t179 * mrSges(6,2) + (t176 * t209 - t177 * t213) * t235;
t129 = -t204 * t151 + (t146 * t213 + t147 * t209) * t207;
t125 = m(5) * t156 + t197 * mrSges(5,1) - t196 * mrSges(5,2) + t129;
t134 = -t209 * t146 + t213 * t147;
t132 = m(5) * t157 - t196 * mrSges(5,1) - t197 * mrSges(5,2) + t134;
t123 = t214 * t125 + t210 * t132;
t120 = m(4) * t164 + t200 * mrSges(4,1) - t199 * mrSges(4,2) + t123;
t224 = -t210 * t125 + t214 * t132;
t121 = m(4) * t165 - t199 * mrSges(4,1) - t200 * mrSges(4,2) + t224;
t114 = t215 * t120 + t211 * t121;
t112 = m(3) * t171 + qJDD(2) * mrSges(3,1) - t217 * mrSges(3,2) + t114;
t225 = -t211 * t120 + t215 * t121;
t113 = m(3) * t172 - t217 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t225;
t108 = -t212 * t112 + t216 * t113;
t237 = pkin(7) * t108;
t128 = t146 * t233 + t147 * t234 + t207 * t151;
t223 = m(5) * t183 + t128;
t222 = m(4) * t183 + t223;
t126 = m(3) * t183 + t222;
t104 = t112 * t229 + t113 * t230 - t205 * t126;
t174 = Ifges(6,6) * t189 + (Ifges(6,4) * t209 + Ifges(6,2) * t213) * t235;
t175 = Ifges(6,5) * t189 + (Ifges(6,1) * t209 + Ifges(6,4) * t213) * t235;
t136 = mrSges(6,1) * t148 - mrSges(6,2) * t149 + Ifges(6,5) * t179 + Ifges(6,6) * t180 + Ifges(6,3) * t188 + (t174 * t209 - t175 * t213) * t235;
t173 = Ifges(6,3) * t189 + (Ifges(6,5) * t209 + Ifges(6,6) * t213) * t235;
t139 = -mrSges(6,1) * t152 + mrSges(6,3) * t149 + Ifges(6,4) * t179 + Ifges(6,2) * t180 + Ifges(6,6) * t188 - t173 * t227 + t189 * t175;
t140 = mrSges(6,2) * t152 - mrSges(6,3) * t148 + Ifges(6,1) * t179 + Ifges(6,4) * t180 + Ifges(6,5) * t188 + t173 * t226 - t189 * t174;
t115 = -mrSges(5,1) * t183 + mrSges(5,3) * t157 + t196 * Ifges(5,5) + Ifges(5,6) * t197 - pkin(4) * t128 - t204 * t136 + (pkin(10) * t134 + t139 * t213 + t140 * t209) * t207;
t116 = mrSges(5,2) * t183 - mrSges(5,3) * t156 + Ifges(5,5) * t197 - t196 * Ifges(5,6) - t209 * t139 + t213 * t140 + (-t128 * t204 - t129 * t207) * pkin(10);
t100 = -mrSges(4,1) * t183 + mrSges(4,3) * t165 + t199 * Ifges(4,5) + Ifges(4,6) * t200 - pkin(3) * t223 + pkin(9) * t224 + t214 * t115 + t210 * t116;
t105 = mrSges(4,2) * t183 - mrSges(4,3) * t164 + Ifges(4,5) * t200 - t199 * Ifges(4,6) - pkin(9) * t123 - t210 * t115 + t214 * t116;
t95 = -mrSges(3,1) * t183 + mrSges(3,3) * t172 + t217 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t222 + pkin(8) * t225 + t215 * t100 + t211 * t105;
t97 = mrSges(3,2) * t183 - mrSges(3,3) * t171 + Ifges(3,5) * qJDD(2) - t217 * Ifges(3,6) - pkin(8) * t114 - t211 * t100 + t215 * t105;
t219 = mrSges(5,1) * t156 - mrSges(5,2) * t157 + Ifges(5,3) * t197 + pkin(4) * t129 + t134 * t236 + t207 * t136 + t139 * t233 + t140 * t234;
t218 = mrSges(4,1) * t164 - mrSges(4,2) * t165 + Ifges(4,3) * t200 + pkin(3) * t123 + t219;
t99 = mrSges(3,1) * t171 - mrSges(3,2) * t172 + Ifges(3,3) * qJDD(2) + pkin(2) * t114 + t218;
t220 = mrSges(2,1) * t190 - mrSges(2,2) * t191 + pkin(1) * t104 + t205 * t237 + t208 * t99 + t95 * t231 + t97 * t232;
t106 = m(2) * t191 + t108;
t103 = t208 * t126 + (t112 * t216 + t113 * t212) * t205;
t101 = m(2) * t190 + t104;
t93 = mrSges(2,2) * t202 - mrSges(2,3) * t190 - t212 * t95 + t216 * t97 + (-t103 * t205 - t104 * t208) * pkin(7);
t92 = -mrSges(2,1) * t202 + mrSges(2,3) * t191 - pkin(1) * t103 - t205 * t99 + (t212 * t97 + t216 * t95 + t237) * t208;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t206 * t93 - t203 * t92 - qJ(1) * (t206 * t101 + t203 * t106), t93, t97, t105, t116, t140; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t203 * t93 + t206 * t92 + qJ(1) * (-t203 * t101 + t206 * t106), t92, t95, t100, t115, t139; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t220, t220, t99, t218, t219, t136;];
m_new = t1;
