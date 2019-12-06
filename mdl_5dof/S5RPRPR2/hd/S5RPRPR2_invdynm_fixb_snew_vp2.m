% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:26
% EndTime: 2019-12-05 17:49:30
% DurationCPUTime: 3.79s
% Computational Cost: add. (62114->204), mult. (89300->259), div. (0->0), fcn. (51705->10), ass. (0->97)
t195 = qJD(1) + qJD(3);
t191 = t195 ^ 2;
t206 = sin(qJ(1));
t209 = cos(qJ(1));
t180 = t209 * g(2) + t206 * g(3);
t174 = qJDD(1) * pkin(1) + t180;
t179 = t206 * g(2) - t209 * g(3);
t210 = qJD(1) ^ 2;
t175 = -t210 * pkin(1) + t179;
t201 = sin(pkin(8));
t203 = cos(pkin(8));
t161 = t203 * t174 - t201 * t175;
t158 = qJDD(1) * pkin(2) + t161;
t162 = t201 * t174 + t203 * t175;
t159 = -t210 * pkin(2) + t162;
t205 = sin(qJ(3));
t208 = cos(qJ(3));
t143 = t205 * t158 + t208 * t159;
t192 = qJDD(1) + qJDD(3);
t140 = -t191 * pkin(3) + t192 * qJ(4) + t143;
t200 = sin(pkin(9));
t199 = -g(1) + qJDD(2);
t202 = cos(pkin(9));
t231 = qJD(4) * t195;
t232 = t202 * t199 - 0.2e1 * t200 * t231;
t237 = pkin(4) * t202;
t133 = (-pkin(7) * t192 + t191 * t237 - t140) * t200 + t232;
t137 = t200 * t199 + (t140 + 0.2e1 * t231) * t202;
t234 = t192 * t202;
t194 = t202 ^ 2;
t235 = t191 * t194;
t134 = -pkin(4) * t235 + pkin(7) * t234 + t137;
t204 = sin(qJ(5));
t207 = cos(qJ(5));
t131 = t207 * t133 - t204 * t134;
t218 = -t200 * t204 + t202 * t207;
t165 = t218 * t195;
t219 = t200 * t207 + t202 * t204;
t166 = t219 * t195;
t152 = -t165 * mrSges(6,1) + t166 * mrSges(6,2);
t154 = t165 * qJD(5) + t219 * t192;
t163 = -qJD(5) * mrSges(6,2) + t165 * mrSges(6,3);
t128 = m(6) * t131 + qJDD(5) * mrSges(6,1) - t154 * mrSges(6,3) + qJD(5) * t163 - t166 * t152;
t132 = t204 * t133 + t207 * t134;
t153 = -t166 * qJD(5) + t218 * t192;
t164 = qJD(5) * mrSges(6,1) - t166 * mrSges(6,3);
t129 = m(6) * t132 - qJDD(5) * mrSges(6,2) + t153 * mrSges(6,3) - qJD(5) * t164 + t165 * t152;
t119 = t207 * t128 + t204 * t129;
t136 = -t200 * t140 + t232;
t146 = Ifges(6,4) * t166 + Ifges(6,2) * t165 + Ifges(6,6) * qJD(5);
t147 = Ifges(6,1) * t166 + Ifges(6,4) * t165 + Ifges(6,5) * qJD(5);
t215 = -mrSges(6,1) * t131 + mrSges(6,2) * t132 - Ifges(6,5) * t154 - Ifges(6,6) * t153 - Ifges(6,3) * qJDD(5) - t166 * t146 + t165 * t147;
t224 = Ifges(5,4) * t200 + Ifges(5,2) * t202;
t225 = Ifges(5,1) * t200 + Ifges(5,4) * t202;
t238 = -mrSges(5,1) * t136 + mrSges(5,2) * t137 - pkin(4) * t119 - (t200 * t224 - t202 * t225) * t191 + t215;
t236 = mrSges(5,2) * t200;
t221 = mrSges(5,3) * t192 + (-mrSges(5,1) * t202 + t236) * t191;
t117 = m(5) * t136 - t221 * t200 + t119;
t226 = -t204 * t128 + t207 * t129;
t118 = m(5) * t137 + t221 * t202 + t226;
t227 = -t200 * t117 + t202 * t118;
t111 = m(4) * t143 - t191 * mrSges(4,1) - t192 * mrSges(4,2) + t227;
t142 = t208 * t158 - t205 * t159;
t222 = qJDD(4) - t142;
t139 = -t192 * pkin(3) - t191 * qJ(4) + t222;
t193 = t200 ^ 2;
t135 = (-pkin(3) - t237) * t192 + (-qJ(4) + (-t193 - t194) * pkin(7)) * t191 + t222;
t216 = m(6) * t135 - t153 * mrSges(6,1) + t154 * mrSges(6,2) - t165 * t163 + t166 * t164;
t213 = -m(5) * t139 + mrSges(5,1) * t234 - t216 + (t191 * t193 + t235) * mrSges(5,3);
t123 = m(4) * t142 - t191 * mrSges(4,2) + (mrSges(4,1) - t236) * t192 + t213;
t106 = t205 * t111 + t208 * t123;
t101 = m(3) * t161 + qJDD(1) * mrSges(3,1) - t210 * mrSges(3,2) + t106;
t228 = t208 * t111 - t205 * t123;
t102 = m(3) * t162 - t210 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t228;
t96 = t203 * t101 + t201 * t102;
t223 = Ifges(5,5) * t200 + Ifges(5,6) * t202;
t233 = t191 * t223;
t113 = t202 * t117 + t200 * t118;
t230 = m(4) * t199 + t113;
t229 = -t201 * t101 + t203 * t102;
t145 = Ifges(6,5) * t166 + Ifges(6,6) * t165 + Ifges(6,3) * qJD(5);
t120 = -mrSges(6,1) * t135 + mrSges(6,3) * t132 + Ifges(6,4) * t154 + Ifges(6,2) * t153 + Ifges(6,6) * qJDD(5) + qJD(5) * t147 - t166 * t145;
t121 = mrSges(6,2) * t135 - mrSges(6,3) * t131 + Ifges(6,1) * t154 + Ifges(6,4) * t153 + Ifges(6,5) * qJDD(5) - qJD(5) * t146 + t165 * t145;
t104 = -mrSges(5,1) * t139 + mrSges(5,3) * t137 - pkin(4) * t216 + pkin(7) * t226 + t207 * t120 + t204 * t121 + t224 * t192 - t200 * t233;
t108 = mrSges(5,2) * t139 - mrSges(5,3) * t136 - pkin(7) * t119 - t204 * t120 + t207 * t121 + t225 * t192 + t202 * t233;
t217 = -mrSges(4,2) * t143 + qJ(4) * t227 + t202 * t104 + t200 * t108 + pkin(3) * (-t192 * t236 + t213) + mrSges(4,1) * t142 + Ifges(4,3) * t192;
t214 = mrSges(3,1) * t161 - mrSges(3,2) * t162 + Ifges(3,3) * qJDD(1) + pkin(2) * t106 + t217;
t211 = mrSges(2,1) * t180 - mrSges(2,2) * t179 + Ifges(2,3) * qJDD(1) + pkin(1) * t96 + t214;
t97 = (Ifges(4,6) - t223) * t192 + t191 * Ifges(4,5) - mrSges(4,1) * t199 + mrSges(4,3) * t143 - pkin(3) * t113 + t238;
t94 = m(2) * t180 + qJDD(1) * mrSges(2,1) - t210 * mrSges(2,2) + t96;
t93 = m(2) * t179 - t210 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t229;
t92 = mrSges(4,2) * t199 - mrSges(4,3) * t142 + Ifges(4,5) * t192 - t191 * Ifges(4,6) - qJ(4) * t113 - t200 * t104 + t202 * t108;
t91 = mrSges(3,2) * t199 - mrSges(3,3) * t161 + Ifges(3,5) * qJDD(1) - t210 * Ifges(3,6) - pkin(6) * t106 - t205 * t97 + t208 * t92;
t90 = -mrSges(3,1) * t199 + mrSges(3,3) * t162 + t210 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t230 + pkin(6) * t228 + t205 * t92 + t208 * t97;
t89 = -mrSges(2,2) * g(1) - mrSges(2,3) * t180 + Ifges(2,5) * qJDD(1) - t210 * Ifges(2,6) - qJ(2) * t96 - t201 * t90 + t203 * t91;
t88 = Ifges(2,6) * qJDD(1) + t210 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t179 + t201 * t91 + t203 * t90 - pkin(1) * (m(3) * t199 + t230) + qJ(2) * t229;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t211, t89, t91, t92, t108, t121; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t206 * t89 - t209 * t88 - pkin(5) * (-t206 * t94 + t209 * t93), t88, t90, t97, t104, t120; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209 * t89 - t206 * t88 + pkin(5) * (-t206 * t93 - t209 * t94), t211, t214, t217, t223 * t192 - t238, -t215;];
m_new = t1;
