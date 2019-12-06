% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:16
% EndTime: 2019-12-05 18:18:20
% DurationCPUTime: 3.71s
% Computational Cost: add. (65636->204), mult. (89300->259), div. (0->0), fcn. (51705->10), ass. (0->97)
t195 = qJD(1) + qJD(2);
t191 = t195 ^ 2;
t205 = sin(qJ(1));
t208 = cos(qJ(1));
t179 = t208 * g(2) + t205 * g(3);
t173 = qJDD(1) * pkin(1) + t179;
t178 = t205 * g(2) - g(3) * t208;
t209 = qJD(1) ^ 2;
t174 = -pkin(1) * t209 + t178;
t204 = sin(qJ(2));
t207 = cos(qJ(2));
t160 = t207 * t173 - t174 * t204;
t192 = qJDD(1) + qJDD(2);
t157 = pkin(2) * t192 + t160;
t161 = t204 * t173 + t207 * t174;
t158 = -pkin(2) * t191 + t161;
t200 = sin(pkin(8));
t202 = cos(pkin(8));
t142 = t200 * t157 + t202 * t158;
t139 = -pkin(3) * t191 + qJ(4) * t192 + t142;
t199 = sin(pkin(9));
t198 = -g(1) + qJDD(3);
t201 = cos(pkin(9));
t230 = qJD(4) * t195;
t231 = t201 * t198 - 0.2e1 * t199 * t230;
t236 = pkin(4) * t201;
t132 = (-pkin(7) * t192 + t191 * t236 - t139) * t199 + t231;
t136 = t199 * t198 + (t139 + 0.2e1 * t230) * t201;
t233 = t192 * t201;
t194 = t201 ^ 2;
t234 = t191 * t194;
t133 = -pkin(4) * t234 + pkin(7) * t233 + t136;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t130 = t132 * t206 - t133 * t203;
t217 = -t199 * t203 + t201 * t206;
t164 = t217 * t195;
t218 = t199 * t206 + t201 * t203;
t165 = t218 * t195;
t151 = -mrSges(6,1) * t164 + mrSges(6,2) * t165;
t153 = t164 * qJD(5) + t192 * t218;
t162 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t164;
t127 = m(6) * t130 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t153 + qJD(5) * t162 - t151 * t165;
t131 = t132 * t203 + t133 * t206;
t152 = -t165 * qJD(5) + t192 * t217;
t163 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t165;
t128 = m(6) * t131 - qJDD(5) * mrSges(6,2) + t152 * mrSges(6,3) - qJD(5) * t163 + t151 * t164;
t118 = t127 * t206 + t128 * t203;
t135 = -t139 * t199 + t231;
t145 = Ifges(6,4) * t165 + Ifges(6,2) * t164 + Ifges(6,6) * qJD(5);
t146 = Ifges(6,1) * t165 + Ifges(6,4) * t164 + Ifges(6,5) * qJD(5);
t214 = -mrSges(6,1) * t130 + mrSges(6,2) * t131 - Ifges(6,5) * t153 - Ifges(6,6) * t152 - Ifges(6,3) * qJDD(5) - t165 * t145 + t164 * t146;
t223 = Ifges(5,4) * t199 + Ifges(5,2) * t201;
t224 = Ifges(5,1) * t199 + Ifges(5,4) * t201;
t237 = -mrSges(5,1) * t135 + mrSges(5,2) * t136 - pkin(4) * t118 - (t199 * t223 - t201 * t224) * t191 + t214;
t235 = mrSges(5,2) * t199;
t220 = mrSges(5,3) * t192 + (-mrSges(5,1) * t201 + t235) * t191;
t116 = m(5) * t135 - t199 * t220 + t118;
t225 = -t203 * t127 + t128 * t206;
t117 = m(5) * t136 + t201 * t220 + t225;
t226 = -t116 * t199 + t117 * t201;
t110 = m(4) * t142 - mrSges(4,1) * t191 - mrSges(4,2) * t192 + t226;
t141 = t202 * t157 - t200 * t158;
t221 = qJDD(4) - t141;
t138 = -t192 * pkin(3) - t191 * qJ(4) + t221;
t193 = t199 ^ 2;
t134 = (-pkin(3) - t236) * t192 + (-qJ(4) + (-t193 - t194) * pkin(7)) * t191 + t221;
t215 = m(6) * t134 - t152 * mrSges(6,1) + t153 * mrSges(6,2) - t164 * t162 + t165 * t163;
t212 = -m(5) * t138 + mrSges(5,1) * t233 - t215 + (t191 * t193 + t234) * mrSges(5,3);
t122 = m(4) * t141 - t191 * mrSges(4,2) + (mrSges(4,1) - t235) * t192 + t212;
t105 = t110 * t200 + t122 * t202;
t100 = m(3) * t160 + mrSges(3,1) * t192 - mrSges(3,2) * t191 + t105;
t227 = t110 * t202 - t122 * t200;
t101 = m(3) * t161 - mrSges(3,1) * t191 - mrSges(3,2) * t192 + t227;
t95 = t100 * t207 + t101 * t204;
t222 = Ifges(5,5) * t199 + Ifges(5,6) * t201;
t232 = t191 * t222;
t112 = t116 * t201 + t117 * t199;
t229 = m(4) * t198 + t112;
t228 = -t100 * t204 + t101 * t207;
t144 = Ifges(6,5) * t165 + Ifges(6,6) * t164 + Ifges(6,3) * qJD(5);
t119 = -mrSges(6,1) * t134 + mrSges(6,3) * t131 + Ifges(6,4) * t153 + Ifges(6,2) * t152 + Ifges(6,6) * qJDD(5) + qJD(5) * t146 - t144 * t165;
t120 = mrSges(6,2) * t134 - mrSges(6,3) * t130 + Ifges(6,1) * t153 + Ifges(6,4) * t152 + Ifges(6,5) * qJDD(5) - qJD(5) * t145 + t144 * t164;
t103 = -mrSges(5,1) * t138 + mrSges(5,3) * t136 - pkin(4) * t215 + pkin(7) * t225 + t206 * t119 + t203 * t120 + t192 * t223 - t199 * t232;
t107 = mrSges(5,2) * t138 - mrSges(5,3) * t135 - pkin(7) * t118 - t203 * t119 + t206 * t120 + t192 * t224 + t201 * t232;
t216 = -mrSges(4,2) * t142 + qJ(4) * t226 + t201 * t103 + t199 * t107 + pkin(3) * (-t192 * t235 + t212) + mrSges(4,1) * t141 + Ifges(4,3) * t192;
t213 = mrSges(3,1) * t160 - mrSges(3,2) * t161 + Ifges(3,3) * t192 + pkin(2) * t105 + t216;
t210 = mrSges(2,1) * t179 - mrSges(2,2) * t178 + Ifges(2,3) * qJDD(1) + pkin(1) * t95 + t213;
t96 = (Ifges(4,6) - t222) * t192 - mrSges(4,1) * t198 + t191 * Ifges(4,5) + mrSges(4,3) * t142 - pkin(3) * t112 + t237;
t93 = m(2) * t179 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t209 + t95;
t92 = m(2) * t178 - mrSges(2,1) * t209 - qJDD(1) * mrSges(2,2) + t228;
t91 = mrSges(4,2) * t198 - mrSges(4,3) * t141 + Ifges(4,5) * t192 - Ifges(4,6) * t191 - qJ(4) * t112 - t103 * t199 + t107 * t201;
t90 = -mrSges(3,2) * g(1) - mrSges(3,3) * t160 + Ifges(3,5) * t192 - Ifges(3,6) * t191 - qJ(3) * t105 - t200 * t96 + t202 * t91;
t89 = mrSges(3,1) * g(1) + mrSges(3,3) * t161 + t191 * Ifges(3,5) + Ifges(3,6) * t192 - pkin(2) * t229 + qJ(3) * t227 + t200 * t91 + t202 * t96;
t88 = -mrSges(2,2) * g(1) - mrSges(2,3) * t179 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t209 - pkin(6) * t95 - t204 * t89 + t207 * t90;
t87 = Ifges(2,6) * qJDD(1) + t209 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t178 + t204 * t90 + t207 * t89 - pkin(1) * (-m(3) * g(1) + t229) + pkin(6) * t228;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t210, t88, t90, t91, t107, t120; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t205 * t88 - t208 * t87 - pkin(5) * (-t205 * t93 + t208 * t92), t87, t89, t96, t103, t119; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t208 * t88 - t205 * t87 + pkin(5) * (-t205 * t92 - t208 * t93), t210, t213, t216, t192 * t222 - t237, -t214;];
m_new = t1;
