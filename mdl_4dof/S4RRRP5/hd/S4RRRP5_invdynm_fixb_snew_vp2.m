% Calculate vector of cutting torques with Newton-Euler for
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
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:40
% EndTime: 2019-12-31 17:16:43
% DurationCPUTime: 1.47s
% Computational Cost: add. (14324->234), mult. (29773->287), div. (0->0), fcn. (17540->6), ass. (0->86)
t194 = cos(qJ(2));
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t178 = -t195 * g(1) - t193 * g(2);
t196 = qJD(1) ^ 2;
t165 = -t196 * pkin(1) + qJDD(1) * pkin(5) + t178;
t192 = sin(qJ(2));
t216 = t192 * t165;
t152 = -t194 * g(3) - t216;
t153 = -t192 * g(3) + t194 * t165;
t160 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t192 + Ifges(3,2) * t194) * qJD(1);
t161 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t192 + Ifges(3,4) * t194) * qJD(1);
t211 = qJD(1) * qJD(2);
t171 = t192 * qJDD(1) + t194 * t211;
t172 = t194 * qJDD(1) - t192 * t211;
t218 = pkin(2) * t196;
t123 = qJDD(2) * pkin(2) - t171 * pkin(6) - t216 + (pkin(6) * t211 + t192 * t218 - g(3)) * t194;
t213 = qJD(1) * t192;
t176 = qJD(2) * pkin(2) - pkin(6) * t213;
t190 = t194 ^ 2;
t124 = t172 * pkin(6) - qJD(2) * t176 - t190 * t218 + t153;
t191 = sin(qJ(3));
t219 = cos(qJ(3));
t117 = t219 * t123 - t191 * t124;
t118 = t191 * t123 + t219 * t124;
t163 = (t191 * t194 + t219 * t192) * qJD(1);
t134 = t163 * qJD(3) + t191 * t171 - t219 * t172;
t212 = qJD(1) * t194;
t162 = t191 * t213 - t219 * t212;
t135 = -t162 * qJD(3) + t219 * t171 + t191 * t172;
t188 = qJD(2) + qJD(3);
t141 = Ifges(4,4) * t163 - Ifges(4,2) * t162 + Ifges(4,6) * t188;
t143 = Ifges(4,1) * t163 - Ifges(4,4) * t162 + Ifges(4,5) * t188;
t148 = t162 * mrSges(5,1) - t163 * mrSges(5,3);
t187 = qJDD(2) + qJDD(3);
t147 = t162 * pkin(3) - t163 * qJ(4);
t186 = t188 ^ 2;
t113 = -t186 * pkin(3) + t187 * qJ(4) + 0.2e1 * qJD(4) * t188 - t162 * t147 + t118;
t115 = -t187 * pkin(3) - t186 * qJ(4) + t163 * t147 + qJDD(4) - t117;
t138 = Ifges(5,5) * t163 + Ifges(5,6) * t188 + Ifges(5,3) * t162;
t142 = Ifges(5,1) * t163 + Ifges(5,4) * t188 + Ifges(5,5) * t162;
t201 = mrSges(5,1) * t115 - mrSges(5,3) * t113 - Ifges(5,4) * t135 - Ifges(5,2) * t187 - Ifges(5,6) * t134 + t163 * t138 - t162 * t142;
t157 = -t162 * mrSges(5,2) + t188 * mrSges(5,3);
t207 = -m(5) * t115 + t187 * mrSges(5,1) + t188 * t157;
t156 = -t188 * mrSges(5,1) + t163 * mrSges(5,2);
t210 = m(5) * t113 + t187 * mrSges(5,3) + t188 * t156;
t199 = mrSges(4,2) * t118 - t162 * t143 - qJ(4) * (-t134 * mrSges(5,2) - t162 * t148 + t210) - pkin(3) * (-t135 * mrSges(5,2) - t163 * t148 + t207) - mrSges(4,1) * t117 - t163 * t141 + Ifges(4,6) * t134 - Ifges(4,5) * t135 - Ifges(4,3) * t187 + t201;
t155 = t188 * mrSges(4,1) - t163 * mrSges(4,3);
t214 = -t162 * mrSges(4,1) - t163 * mrSges(4,2) - t148;
t217 = -mrSges(4,3) - mrSges(5,2);
t101 = m(4) * t118 - t187 * mrSges(4,2) + t217 * t134 - t188 * t155 + t214 * t162 + t210;
t154 = -t188 * mrSges(4,2) - t162 * mrSges(4,3);
t103 = m(4) * t117 + t187 * mrSges(4,1) + t217 * t135 + t188 * t154 + t214 * t163 + t207;
t96 = t191 * t101 + t219 * t103;
t220 = mrSges(3,1) * t152 - mrSges(3,2) * t153 + Ifges(3,5) * t171 + Ifges(3,6) * t172 + Ifges(3,3) * qJDD(2) + pkin(2) * t96 + (t192 * t160 - t194 * t161) * qJD(1) - t199;
t140 = Ifges(5,4) * t163 + Ifges(5,2) * t188 + Ifges(5,6) * t162;
t215 = -Ifges(4,5) * t163 + Ifges(4,6) * t162 - Ifges(4,3) * t188 - t140;
t177 = t193 * g(1) - t195 * g(2);
t170 = (-mrSges(3,1) * t194 + mrSges(3,2) * t192) * qJD(1);
t175 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t212;
t94 = m(3) * t152 + qJDD(2) * mrSges(3,1) - t171 * mrSges(3,3) + qJD(2) * t175 - t170 * t213 + t96;
t174 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t213;
t208 = t219 * t101 - t191 * t103;
t95 = m(3) * t153 - qJDD(2) * mrSges(3,2) + t172 * mrSges(3,3) - qJD(2) * t174 + t170 * t212 + t208;
t209 = -t192 * t94 + t194 * t95;
t204 = -qJDD(1) * pkin(1) - t177;
t136 = -t172 * pkin(2) + t176 * t213 + (-pkin(6) * t190 - pkin(5)) * t196 + t204;
t110 = -0.2e1 * qJD(4) * t163 + (t162 * t188 - t135) * qJ(4) + (t163 * t188 + t134) * pkin(3) + t136;
t206 = -mrSges(5,1) * t110 + mrSges(5,2) * t113;
t104 = m(5) * t110 + t134 * mrSges(5,1) - t135 * mrSges(5,3) - t163 * t156 + t162 * t157;
t203 = mrSges(5,2) * t115 - mrSges(5,3) * t110 + Ifges(5,1) * t135 + Ifges(5,4) * t187 + Ifges(5,5) * t134 + t188 * t138;
t164 = -t196 * pkin(5) + t204;
t200 = m(4) * t136 + t134 * mrSges(4,1) + t135 * mrSges(4,2) + t162 * t154 + t163 * t155 + t104;
t198 = -m(3) * t164 + t172 * mrSges(3,1) - t171 * mrSges(3,2) - t174 * t213 + t175 * t212 - t200;
t159 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t192 + Ifges(3,6) * t194) * qJD(1);
t91 = -mrSges(4,1) * t136 + mrSges(4,3) * t118 - pkin(3) * t104 + (t142 + t143) * t188 + (Ifges(4,6) - Ifges(5,6)) * t187 + t215 * t163 + (Ifges(4,4) - Ifges(5,5)) * t135 + (-Ifges(4,2) - Ifges(5,3)) * t134 + t206;
t92 = mrSges(4,2) * t136 - mrSges(4,3) * t117 + Ifges(4,1) * t135 - Ifges(4,4) * t134 + Ifges(4,5) * t187 - qJ(4) * t104 - t188 * t141 + t215 * t162 + t203;
t85 = -mrSges(3,1) * t164 + mrSges(3,3) * t153 + Ifges(3,4) * t171 + Ifges(3,2) * t172 + Ifges(3,6) * qJDD(2) - pkin(2) * t200 + pkin(6) * t208 + qJD(2) * t161 - t159 * t213 + t191 * t92 + t219 * t91;
t87 = mrSges(3,2) * t164 - mrSges(3,3) * t152 + Ifges(3,1) * t171 + Ifges(3,4) * t172 + Ifges(3,5) * qJDD(2) - pkin(6) * t96 - qJD(2) * t160 + t159 * t212 - t191 * t91 + t219 * t92;
t202 = mrSges(2,1) * t177 - mrSges(2,2) * t178 + Ifges(2,3) * qJDD(1) + pkin(1) * t198 + pkin(5) * t209 + t192 * t87 + t194 * t85;
t97 = m(2) * t177 + qJDD(1) * mrSges(2,1) - t196 * mrSges(2,2) + t198;
t90 = t192 * t95 + t194 * t94;
t88 = m(2) * t178 - t196 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t209;
t83 = mrSges(2,1) * g(3) + mrSges(2,3) * t178 + t196 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t90 - t220;
t82 = -mrSges(2,2) * g(3) - mrSges(2,3) * t177 + Ifges(2,5) * qJDD(1) - t196 * Ifges(2,6) - pkin(5) * t90 - t192 * t85 + t194 * t87;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t195 * t82 - t193 * t83 - pkin(4) * (t193 * t88 + t195 * t97), t82, t87, t92, -t162 * t140 + t203; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t193 * t82 + t195 * t83 + pkin(4) * (-t193 * t97 + t195 * t88), t83, t85, t91, -t201; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t202, t202, t220, -t199, Ifges(5,5) * t135 + Ifges(5,6) * t187 + Ifges(5,3) * t134 + t163 * t140 - t188 * t142 - t206;];
m_new = t1;
