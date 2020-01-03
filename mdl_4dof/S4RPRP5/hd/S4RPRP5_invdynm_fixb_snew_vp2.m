% Calculate vector of cutting torques with Newton-Euler for
% S4RPRP5
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:48
% EndTime: 2019-12-31 16:44:50
% DurationCPUTime: 1.26s
% Computational Cost: add. (10916->213), mult. (25861->260), div. (0->0), fcn. (15800->6), ass. (0->89)
t194 = qJD(1) ^ 2;
t188 = sin(pkin(6));
t221 = qJD(1) * t188;
t189 = cos(pkin(6));
t230 = qJD(1) * t189;
t191 = sin(qJ(1));
t192 = cos(qJ(1));
t171 = -g(1) * t192 - g(2) * t191;
t164 = -pkin(1) * t194 + qJDD(1) * qJ(2) + t171;
t218 = qJD(1) * qJD(2);
t213 = -t189 * g(3) - 0.2e1 * t188 * t218;
t152 = -t188 * t164 + t213;
t153 = -g(3) * t188 + (t164 + 0.2e1 * t218) * t189;
t227 = pkin(2) * t189;
t125 = (-pkin(5) * qJDD(1) + t194 * t227 - t164) * t188 + t213;
t216 = qJDD(1) * t189;
t181 = t189 ^ 2;
t224 = t181 * t194;
t126 = -pkin(2) * t224 + pkin(5) * t216 + t153;
t190 = sin(qJ(3));
t228 = cos(qJ(3));
t119 = t228 * t125 - t190 * t126;
t120 = t190 * t125 + t228 * t126;
t214 = t189 * t228;
t162 = -qJD(1) * t214 + t190 * t221;
t202 = t228 * t188 + t189 * t190;
t163 = t202 * qJD(1);
t131 = Ifges(4,4) * t163 - Ifges(4,2) * t162 + Ifges(4,6) * qJD(3);
t133 = Ifges(4,1) * t163 - Ifges(4,4) * t162 + Ifges(4,5) * qJD(3);
t138 = mrSges(5,1) * t162 - mrSges(5,3) * t163;
t217 = qJDD(1) * t188;
t220 = qJD(3) * t163;
t150 = -qJDD(1) * t214 + t190 * t217 + t220;
t219 = t162 * qJD(3);
t151 = t202 * qJDD(1) - t219;
t137 = pkin(3) * t162 - qJ(4) * t163;
t193 = qJD(3) ^ 2;
t115 = -pkin(3) * t193 + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t137 * t162 + t120;
t117 = -qJDD(3) * pkin(3) - t193 * qJ(4) + t163 * t137 + qJDD(4) - t119;
t128 = Ifges(5,5) * t163 + Ifges(5,6) * qJD(3) + Ifges(5,3) * t162;
t132 = Ifges(5,1) * t163 + Ifges(5,4) * qJD(3) + Ifges(5,5) * t162;
t199 = mrSges(5,1) * t117 - mrSges(5,3) * t115 - Ifges(5,4) * t151 - Ifges(5,2) * qJDD(3) - Ifges(5,6) * t150 + t163 * t128 - t162 * t132;
t159 = -mrSges(5,2) * t162 + qJD(3) * mrSges(5,3);
t210 = -m(5) * t117 + qJDD(3) * mrSges(5,1) + qJD(3) * t159;
t158 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t163;
t215 = m(5) * t115 + qJDD(3) * mrSges(5,3) + qJD(3) * t158;
t196 = mrSges(4,2) * t120 - t162 * t133 - qJ(4) * (-t150 * mrSges(5,2) - t162 * t138 + t215) - pkin(3) * (-t151 * mrSges(5,2) - t163 * t138 + t210) - mrSges(4,1) * t119 - t163 * t131 + Ifges(4,6) * t150 - Ifges(4,5) * t151 - Ifges(4,3) * qJDD(3) + t199;
t206 = Ifges(3,4) * t188 + Ifges(3,2) * t189;
t207 = Ifges(3,1) * t188 + Ifges(3,4) * t189;
t157 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t163;
t222 = -mrSges(4,1) * t162 - mrSges(4,2) * t163 - t138;
t226 = -mrSges(4,3) - mrSges(5,2);
t104 = m(4) * t120 - qJDD(3) * mrSges(4,2) - qJD(3) * t157 + t226 * t150 + t222 * t162 + t215;
t156 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t162;
t105 = m(4) * t119 + qJDD(3) * mrSges(4,1) + qJD(3) * t156 + t226 * t151 + t222 * t163 + t210;
t98 = t190 * t104 + t228 * t105;
t229 = -mrSges(3,1) * t152 + mrSges(3,2) * t153 - pkin(2) * t98 - (t206 * t221 - t207 * t230) * qJD(1) + t196;
t225 = mrSges(3,2) * t188;
t130 = Ifges(5,4) * t163 + Ifges(5,2) * qJD(3) + Ifges(5,6) * t162;
t223 = -Ifges(4,5) * t163 + Ifges(4,6) * t162 - Ifges(4,3) * qJD(3) - t130;
t170 = t191 * g(1) - t192 * g(2);
t203 = mrSges(3,3) * qJDD(1) + t194 * (-mrSges(3,1) * t189 + t225);
t96 = m(3) * t152 - t203 * t188 + t98;
t211 = t228 * t104 - t190 * t105;
t97 = m(3) * t153 + t203 * t189 + t211;
t212 = -t188 * t96 + t189 * t97;
t209 = qJDD(2) - t170;
t180 = t188 ^ 2;
t149 = (-pkin(1) - t227) * qJDD(1) + (-qJ(2) + (-t180 - t181) * pkin(5)) * t194 + t209;
t112 = -0.2e1 * qJD(4) * t163 + (-t151 + t219) * qJ(4) + (t150 + t220) * pkin(3) + t149;
t208 = -mrSges(5,1) * t112 + mrSges(5,2) * t115;
t205 = Ifges(3,5) * t188 + Ifges(3,6) * t189;
t106 = m(5) * t112 + t150 * mrSges(5,1) - t151 * mrSges(5,3) - t163 * t158 + t162 * t159;
t201 = mrSges(5,2) * t117 - mrSges(5,3) * t112 + Ifges(5,1) * t151 + Ifges(5,4) * qJDD(3) + Ifges(5,5) * t150 + qJD(3) * t128;
t161 = -qJDD(1) * pkin(1) - t194 * qJ(2) + t209;
t198 = m(4) * t149 + t150 * mrSges(4,1) + mrSges(4,2) * t151 + t162 * t156 + t157 * t163 + t106;
t197 = -m(3) * t161 + mrSges(3,1) * t216 - t198 + (t180 * t194 + t224) * mrSges(3,3);
t166 = t205 * qJD(1);
t93 = -mrSges(4,1) * t149 + mrSges(4,3) * t120 - pkin(3) * t106 + t223 * t163 + (Ifges(4,4) - Ifges(5,5)) * t151 + (-Ifges(4,2) - Ifges(5,3)) * t150 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + (t132 + t133) * qJD(3) + t208;
t94 = mrSges(4,2) * t149 - mrSges(4,3) * t119 + Ifges(4,1) * t151 - Ifges(4,4) * t150 + Ifges(4,5) * qJDD(3) - qJ(4) * t106 - qJD(3) * t131 + t223 * t162 + t201;
t87 = -mrSges(3,1) * t161 + mrSges(3,3) * t153 - pkin(2) * t198 + pkin(5) * t211 + t206 * qJDD(1) - t166 * t221 + t190 * t94 + t228 * t93;
t89 = mrSges(3,2) * t161 - mrSges(3,3) * t152 - pkin(5) * t98 + t207 * qJDD(1) + t166 * t230 - t190 * t93 + t228 * t94;
t200 = -mrSges(2,2) * t171 + qJ(2) * t212 + t188 * t89 + t189 * t87 + pkin(1) * (-mrSges(3,2) * t217 + t197) + mrSges(2,1) * t170 + Ifges(2,3) * qJDD(1);
t99 = t197 + (mrSges(2,1) - t225) * qJDD(1) + m(2) * t170 - mrSges(2,2) * t194;
t92 = t188 * t97 + t189 * t96;
t90 = m(2) * t171 - mrSges(2,1) * t194 - qJDD(1) * mrSges(2,2) + t212;
t85 = mrSges(2,1) * g(3) + (Ifges(2,6) - t205) * qJDD(1) + t194 * Ifges(2,5) + mrSges(2,3) * t171 - pkin(1) * t92 + t229;
t84 = -mrSges(2,2) * g(3) - mrSges(2,3) * t170 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t194 - qJ(2) * t92 - t188 * t87 + t189 * t89;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t192 * t84 - t191 * t85 - pkin(4) * (t191 * t90 + t192 * t99), t84, t89, t94, -t130 * t162 + t201; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t191 * t84 + t192 * t85 + pkin(4) * (-t191 * t99 + t192 * t90), t85, t87, t93, -t199; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t200, t200, t205 * qJDD(1) - t229, -t196, Ifges(5,5) * t151 + Ifges(5,6) * qJDD(3) + Ifges(5,3) * t150 - qJD(3) * t132 + t130 * t163 - t208;];
m_new = t1;
