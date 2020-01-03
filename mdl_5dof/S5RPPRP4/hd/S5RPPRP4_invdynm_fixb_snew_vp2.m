% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:05
% EndTime: 2019-12-31 17:52:07
% DurationCPUTime: 1.34s
% Computational Cost: add. (12864->221), mult. (22516->263), div. (0->0), fcn. (8150->6), ass. (0->87)
t190 = sin(qJ(4));
t192 = cos(qJ(4));
t214 = qJD(1) * qJD(4);
t212 = t192 * t214;
t157 = -qJDD(1) * t190 - t212;
t195 = qJD(1) ^ 2;
t191 = sin(qJ(1));
t193 = cos(qJ(1));
t167 = -t193 * g(1) - t191 * g(2);
t206 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t167;
t223 = -pkin(1) - pkin(2);
t127 = t195 * t223 + t206;
t166 = t191 * g(1) - t193 * g(2);
t205 = -t195 * qJ(2) + qJDD(2) - t166;
t130 = qJDD(1) * t223 + t205;
t188 = sin(pkin(7));
t189 = cos(pkin(7));
t123 = t127 * t189 + t130 * t188;
t120 = -pkin(3) * t195 - qJDD(1) * pkin(6) + t123;
t184 = g(3) + qJDD(3);
t169 = t192 * t184;
t213 = qJD(1) * qJD(5);
t222 = pkin(4) * t195;
t112 = qJDD(4) * pkin(4) + t169 + (-t157 - t212) * qJ(5) + (t192 * t222 - t120 + (2 * t213)) * t190;
t155 = (mrSges(6,1) * t192 - mrSges(6,2) * t190) * qJD(1);
t215 = qJD(1) * t192;
t164 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t215;
t216 = qJD(1) * t190;
t211 = m(6) * t112 + qJDD(4) * mrSges(6,1) + qJD(4) * t164 + t155 * t216;
t107 = -mrSges(6,3) * t157 + t211;
t116 = -t190 * t120 + t169;
t117 = t120 * t192 + t190 * t184;
t143 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t190 - Ifges(5,4) * t192) * qJD(1);
t158 = -qJDD(1) * t192 + t190 * t214;
t161 = qJD(4) * pkin(4) + qJ(5) * t216;
t183 = t192 ^ 2;
t113 = qJ(5) * t158 - qJD(4) * t161 - t183 * t222 - 0.2e1 * t192 * t213 + t117;
t142 = Ifges(6,5) * qJD(4) + (-Ifges(6,1) * t190 - Ifges(6,4) * t192) * qJD(1);
t204 = -mrSges(6,1) * t112 + mrSges(6,2) * t113 - Ifges(6,5) * t157 - Ifges(6,6) * t158 - Ifges(6,3) * qJDD(4) - t142 * t215;
t140 = Ifges(6,6) * qJD(4) + (-Ifges(6,4) * t190 - Ifges(6,2) * t192) * qJD(1);
t217 = t140 + Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t190 - Ifges(5,2) * t192) * qJD(1);
t225 = mrSges(5,1) * t116 - mrSges(5,2) * t117 + Ifges(5,5) * t157 + Ifges(5,6) * t158 + Ifges(5,3) * qJDD(4) + pkin(4) * t107 - (-t192 * t143 + t190 * t217) * qJD(1) - t204;
t221 = mrSges(2,1) + mrSges(3,1);
t220 = -mrSges(5,2) - mrSges(6,2);
t218 = m(6) * t113 + t158 * mrSges(6,3);
t122 = -t127 * t188 + t189 * t130;
t209 = qJDD(1) * pkin(3) - t122;
t119 = -pkin(6) * t195 + t209;
t163 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t216;
t165 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t215;
t115 = -t161 * t216 - t158 * pkin(4) + qJDD(5) + (-qJ(5) * t183 - pkin(6)) * t195 + t209;
t162 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t216;
t210 = -m(6) * t115 + t158 * mrSges(6,1) + t162 * t216;
t106 = -m(5) * t119 + t163 * t216 + t158 * mrSges(5,1) + t220 * t157 + (-t164 - t165) * t215 + t210;
t102 = m(4) * t122 - qJDD(1) * mrSges(4,1) - t195 * mrSges(4,2) + t106;
t156 = (mrSges(5,1) * t192 - mrSges(5,2) * t190) * qJD(1);
t104 = t156 * t216 + m(5) * t116 + qJDD(4) * mrSges(5,1) + qJD(4) * t165 + (-mrSges(5,3) - mrSges(6,3)) * t157 + t211;
t105 = m(5) * t117 + t158 * mrSges(5,3) + t220 * qJDD(4) + (-t162 - t163) * qJD(4) + (-t155 - t156) * t215 + t218;
t101 = -t104 * t190 + t105 * t192;
t97 = m(4) * t123 - mrSges(4,1) * t195 + qJDD(1) * mrSges(4,2) + t101;
t93 = -t188 * t102 + t189 * t97;
t92 = t102 * t189 + t188 * t97;
t100 = t104 * t192 + t105 * t190;
t133 = -pkin(1) * t195 + t206;
t208 = m(3) * t133 + qJDD(1) * mrSges(3,3) + t93;
t207 = mrSges(6,2) * t115 - mrSges(6,3) * t112 + Ifges(6,1) * t157 + Ifges(6,4) * t158 + Ifges(6,5) * qJDD(4);
t99 = -m(4) * t184 - t100;
t138 = Ifges(6,3) * qJD(4) + (-Ifges(6,5) * t190 - Ifges(6,6) * t192) * qJD(1);
t203 = -mrSges(6,1) * t115 + mrSges(6,3) * t113 + Ifges(6,4) * t157 + Ifges(6,2) * t158 + Ifges(6,6) * qJDD(4) + qJD(4) * t142 + t138 * t216;
t136 = -qJDD(1) * pkin(1) + t205;
t202 = -m(3) * t136 + qJDD(1) * mrSges(3,1) + t195 * mrSges(3,3) - t92;
t139 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t190 - Ifges(5,6) * t192) * qJD(1);
t94 = Ifges(5,4) * t157 + Ifges(5,2) * t158 + Ifges(5,6) * qJDD(4) + t139 * t216 + qJD(4) * t143 - mrSges(5,1) * t119 + mrSges(5,3) * t117 - pkin(4) * (t157 * mrSges(6,2) + t164 * t215 - t210) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t162 - t155 * t215 + t218) + t203;
t95 = mrSges(5,2) * t119 - mrSges(5,3) * t116 + Ifges(5,1) * t157 + Ifges(5,4) * t158 + Ifges(5,5) * qJDD(4) - qJ(5) * t107 - t217 * qJD(4) + (-t138 - t139) * t215 + t207;
t86 = mrSges(4,2) * t184 - mrSges(4,3) * t122 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t195 - pkin(6) * t100 - t190 * t94 + t192 * t95;
t87 = -mrSges(4,1) * t184 + mrSges(4,3) * t123 + t195 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t100 - t225;
t201 = mrSges(3,2) * t136 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t195 * Ifges(3,6) - qJ(3) * t92 - t188 * t87 + t189 * t86;
t200 = mrSges(3,2) * t133 - pkin(2) * t99 - qJ(3) * t93 - t188 * t86 - t189 * t87;
t199 = mrSges(4,1) * t122 - mrSges(4,2) * t123 - Ifges(4,3) * qJDD(1) + pkin(3) * t106 + pkin(6) * t101 + t190 * t95 + t192 * t94;
t198 = -mrSges(3,1) * t136 + mrSges(3,3) * t133 + Ifges(3,2) * qJDD(1) - pkin(2) * t92 - t199;
t196 = -mrSges(2,2) * t167 + mrSges(2,1) * t166 + Ifges(2,3) * qJDD(1) + t198 + qJ(2) * (-mrSges(3,1) * t195 + t208) + pkin(1) * t202;
t98 = -m(3) * g(3) + t99;
t89 = m(2) * t166 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t195 + t202;
t88 = m(2) * t167 - qJDD(1) * mrSges(2,2) - t195 * t221 + t208;
t84 = -mrSges(2,2) * g(3) - mrSges(2,3) * t166 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t195 - qJ(2) * t98 + t201;
t83 = mrSges(2,3) * t167 - pkin(1) * t98 + (Ifges(3,4) + Ifges(2,5)) * t195 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t221 * g(3) + t200;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t193 * t84 - t191 * t83 - pkin(5) * (t191 * t88 + t193 * t89), t84, t201, t86, t95, -qJD(4) * t140 - t138 * t215 + t207; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t191 * t84 + t193 * t83 + pkin(5) * (-t191 * t89 + t193 * t88), t83, t198, t87, t94, t203; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t196, t196, -mrSges(3,1) * g(3) - t195 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t200, t199, t225, -t140 * t216 - t204;];
m_new = t1;
