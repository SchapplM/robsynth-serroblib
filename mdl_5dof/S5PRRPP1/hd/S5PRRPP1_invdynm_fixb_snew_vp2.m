% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:17
% EndTime: 2019-12-05 16:06:25
% DurationCPUTime: 3.02s
% Computational Cost: add. (32126->253), mult. (69273->313), div. (0->0), fcn. (42162->8), ass. (0->97)
t223 = sin(pkin(7));
t224 = cos(pkin(7));
t206 = t223 * g(1) - t224 * g(2);
t207 = -t224 * g(1) - t223 * g(2);
t226 = sin(qJ(2));
t228 = cos(qJ(2));
t182 = t226 * t206 + t228 * t207;
t230 = qJD(2) ^ 2;
t177 = -t230 * pkin(2) + qJDD(2) * pkin(6) + t182;
t221 = -g(3) + qJDD(1);
t225 = sin(qJ(3));
t227 = cos(qJ(3));
t151 = -t225 * t177 + t227 * t221;
t249 = qJD(2) * qJD(3);
t247 = t227 * t249;
t203 = t225 * qJDD(2) + t247;
t146 = (-t203 + t247) * qJ(4) + (t225 * t227 * t230 + qJDD(3)) * pkin(3) + t151;
t152 = t227 * t177 + t225 * t221;
t204 = t227 * qJDD(2) - t225 * t249;
t251 = qJD(2) * t225;
t208 = qJD(3) * pkin(3) - qJ(4) * t251;
t220 = t227 ^ 2;
t147 = -t220 * t230 * pkin(3) + t204 * qJ(4) - qJD(3) * t208 + t152;
t222 = sin(pkin(8));
t250 = qJD(2) * t227;
t254 = cos(pkin(8));
t190 = t222 * t251 - t254 * t250;
t256 = -2 * qJD(4);
t143 = t222 * t146 + t254 * t147 + t190 * t256;
t178 = t222 * t203 - t254 * t204;
t191 = (t222 * t227 + t254 * t225) * qJD(2);
t186 = qJD(3) * mrSges(5,1) - t191 * mrSges(5,3);
t163 = t190 * pkin(4) - t191 * qJ(5);
t229 = qJD(3) ^ 2;
t136 = -t229 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t190 * t163 + t143;
t187 = -qJD(3) * mrSges(6,1) + t191 * mrSges(6,2);
t248 = m(6) * t136 + qJDD(3) * mrSges(6,3) + qJD(3) * t187;
t164 = t190 * mrSges(6,1) - t191 * mrSges(6,3);
t252 = -t190 * mrSges(5,1) - t191 * mrSges(5,2) - t164;
t255 = -mrSges(5,3) - mrSges(6,2);
t126 = m(5) * t143 - qJDD(3) * mrSges(5,2) - qJD(3) * t186 + t255 * t178 + t252 * t190 + t248;
t239 = t254 * t146 - t222 * t147;
t142 = t191 * t256 + t239;
t179 = t254 * t203 + t222 * t204;
t185 = -qJD(3) * mrSges(5,2) - t190 * mrSges(5,3);
t138 = -qJDD(3) * pkin(4) - t229 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t163) * t191 - t239;
t188 = -t190 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t243 = -m(6) * t138 + qJDD(3) * mrSges(6,1) + qJD(3) * t188;
t127 = m(5) * t142 + qJDD(3) * mrSges(5,1) + qJD(3) * t185 + t255 * t179 + t252 * t191 + t243;
t120 = t222 * t126 + t254 * t127;
t193 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t225 + Ifges(4,2) * t227) * qJD(2);
t194 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t225 + Ifges(4,4) * t227) * qJD(2);
t157 = Ifges(5,4) * t191 - Ifges(5,2) * t190 + Ifges(5,6) * qJD(3);
t159 = Ifges(5,1) * t191 - Ifges(5,4) * t190 + Ifges(5,5) * qJD(3);
t154 = Ifges(6,5) * t191 + Ifges(6,6) * qJD(3) + Ifges(6,3) * t190;
t158 = Ifges(6,1) * t191 + Ifges(6,4) * qJD(3) + Ifges(6,5) * t190;
t236 = mrSges(6,1) * t138 - mrSges(6,3) * t136 - Ifges(6,4) * t179 - Ifges(6,2) * qJDD(3) - Ifges(6,6) * t178 + t191 * t154 - t190 * t158;
t233 = mrSges(5,2) * t143 - t190 * t159 - qJ(5) * (-t178 * mrSges(6,2) - t190 * t164 + t248) - pkin(4) * (-t179 * mrSges(6,2) - t191 * t164 + t243) - mrSges(5,1) * t142 - t191 * t157 + Ifges(5,6) * t178 - Ifges(5,5) * t179 - Ifges(5,3) * qJDD(3) + t236;
t257 = mrSges(4,1) * t151 - mrSges(4,2) * t152 + Ifges(4,5) * t203 + Ifges(4,6) * t204 + Ifges(4,3) * qJDD(3) + pkin(3) * t120 + (t225 * t193 - t227 * t194) * qJD(2) - t233;
t202 = (-mrSges(4,1) * t227 + mrSges(4,2) * t225) * qJD(2);
t210 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t250;
t116 = m(4) * t151 + qJDD(3) * mrSges(4,1) - t203 * mrSges(4,3) + qJD(3) * t210 - t202 * t251 + t120;
t209 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t251;
t244 = t254 * t126 - t222 * t127;
t117 = m(4) * t152 - qJDD(3) * mrSges(4,2) + t204 * mrSges(4,3) - qJD(3) * t209 + t202 * t250 + t244;
t245 = -t225 * t116 + t227 * t117;
t110 = m(3) * t182 - t230 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t245;
t181 = t228 * t206 - t226 * t207;
t240 = -qJDD(2) * pkin(2) - t181;
t176 = -t230 * pkin(6) + t240;
t148 = -t204 * pkin(3) + qJDD(4) + t208 * t251 + (-qJ(4) * t220 - pkin(6)) * t230 + t240;
t140 = -0.2e1 * qJD(5) * t191 + (qJD(3) * t190 - t179) * qJ(5) + (qJD(3) * t191 + t178) * pkin(4) + t148;
t133 = m(6) * t140 + t178 * mrSges(6,1) - t179 * mrSges(6,3) - t191 * t187 + t190 * t188;
t234 = m(5) * t148 + t178 * mrSges(5,1) + t179 * mrSges(5,2) + t190 * t185 + t191 * t186 + t133;
t232 = -m(4) * t176 + t204 * mrSges(4,1) - t203 * mrSges(4,2) - t209 * t251 + t210 * t250 - t234;
t122 = m(3) * t181 + qJDD(2) * mrSges(3,1) - t230 * mrSges(3,2) + t232;
t107 = t226 * t110 + t228 * t122;
t112 = t227 * t116 + t225 * t117;
t156 = Ifges(6,4) * t191 + Ifges(6,2) * qJD(3) + Ifges(6,6) * t190;
t253 = -Ifges(5,5) * t191 + Ifges(5,6) * t190 - Ifges(5,3) * qJD(3) - t156;
t246 = t228 * t110 - t226 * t122;
t242 = -mrSges(6,1) * t140 + mrSges(6,2) * t136;
t238 = mrSges(6,2) * t138 - mrSges(6,3) * t140 + Ifges(6,1) * t179 + Ifges(6,4) * qJDD(3) + Ifges(6,5) * t178 + qJD(3) * t154;
t118 = -mrSges(5,1) * t148 + mrSges(5,3) * t143 - pkin(4) * t133 + t253 * t191 + (Ifges(5,4) - Ifges(6,5)) * t179 + (-Ifges(5,2) - Ifges(6,3)) * t178 + (Ifges(5,6) - Ifges(6,6)) * qJDD(3) + (t158 + t159) * qJD(3) + t242;
t119 = mrSges(5,2) * t148 - mrSges(5,3) * t142 + Ifges(5,1) * t179 - Ifges(5,4) * t178 + Ifges(5,5) * qJDD(3) - qJ(5) * t133 - qJD(3) * t157 + t253 * t190 + t238;
t192 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t225 + Ifges(4,6) * t227) * qJD(2);
t101 = -mrSges(4,1) * t176 + mrSges(4,3) * t152 + Ifges(4,4) * t203 + Ifges(4,2) * t204 + Ifges(4,6) * qJDD(3) - pkin(3) * t234 + qJ(4) * t244 + qJD(3) * t194 + t254 * t118 + t222 * t119 - t192 * t251;
t103 = mrSges(4,2) * t176 - mrSges(4,3) * t151 + Ifges(4,1) * t203 + Ifges(4,4) * t204 + Ifges(4,5) * qJDD(3) - qJ(4) * t120 - qJD(3) * t193 - t222 * t118 + t254 * t119 + t192 * t250;
t237 = mrSges(3,1) * t181 - mrSges(3,2) * t182 + Ifges(3,3) * qJDD(2) + pkin(2) * t232 + pkin(6) * t245 + t227 * t101 + t225 * t103;
t235 = mrSges(2,1) * t206 - mrSges(2,2) * t207 + pkin(1) * t107 + t237;
t105 = m(2) * t207 + t246;
t104 = m(2) * t206 + t107;
t99 = -mrSges(3,1) * t221 + mrSges(3,3) * t182 + t230 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t112 - t257;
t98 = mrSges(3,2) * t221 - mrSges(3,3) * t181 + Ifges(3,5) * qJDD(2) - t230 * Ifges(3,6) - pkin(6) * t112 - t225 * t101 + t227 * t103;
t97 = mrSges(2,2) * t221 - mrSges(2,3) * t206 - pkin(5) * t107 - t226 * t99 + t228 * t98;
t96 = -mrSges(2,1) * t221 + mrSges(2,3) * t207 + t226 * t98 + t228 * t99 - pkin(1) * (m(3) * t221 + t112) + pkin(5) * t246;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t224 * t97 - t223 * t96 - qJ(1) * (t224 * t104 + t223 * t105), t97, t98, t103, t119, -t190 * t156 + t238; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t223 * t97 + t224 * t96 + qJ(1) * (-t223 * t104 + t224 * t105), t96, t99, t101, t118, -t236; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t235, t235, t237, t257, -t233, Ifges(6,5) * t179 + Ifges(6,6) * qJDD(3) + Ifges(6,3) * t178 - qJD(3) * t158 + t191 * t156 - t242;];
m_new = t1;
