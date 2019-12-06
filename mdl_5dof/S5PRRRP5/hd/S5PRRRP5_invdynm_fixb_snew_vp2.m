% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:50
% EndTime: 2019-12-05 16:47:58
% DurationCPUTime: 3.37s
% Computational Cost: add. (36009->252), mult. (72634->312), div. (0->0), fcn. (45317->8), ass. (0->97)
t228 = sin(qJ(4));
t229 = sin(qJ(3));
t231 = cos(qJ(4));
t232 = cos(qJ(3));
t198 = (-t228 * t229 + t231 * t232) * qJD(2);
t252 = qJD(2) * qJD(3);
t249 = t232 * t252;
t208 = qJDD(2) * t229 + t249;
t209 = qJDD(2) * t232 - t229 * t252;
t168 = qJD(4) * t198 + t208 * t231 + t209 * t228;
t199 = (t228 * t232 + t229 * t231) * qJD(2);
t181 = -mrSges(6,1) * t198 + mrSges(6,2) * t199;
t227 = sin(pkin(8));
t256 = cos(pkin(8));
t212 = -t256 * g(1) - t227 * g(2);
t226 = -g(3) + qJDD(1);
t230 = sin(qJ(2));
t233 = cos(qJ(2));
t193 = t233 * t212 + t230 * t226;
t234 = qJD(2) ^ 2;
t186 = -pkin(2) * t234 + qJDD(2) * pkin(6) + t193;
t211 = g(1) * t227 - t256 * g(2);
t170 = -t186 * t229 - t232 * t211;
t148 = (-t208 + t249) * pkin(7) + (t229 * t232 * t234 + qJDD(3)) * pkin(3) + t170;
t171 = t232 * t186 - t229 * t211;
t254 = qJD(2) * t229;
t215 = qJD(3) * pkin(3) - pkin(7) * t254;
t225 = t232 ^ 2;
t149 = -pkin(3) * t225 * t234 + pkin(7) * t209 - qJD(3) * t215 + t171;
t143 = t231 * t148 - t149 * t228;
t223 = qJDD(3) + qJDD(4);
t224 = qJD(3) + qJD(4);
t135 = -0.2e1 * qJD(5) * t199 + (t198 * t224 - t168) * qJ(5) + (t198 * t199 + t223) * pkin(4) + t143;
t187 = -mrSges(6,2) * t224 + mrSges(6,3) * t198;
t251 = m(6) * t135 + t223 * mrSges(6,1) + t224 * t187;
t132 = -t168 * mrSges(6,3) - t181 * t199 + t251;
t144 = t228 * t148 + t231 * t149;
t167 = -qJD(4) * t199 - t208 * t228 + t209 * t231;
t175 = Ifges(5,4) * t199 + Ifges(5,2) * t198 + Ifges(5,6) * t224;
t176 = Ifges(6,1) * t199 + Ifges(6,4) * t198 + Ifges(6,5) * t224;
t177 = Ifges(5,1) * t199 + Ifges(5,4) * t198 + Ifges(5,5) * t224;
t189 = pkin(4) * t224 - qJ(5) * t199;
t194 = t198 ^ 2;
t138 = -pkin(4) * t194 + t167 * qJ(5) + 0.2e1 * qJD(5) * t198 - t189 * t224 + t144;
t174 = Ifges(6,4) * t199 + Ifges(6,2) * t198 + Ifges(6,6) * t224;
t240 = -mrSges(6,1) * t135 + mrSges(6,2) * t138 - Ifges(6,5) * t168 - Ifges(6,6) * t167 - Ifges(6,3) * t223 - t199 * t174;
t259 = mrSges(5,1) * t143 - mrSges(5,2) * t144 + Ifges(5,5) * t168 + Ifges(5,6) * t167 + Ifges(5,3) * t223 + pkin(4) * t132 + t199 * t175 - t240 - (t177 + t176) * t198;
t182 = -mrSges(5,1) * t198 + mrSges(5,2) * t199;
t188 = -mrSges(5,2) * t224 + mrSges(5,3) * t198;
t125 = m(5) * t143 + mrSges(5,1) * t223 + t188 * t224 + (-t181 - t182) * t199 + (-mrSges(5,3) - mrSges(6,3)) * t168 + t251;
t190 = mrSges(6,1) * t224 - mrSges(6,3) * t199;
t191 = mrSges(5,1) * t224 - mrSges(5,3) * t199;
t250 = m(6) * t138 + t167 * mrSges(6,3) + t198 * t181;
t129 = m(5) * t144 + t167 * mrSges(5,3) + t182 * t198 + (-t190 - t191) * t224 + (-mrSges(5,2) - mrSges(6,2)) * t223 + t250;
t123 = t231 * t125 + t228 * t129;
t196 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t229 + Ifges(4,2) * t232) * qJD(2);
t197 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t229 + Ifges(4,4) * t232) * qJD(2);
t258 = mrSges(4,1) * t170 - mrSges(4,2) * t171 + Ifges(4,5) * t208 + Ifges(4,6) * t209 + Ifges(4,3) * qJDD(3) + pkin(3) * t123 + (t196 * t229 - t197 * t232) * qJD(2) + t259;
t253 = qJD(2) * t232;
t207 = (-mrSges(4,1) * t232 + mrSges(4,2) * t229) * qJD(2);
t214 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t253;
t121 = m(4) * t170 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t208 + qJD(3) * t214 - t207 * t254 + t123;
t213 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t254;
t246 = -t228 * t125 + t231 * t129;
t122 = m(4) * t171 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t209 - qJD(3) * t213 + t207 * t253 + t246;
t117 = -t121 * t229 + t232 * t122;
t113 = m(3) * t193 - mrSges(3,1) * t234 - qJDD(2) * mrSges(3,2) + t117;
t192 = -t230 * t212 + t226 * t233;
t243 = -qJDD(2) * pkin(2) - t192;
t185 = -pkin(6) * t234 + t243;
t150 = -pkin(3) * t209 + t215 * t254 + (-pkin(7) * t225 - pkin(6)) * t234 + t243;
t141 = -t167 * pkin(4) - qJ(5) * t194 + t189 * t199 + qJDD(5) + t150;
t245 = m(6) * t141 - t167 * mrSges(6,1) + t168 * mrSges(6,2) - t198 * t187 + t199 * t190;
t238 = m(5) * t150 - t167 * mrSges(5,1) + t168 * mrSges(5,2) - t198 * t188 + t191 * t199 + t245;
t130 = -m(4) * t185 + t209 * mrSges(4,1) - mrSges(4,2) * t208 - t213 * t254 + t214 * t253 - t238;
t128 = m(3) * t192 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t234 + t130;
t247 = t233 * t113 - t128 * t230;
t116 = t121 * t232 + t122 * t229;
t242 = -mrSges(6,1) * t141 + mrSges(6,3) * t138 + Ifges(6,4) * t168 + Ifges(6,2) * t167 + Ifges(6,6) * t223 + t224 * t176;
t172 = Ifges(6,5) * t199 + Ifges(6,6) * t198 + Ifges(6,3) * t224;
t173 = Ifges(5,5) * t199 + Ifges(5,6) * t198 + Ifges(5,3) * t224;
t118 = Ifges(5,4) * t168 + Ifges(5,2) * t167 + Ifges(5,6) * t223 + t224 * t177 - mrSges(5,1) * t150 + mrSges(5,3) * t144 - pkin(4) * t245 + qJ(5) * (-mrSges(6,2) * t223 - t190 * t224 + t250) + (-t173 - t172) * t199 + t242;
t239 = mrSges(6,2) * t141 - mrSges(6,3) * t135 + Ifges(6,1) * t168 + Ifges(6,4) * t167 + Ifges(6,5) * t223 + t198 * t172;
t120 = mrSges(5,2) * t150 - mrSges(5,3) * t143 + Ifges(5,1) * t168 + Ifges(5,4) * t167 + Ifges(5,5) * t223 - qJ(5) * t132 + t173 * t198 + (-t174 - t175) * t224 + t239;
t195 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t229 + Ifges(4,6) * t232) * qJD(2);
t107 = -mrSges(4,1) * t185 + mrSges(4,3) * t171 + Ifges(4,4) * t208 + Ifges(4,2) * t209 + Ifges(4,6) * qJDD(3) - pkin(3) * t238 + pkin(7) * t246 + qJD(3) * t197 + t231 * t118 + t228 * t120 - t195 * t254;
t108 = mrSges(4,2) * t185 - mrSges(4,3) * t170 + Ifges(4,1) * t208 + Ifges(4,4) * t209 + Ifges(4,5) * qJDD(3) - pkin(7) * t123 - qJD(3) * t196 - t118 * t228 + t120 * t231 + t195 * t253;
t104 = -mrSges(3,2) * t211 - mrSges(3,3) * t192 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t234 - pkin(6) * t116 - t107 * t229 + t108 * t232;
t106 = mrSges(3,1) * t211 + mrSges(3,3) * t193 + t234 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t116 - t258;
t241 = -mrSges(2,2) * t212 + pkin(5) * t247 + t230 * t104 + t233 * t106 + pkin(1) * (m(3) * t211 - t116) + mrSges(2,1) * t211;
t236 = mrSges(3,1) * t192 - mrSges(3,2) * t193 + Ifges(3,3) * qJDD(2) + pkin(2) * t130 + pkin(6) * t117 + t107 * t232 + t108 * t229;
t114 = (m(2) + m(3)) * t211 - t116;
t111 = t113 * t230 + t128 * t233;
t109 = m(2) * t212 + t247;
t102 = -mrSges(2,1) * t226 + mrSges(2,3) * t212 - pkin(1) * t111 - t236;
t101 = mrSges(2,2) * t226 - mrSges(2,3) * t211 - pkin(5) * t111 + t104 * t233 - t106 * t230;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t256 * t101 - t227 * t102 - qJ(1) * (t227 * t109 + t256 * t114), t101, t104, t108, t120, -t174 * t224 + t239; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t227 * t101 + t256 * t102 + qJ(1) * (t256 * t109 - t227 * t114), t102, t106, t107, t118, -t199 * t172 + t242; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t241, t241, t236, t258, t259, -t198 * t176 - t240;];
m_new = t1;
