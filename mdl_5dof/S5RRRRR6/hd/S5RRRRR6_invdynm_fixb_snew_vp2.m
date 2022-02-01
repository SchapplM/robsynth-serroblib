% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:51
% EndTime: 2022-01-20 12:08:04
% DurationCPUTime: 6.75s
% Computational Cost: add. (160198->270), mult. (205381->343), div. (0->0), fcn. (133069->10), ass. (0->113)
t228 = qJDD(1) + qJDD(2);
t235 = sin(qJ(3));
t240 = cos(qJ(3));
t230 = qJD(1) + qJD(2);
t258 = qJD(3) * t230;
t207 = t235 * t228 + t240 * t258;
t237 = sin(qJ(1));
t242 = cos(qJ(1));
t218 = t237 * g(1) - t242 * g(2);
t212 = qJDD(1) * pkin(1) + t218;
t219 = -t242 * g(1) - t237 * g(2);
t243 = qJD(1) ^ 2;
t213 = -t243 * pkin(1) + t219;
t236 = sin(qJ(2));
t241 = cos(qJ(2));
t193 = t236 * t212 + t241 * t213;
t226 = t230 ^ 2;
t190 = -t226 * pkin(2) + t228 * pkin(7) + t193;
t259 = t235 * t190;
t262 = pkin(3) * t226;
t168 = qJDD(3) * pkin(3) - t207 * pkin(8) - t259 + (pkin(8) * t258 + t235 * t262 - g(3)) * t240;
t180 = -t235 * g(3) + t240 * t190;
t208 = t240 * t228 - t235 * t258;
t261 = t230 * t235;
t216 = qJD(3) * pkin(3) - pkin(8) * t261;
t232 = t240 ^ 2;
t169 = t208 * pkin(8) - qJD(3) * t216 - t232 * t262 + t180;
t234 = sin(qJ(4));
t239 = cos(qJ(4));
t150 = t239 * t168 - t234 * t169;
t201 = (-t234 * t235 + t239 * t240) * t230;
t176 = t201 * qJD(4) + t239 * t207 + t234 * t208;
t202 = (t234 * t240 + t235 * t239) * t230;
t227 = qJDD(3) + qJDD(4);
t229 = qJD(3) + qJD(4);
t145 = (t201 * t229 - t176) * pkin(9) + (t201 * t202 + t227) * pkin(4) + t150;
t151 = t234 * t168 + t239 * t169;
t175 = -t202 * qJD(4) - t234 * t207 + t239 * t208;
t196 = t229 * pkin(4) - t202 * pkin(9);
t197 = t201 ^ 2;
t146 = -t197 * pkin(4) + t175 * pkin(9) - t229 * t196 + t151;
t233 = sin(qJ(5));
t238 = cos(qJ(5));
t143 = t238 * t145 - t233 * t146;
t185 = t238 * t201 - t233 * t202;
t157 = t185 * qJD(5) + t233 * t175 + t238 * t176;
t186 = t233 * t201 + t238 * t202;
t164 = -t185 * mrSges(6,1) + t186 * mrSges(6,2);
t224 = qJD(5) + t229;
t177 = -t224 * mrSges(6,2) + t185 * mrSges(6,3);
t223 = qJDD(5) + t227;
t140 = m(6) * t143 + t223 * mrSges(6,1) - t157 * mrSges(6,3) - t186 * t164 + t224 * t177;
t144 = t233 * t145 + t238 * t146;
t156 = -t186 * qJD(5) + t238 * t175 - t233 * t176;
t178 = t224 * mrSges(6,1) - t186 * mrSges(6,3);
t141 = m(6) * t144 - t223 * mrSges(6,2) + t156 * mrSges(6,3) + t185 * t164 - t224 * t178;
t131 = t238 * t140 + t233 * t141;
t188 = -t201 * mrSges(5,1) + t202 * mrSges(5,2);
t194 = -t229 * mrSges(5,2) + t201 * mrSges(5,3);
t128 = m(5) * t150 + t227 * mrSges(5,1) - t176 * mrSges(5,3) - t202 * t188 + t229 * t194 + t131;
t195 = t229 * mrSges(5,1) - t202 * mrSges(5,3);
t254 = -t233 * t140 + t238 * t141;
t129 = m(5) * t151 - t227 * mrSges(5,2) + t175 * mrSges(5,3) + t201 * t188 - t229 * t195 + t254;
t124 = t239 * t128 + t234 * t129;
t179 = -t240 * g(3) - t259;
t199 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t235 + Ifges(4,2) * t240) * t230;
t200 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t235 + Ifges(4,4) * t240) * t230;
t182 = Ifges(5,4) * t202 + Ifges(5,2) * t201 + Ifges(5,6) * t229;
t183 = Ifges(5,1) * t202 + Ifges(5,4) * t201 + Ifges(5,5) * t229;
t160 = Ifges(6,4) * t186 + Ifges(6,2) * t185 + Ifges(6,6) * t224;
t161 = Ifges(6,1) * t186 + Ifges(6,4) * t185 + Ifges(6,5) * t224;
t249 = -mrSges(6,1) * t143 + mrSges(6,2) * t144 - Ifges(6,5) * t157 - Ifges(6,6) * t156 - Ifges(6,3) * t223 - t186 * t160 + t185 * t161;
t246 = -mrSges(5,1) * t150 + mrSges(5,2) * t151 - Ifges(5,5) * t176 - Ifges(5,6) * t175 - Ifges(5,3) * t227 - pkin(4) * t131 - t202 * t182 + t201 * t183 + t249;
t263 = mrSges(4,1) * t179 - mrSges(4,2) * t180 + Ifges(4,5) * t207 + Ifges(4,6) * t208 + Ifges(4,3) * qJDD(3) + pkin(3) * t124 + (t235 * t199 - t240 * t200) * t230 - t246;
t260 = t230 * t240;
t206 = (-mrSges(4,1) * t240 + mrSges(4,2) * t235) * t230;
t215 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t260;
t122 = m(4) * t179 + qJDD(3) * mrSges(4,1) - t207 * mrSges(4,3) + qJD(3) * t215 - t206 * t261 + t124;
t214 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t261;
t255 = -t234 * t128 + t239 * t129;
t123 = m(4) * t180 - qJDD(3) * mrSges(4,2) + t208 * mrSges(4,3) - qJD(3) * t214 + t206 * t260 + t255;
t256 = -t235 * t122 + t240 * t123;
t114 = m(3) * t193 - t226 * mrSges(3,1) - t228 * mrSges(3,2) + t256;
t192 = t241 * t212 - t236 * t213;
t251 = -t228 * pkin(2) - t192;
t189 = -t226 * pkin(7) + t251;
t170 = -t208 * pkin(3) + t216 * t261 + (-pkin(8) * t232 - pkin(7)) * t226 + t251;
t148 = -t175 * pkin(4) - t197 * pkin(9) + t202 * t196 + t170;
t253 = m(6) * t148 - t156 * mrSges(6,1) + t157 * mrSges(6,2) - t185 * t177 + t186 * t178;
t247 = m(5) * t170 - t175 * mrSges(5,1) + t176 * mrSges(5,2) - t201 * t194 + t202 * t195 + t253;
t245 = -m(4) * t189 + t208 * mrSges(4,1) - t207 * mrSges(4,2) - t214 * t261 + t215 * t260 - t247;
t135 = m(3) * t192 + t228 * mrSges(3,1) - t226 * mrSges(3,2) + t245;
t111 = t236 * t114 + t241 * t135;
t116 = t240 * t122 + t235 * t123;
t257 = t241 * t114 - t236 * t135;
t159 = Ifges(6,5) * t186 + Ifges(6,6) * t185 + Ifges(6,3) * t224;
t132 = -mrSges(6,1) * t148 + mrSges(6,3) * t144 + Ifges(6,4) * t157 + Ifges(6,2) * t156 + Ifges(6,6) * t223 - t186 * t159 + t224 * t161;
t133 = mrSges(6,2) * t148 - mrSges(6,3) * t143 + Ifges(6,1) * t157 + Ifges(6,4) * t156 + Ifges(6,5) * t223 + t185 * t159 - t224 * t160;
t181 = Ifges(5,5) * t202 + Ifges(5,6) * t201 + Ifges(5,3) * t229;
t117 = -mrSges(5,1) * t170 + mrSges(5,3) * t151 + Ifges(5,4) * t176 + Ifges(5,2) * t175 + Ifges(5,6) * t227 - pkin(4) * t253 + pkin(9) * t254 + t238 * t132 + t233 * t133 - t202 * t181 + t229 * t183;
t118 = mrSges(5,2) * t170 - mrSges(5,3) * t150 + Ifges(5,1) * t176 + Ifges(5,4) * t175 + Ifges(5,5) * t227 - pkin(9) * t131 - t233 * t132 + t238 * t133 + t201 * t181 - t229 * t182;
t198 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t235 + Ifges(4,6) * t240) * t230;
t105 = -mrSges(4,1) * t189 + mrSges(4,3) * t180 + Ifges(4,4) * t207 + Ifges(4,2) * t208 + Ifges(4,6) * qJDD(3) - pkin(3) * t247 + pkin(8) * t255 + qJD(3) * t200 + t239 * t117 + t234 * t118 - t198 * t261;
t107 = mrSges(4,2) * t189 - mrSges(4,3) * t179 + Ifges(4,1) * t207 + Ifges(4,4) * t208 + Ifges(4,5) * qJDD(3) - pkin(8) * t124 - qJD(3) * t199 - t234 * t117 + t239 * t118 + t198 * t260;
t250 = mrSges(3,1) * t192 - mrSges(3,2) * t193 + Ifges(3,3) * t228 + pkin(2) * t245 + pkin(7) * t256 + t240 * t105 + t235 * t107;
t248 = mrSges(2,1) * t218 - mrSges(2,2) * t219 + Ifges(2,3) * qJDD(1) + pkin(1) * t111 + t250;
t109 = m(2) * t219 - t243 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t257;
t108 = m(2) * t218 + qJDD(1) * mrSges(2,1) - t243 * mrSges(2,2) + t111;
t103 = mrSges(3,1) * g(3) + mrSges(3,3) * t193 + t226 * Ifges(3,5) + Ifges(3,6) * t228 - pkin(2) * t116 - t263;
t102 = -mrSges(3,2) * g(3) - mrSges(3,3) * t192 + Ifges(3,5) * t228 - t226 * Ifges(3,6) - pkin(7) * t116 - t235 * t105 + t240 * t107;
t101 = -mrSges(2,2) * g(3) - mrSges(2,3) * t218 + Ifges(2,5) * qJDD(1) - t243 * Ifges(2,6) - pkin(6) * t111 + t241 * t102 - t236 * t103;
t100 = Ifges(2,6) * qJDD(1) + t243 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t219 + t236 * t102 + t241 * t103 - pkin(1) * (-m(3) * g(3) + t116) + pkin(6) * t257;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t242 * t101 - t237 * t100 - pkin(5) * (t242 * t108 + t237 * t109), t101, t102, t107, t118, t133; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t237 * t101 + t242 * t100 + pkin(5) * (-t237 * t108 + t242 * t109), t100, t103, t105, t117, t132; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t248, t248, t250, t263, -t246, -t249;];
m_new = t1;
