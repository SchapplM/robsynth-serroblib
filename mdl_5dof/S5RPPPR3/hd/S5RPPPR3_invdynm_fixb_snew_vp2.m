% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:54
% EndTime: 2019-12-31 17:43:56
% DurationCPUTime: 2.15s
% Computational Cost: add. (20774->224), mult. (43673->284), div. (0->0), fcn. (24553->8), ass. (0->101)
t222 = sin(qJ(1));
t224 = cos(qJ(1));
t195 = t222 * g(1) - t224 * g(2);
t192 = qJDD(1) * pkin(1) + t195;
t196 = -t224 * g(1) - t222 * g(2);
t225 = qJD(1) ^ 2;
t193 = -t225 * pkin(1) + t196;
t218 = sin(pkin(7));
t220 = cos(pkin(7));
t170 = t218 * t192 + t220 * t193;
t265 = -t225 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t170;
t217 = sin(pkin(8));
t211 = t217 ^ 2;
t219 = cos(pkin(8));
t212 = t219 ^ 2;
t253 = t212 * t225;
t264 = t211 * t225 + t253;
t216 = -g(3) + qJDD(2);
t148 = t219 * t216 - t265 * t217;
t255 = qJ(4) * t217;
t238 = -pkin(3) * t219 - t255;
t183 = t238 * qJD(1);
t249 = qJD(1) * t217;
t144 = t183 * t249 + qJDD(4) - t148;
t138 = (-pkin(4) * t219 * t225 - pkin(6) * qJDD(1)) * t217 + t144;
t149 = t217 * t216 + t265 * t219;
t248 = qJD(1) * t219;
t146 = t183 * t248 + t149;
t245 = qJDD(1) * t219;
t139 = -pkin(4) * t253 - pkin(6) * t245 + t146;
t221 = sin(qJ(5));
t223 = cos(qJ(5));
t136 = t223 * t138 - t221 * t139;
t236 = -t217 * t221 - t219 * t223;
t177 = t236 * qJD(1);
t237 = t217 * t223 - t219 * t221;
t178 = t237 * qJD(1);
t161 = -t177 * mrSges(6,1) + t178 * mrSges(6,2);
t167 = t177 * qJD(5) + t237 * qJDD(1);
t171 = -qJD(5) * mrSges(6,2) + t177 * mrSges(6,3);
t131 = m(6) * t136 + qJDD(5) * mrSges(6,1) - t167 * mrSges(6,3) + qJD(5) * t171 - t178 * t161;
t137 = t221 * t138 + t223 * t139;
t166 = -t178 * qJD(5) + t236 * qJDD(1);
t172 = qJD(5) * mrSges(6,1) - t178 * mrSges(6,3);
t132 = m(6) * t137 - qJDD(5) * mrSges(6,2) + t166 * mrSges(6,3) - qJD(5) * t172 + t177 * t161;
t122 = -t221 * t131 + t223 * t132;
t184 = (-mrSges(5,1) * t219 - mrSges(5,3) * t217) * qJD(1);
t120 = m(5) * t146 + mrSges(5,2) * t245 + t184 * t248 + t122;
t239 = Ifges(5,5) * t217 - Ifges(5,3) * t219;
t186 = t239 * qJD(1);
t190 = (Ifges(5,1) * t217 - Ifges(5,5) * t219) * qJD(1);
t121 = t223 * t131 + t221 * t132;
t156 = Ifges(6,4) * t178 + Ifges(6,2) * t177 + Ifges(6,6) * qJD(5);
t157 = Ifges(6,1) * t178 + Ifges(6,4) * t177 + Ifges(6,5) * qJD(5);
t235 = mrSges(6,1) * t136 - mrSges(6,2) * t137 + Ifges(6,5) * t167 + Ifges(6,6) * t166 + Ifges(6,3) * qJDD(5) + t178 * t156 - t177 * t157;
t246 = qJDD(1) * t217;
t230 = -mrSges(5,1) * t144 + mrSges(5,3) * t146 + Ifges(5,4) * t246 - pkin(4) * t121 - t235;
t233 = -m(5) * t144 - t121;
t257 = Ifges(4,1) * t217;
t263 = ((t186 - (Ifges(4,4) * t217 + Ifges(4,2) * t219) * qJD(1)) * t217 + (t190 + (Ifges(4,4) * t219 + t257) * qJD(1)) * t219) * qJD(1) - mrSges(4,1) * t148 + mrSges(4,2) * t149 - pkin(3) * ((-qJDD(1) * mrSges(5,2) - qJD(1) * t184) * t217 + t233) - qJ(4) * t120 - t230;
t256 = Ifges(4,5) * t217;
t262 = t256 + (Ifges(4,6) - Ifges(5,6)) * t219;
t260 = Ifges(4,4) - Ifges(5,5);
t258 = mrSges(4,2) * t217;
t252 = t225 * qJ(3);
t185 = (-mrSges(4,1) * t219 + t258) * qJD(1);
t117 = m(4) * t148 + ((-mrSges(5,2) - mrSges(4,3)) * qJDD(1) + (-t184 - t185) * qJD(1)) * t217 + t233;
t118 = m(4) * t149 + (qJDD(1) * mrSges(4,3) + qJD(1) * t185) * t219 + t120;
t241 = -t217 * t117 + t219 * t118;
t111 = m(3) * t170 - t225 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t241;
t169 = t220 * t192 - t218 * t193;
t243 = -qJDD(3) + t169;
t234 = -0.2e1 * qJD(4) * t249 - t243;
t142 = (qJ(3) + (-t211 - t212) * pkin(6)) * t225 + (t255 + pkin(2) + (pkin(3) + pkin(4)) * t219) * qJDD(1) - t234;
t133 = -m(6) * t142 + t166 * mrSges(6,1) - t167 * mrSges(6,2) + t177 * t171 - t178 * t172;
t147 = -t252 + (-pkin(2) + t238) * qJDD(1) + t234;
t129 = m(5) * t147 - mrSges(5,1) * t245 - t264 * mrSges(5,2) - mrSges(5,3) * t246 + t133;
t154 = -qJDD(1) * pkin(2) - t243 - t252;
t227 = -m(4) * t154 + mrSges(4,1) * t245 + t264 * mrSges(4,3) - t129;
t127 = (mrSges(3,1) - t258) * qJDD(1) + t227 - t225 * mrSges(3,2) + m(3) * t169;
t108 = t218 * t111 + t220 * t127;
t113 = t219 * t117 + t217 * t118;
t242 = t220 * t111 - t218 * t127;
t187 = (Ifges(4,6) * t219 + t256) * qJD(1);
t188 = (Ifges(5,4) * t217 - Ifges(5,6) * t219) * qJD(1);
t155 = Ifges(6,5) * t178 + Ifges(6,6) * t177 + Ifges(6,3) * qJD(5);
t124 = -mrSges(6,1) * t142 + mrSges(6,3) * t137 + Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * qJDD(5) + qJD(5) * t157 - t178 * t155;
t125 = mrSges(6,2) * t142 - mrSges(6,3) * t136 + Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * qJDD(5) - qJD(5) * t156 + t177 * t155;
t228 = mrSges(5,1) * t147 - mrSges(5,2) * t146 + pkin(4) * t133 + pkin(6) * t122 + t223 * t124 + t221 * t125;
t102 = -mrSges(4,1) * t154 + mrSges(4,3) * t149 - pkin(3) * t129 + (-t187 - t188) * t249 + ((Ifges(4,2) + Ifges(5,3)) * t219 + t260 * t217) * qJDD(1) - t228;
t231 = mrSges(5,2) * t144 - mrSges(5,3) * t147 + Ifges(5,1) * t246 - pkin(6) * t121 - t221 * t124 + t223 * t125 + t188 * t248;
t104 = t187 * t248 + mrSges(4,2) * t154 - mrSges(4,3) * t148 - qJ(4) * t129 + (t260 * t219 + t257) * qJDD(1) + t231;
t232 = -mrSges(3,2) * t170 + qJ(3) * t241 + t219 * t102 + t217 * t104 + pkin(2) * (-mrSges(4,2) * t246 + t227) + mrSges(3,1) * t169 + Ifges(3,3) * qJDD(1);
t229 = mrSges(2,1) * t195 - mrSges(2,2) * t196 + Ifges(2,3) * qJDD(1) + pkin(1) * t108 + t232;
t106 = m(2) * t196 - t225 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t242;
t105 = m(2) * t195 + qJDD(1) * mrSges(2,1) - t225 * mrSges(2,2) + t108;
t100 = (Ifges(3,6) - t262) * qJDD(1) + t225 * Ifges(3,5) - mrSges(3,1) * t216 + mrSges(3,3) * t170 - pkin(2) * t113 + t263;
t99 = mrSges(3,2) * t216 - mrSges(3,3) * t169 + Ifges(3,5) * qJDD(1) - t225 * Ifges(3,6) - qJ(3) * t113 - t217 * t102 + t219 * t104;
t98 = -mrSges(2,2) * g(3) - mrSges(2,3) * t195 + Ifges(2,5) * qJDD(1) - t225 * Ifges(2,6) - qJ(2) * t108 - t218 * t100 + t220 * t99;
t97 = Ifges(2,6) * qJDD(1) + t225 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t196 + t218 * t99 + t220 * t100 - pkin(1) * (m(3) * t216 + t113) + qJ(2) * t242;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t224 * t98 - t222 * t97 - pkin(5) * (t224 * t105 + t222 * t106), t98, t99, t104, -Ifges(5,5) * t245 + t231, t125; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t222 * t98 + t224 * t97 + pkin(5) * (-t222 * t105 + t224 * t106), t97, t100, t102, t230 + (-t217 * t186 - t219 * t190) * qJD(1) - Ifges(5,6) * t245, t124; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t229, t229, t232, t262 * qJDD(1) - t263, t239 * qJDD(1) + t188 * t249 + t228, t235;];
m_new = t1;
