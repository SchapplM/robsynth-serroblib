% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:02
% EndTime: 2019-12-05 17:38:05
% DurationCPUTime: 1.51s
% Computational Cost: add. (15234->225), mult. (29077->265), div. (0->0), fcn. (14228->6), ass. (0->87)
t255 = 2 * qJD(1);
t219 = sin(qJ(1));
t222 = cos(qJ(1));
t190 = -t222 * g(1) - t219 * g(2);
t254 = qJDD(1) * qJ(2) + (qJD(2) * t255) + t190;
t189 = t219 * g(1) - t222 * g(2);
t224 = qJD(1) ^ 2;
t170 = -qJDD(1) * pkin(1) - t224 * qJ(2) + qJDD(2) - t189;
t236 = qJDD(1) * qJ(3) + (qJD(3) * t255) - t170;
t253 = mrSges(3,2) - mrSges(4,3);
t252 = (Ifges(3,4) - Ifges(4,5));
t251 = Ifges(2,6) - Ifges(3,5);
t250 = mrSges(4,3) * t224;
t163 = qJDD(3) + (-pkin(1) - qJ(3)) * t224 + t254;
t158 = -qJDD(1) * pkin(6) + t163;
t218 = sin(qJ(4));
t221 = cos(qJ(4));
t151 = t218 * g(3) + t221 * t158;
t247 = qJD(1) * qJD(4);
t244 = t218 * t247;
t184 = qJDD(1) * t221 - t244;
t132 = (-t184 - t244) * pkin(7) + (-t218 * t221 * t224 + qJDD(4)) * pkin(4) + t151;
t152 = -g(3) * t221 + t218 * t158;
t183 = -qJDD(1) * t218 - t221 * t247;
t248 = qJD(1) * t221;
t188 = qJD(4) * pkin(4) - pkin(7) * t248;
t211 = t218 ^ 2;
t133 = -pkin(4) * t211 * t224 + pkin(7) * t183 - qJD(4) * t188 + t152;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t130 = t132 * t220 - t133 * t217;
t174 = (-t217 * t221 - t218 * t220) * qJD(1);
t143 = qJD(5) * t174 + t183 * t217 + t184 * t220;
t175 = (-t217 * t218 + t220 * t221) * qJD(1);
t153 = -mrSges(6,1) * t174 + mrSges(6,2) * t175;
t198 = qJD(4) + qJD(5);
t164 = -mrSges(6,2) * t198 + mrSges(6,3) * t174;
t197 = qJDD(4) + qJDD(5);
t126 = m(6) * t130 + mrSges(6,1) * t197 - t143 * mrSges(6,3) - t153 * t175 + t164 * t198;
t131 = t132 * t217 + t133 * t220;
t142 = -qJD(5) * t175 + t183 * t220 - t184 * t217;
t165 = mrSges(6,1) * t198 - mrSges(6,3) * t175;
t127 = m(6) * t131 - mrSges(6,2) * t197 + t142 * mrSges(6,3) + t153 * t174 - t165 * t198;
t115 = t220 * t126 + t217 * t127;
t182 = (mrSges(5,1) * t218 + mrSges(5,2) * t221) * qJD(1);
t249 = qJD(1) * t218;
t186 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t249;
t112 = m(5) * t151 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t184 + qJD(4) * t186 - t182 * t248 + t115;
t187 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t248;
t242 = -t217 * t126 + t220 * t127;
t113 = m(5) * t152 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t183 - qJD(4) * t187 - t182 * t249 + t242;
t107 = t221 * t112 + t218 * t113;
t243 = -t218 * t112 + t221 * t113;
t241 = m(4) * t163 + qJDD(1) * mrSges(4,2) + t107;
t135 = t188 * t248 - t183 * pkin(4) + (-pkin(7) * t211 - pkin(6)) * t224 + t236;
t240 = m(6) * t135 - t142 * mrSges(6,1) + t143 * mrSges(6,2) - t174 * t164 + t175 * t165;
t145 = Ifges(6,4) * t175 + Ifges(6,2) * t174 + Ifges(6,6) * t198;
t146 = Ifges(6,1) * t175 + Ifges(6,4) * t174 + Ifges(6,5) * t198;
t237 = mrSges(6,1) * t130 - mrSges(6,2) * t131 + Ifges(6,5) * t143 + Ifges(6,6) * t142 + Ifges(6,3) * t197 + t175 * t145 - t174 * t146;
t144 = Ifges(6,5) * t175 + Ifges(6,6) * t174 + Ifges(6,3) * t198;
t116 = -mrSges(6,1) * t135 + mrSges(6,3) * t131 + Ifges(6,4) * t143 + Ifges(6,2) * t142 + Ifges(6,6) * t197 - t144 * t175 + t146 * t198;
t117 = mrSges(6,2) * t135 - mrSges(6,3) * t130 + Ifges(6,1) * t143 + Ifges(6,4) * t142 + Ifges(6,5) * t197 + t144 * t174 - t145 * t198;
t157 = -pkin(6) * t224 + t236;
t171 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t221 - Ifges(5,6) * t218) * qJD(1);
t172 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t221 - Ifges(5,2) * t218) * qJD(1);
t100 = mrSges(5,2) * t157 - mrSges(5,3) * t151 + Ifges(5,1) * t184 + Ifges(5,4) * t183 + Ifges(5,5) * qJDD(4) - pkin(7) * t115 - qJD(4) * t172 - t116 * t217 + t117 * t220 - t171 * t249;
t233 = -m(5) * t157 + t183 * mrSges(5,1) - t184 * mrSges(5,2) - t186 * t249 - t187 * t248 - t240;
t173 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t221 - Ifges(5,4) * t218) * qJD(1);
t97 = -mrSges(5,1) * t157 + mrSges(5,3) * t152 + Ifges(5,4) * t184 + Ifges(5,2) * t183 + Ifges(5,6) * qJDD(4) - pkin(4) * t240 + pkin(7) * t242 + qJD(4) * t173 + t220 * t116 + t217 * t117 - t171 * t248;
t235 = -mrSges(4,1) * t236 + mrSges(4,2) * g(3) + (t224 * Ifges(4,4)) + Ifges(4,5) * qJDD(1) + pkin(3) * t233 + pkin(6) * t243 + t218 * t100 + t221 * t97;
t168 = pkin(1) * t224 - t254;
t234 = -m(3) * t168 + (t224 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t241;
t232 = mrSges(4,2) * t163 + mrSges(4,3) * t236 + Ifges(4,1) * qJDD(1) - pkin(6) * t107 + t221 * t100 - t218 * t97;
t120 = -m(4) * t236 - t224 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t233;
t231 = mrSges(3,1) * t170 + pkin(2) * t120 + t235;
t230 = mrSges(5,1) * t151 - mrSges(5,2) * t152 + Ifges(5,5) * t184 + Ifges(5,6) * t183 + Ifges(5,3) * qJDD(4) + pkin(4) * t115 + t172 * t248 + t173 * t249 + t237;
t229 = mrSges(3,2) * t170 - mrSges(3,3) * t168 + Ifges(3,1) * qJDD(1) - qJ(3) * t120 + t232;
t228 = -m(3) * t170 + t224 * mrSges(3,3) - t120;
t227 = mrSges(4,1) * t163 - Ifges(4,4) * qJDD(1) + pkin(3) * t107 + t230;
t226 = -mrSges(2,2) * t190 + qJ(2) * (t234 - t250) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t228) + mrSges(2,1) * t189 + Ifges(2,3) * qJDD(1) + t229;
t225 = mrSges(3,1) * t168 + pkin(2) * (-t241 + t250) + qJ(3) * (-m(4) * g(3) + t243) - t227;
t118 = t228 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) + m(2) * t189 - mrSges(2,2) * t224;
t104 = (-m(3) - m(4)) * g(3) + t243;
t101 = m(2) * t190 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t224 + t234;
t95 = -t225 + ((Ifges(2,5) - t252) * t224) + t251 * qJDD(1) + (mrSges(2,1) - t253) * g(3) + mrSges(2,3) * t190 - pkin(1) * t104;
t94 = t231 + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) - t251 * t224 + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t189 - qJ(2) * t104;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t222 * t94 - t219 * t95 - pkin(5) * (t101 * t219 + t118 * t222), t94, t229, t232, t100, t117; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t219 * t94 + t222 * t95 + pkin(5) * (t101 * t222 - t118 * t219), t95, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t224 * Ifges(3,5)) - t231, -mrSges(4,3) * g(3) - t224 * Ifges(4,5) - t227, t97, t116; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t226, t226, Ifges(3,5) * qJDD(1) + t253 * g(3) + (t252 * t224) + t225, t235, t230, t237;];
m_new = t1;
