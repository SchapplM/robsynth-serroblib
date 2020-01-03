% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:24
% EndTime: 2019-12-31 17:46:26
% DurationCPUTime: 2.06s
% Computational Cost: add. (26041->203), mult. (49355->250), div. (0->0), fcn. (23645->8), ass. (0->93)
t191 = qJD(1) ^ 2;
t187 = sin(qJ(1));
t189 = cos(qJ(1));
t158 = -t189 * g(1) - t187 * g(2);
t203 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t158;
t222 = -pkin(1) - pkin(2);
t141 = t222 * t191 + t203;
t157 = t187 * g(1) - t189 * g(2);
t202 = -t191 * qJ(2) + qJDD(2) - t157;
t146 = t222 * qJDD(1) + t202;
t183 = sin(pkin(7));
t185 = cos(pkin(7));
t127 = t185 * t141 + t183 * t146;
t123 = -pkin(3) * t191 - qJDD(1) * qJ(4) + t127;
t182 = sin(pkin(8));
t184 = cos(pkin(8));
t178 = g(3) + qJDD(3);
t215 = qJD(1) * qJD(4);
t217 = t184 * t178 + 0.2e1 * t182 * t215;
t116 = (pkin(4) * t184 * t191 + pkin(6) * qJDD(1) - t123) * t182 + t217;
t120 = t182 * t178 + (t123 - (2 * t215)) * t184;
t214 = qJDD(1) * t184;
t169 = t184 ^ 2;
t218 = t169 * t191;
t117 = -pkin(4) * t218 - pkin(6) * t214 + t120;
t186 = sin(qJ(5));
t188 = cos(qJ(5));
t114 = t116 * t188 - t117 * t186;
t207 = t182 * t186 - t184 * t188;
t150 = t207 * qJD(1);
t208 = -t182 * t188 - t184 * t186;
t151 = t208 * qJD(1);
t132 = -mrSges(6,1) * t150 + mrSges(6,2) * t151;
t137 = t150 * qJD(5) + t208 * qJDD(1);
t142 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t150;
t111 = m(6) * t114 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t137 + qJD(5) * t142 - t132 * t151;
t115 = t116 * t186 + t117 * t188;
t136 = -t151 * qJD(5) + t207 * qJDD(1);
t143 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t151;
t112 = m(6) * t115 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t136 - qJD(5) * t143 + t132 * t150;
t103 = t188 * t111 + t186 * t112;
t119 = -t123 * t182 + t217;
t129 = Ifges(6,4) * t151 + Ifges(6,2) * t150 + Ifges(6,6) * qJD(5);
t130 = Ifges(6,1) * t151 + Ifges(6,4) * t150 + Ifges(6,5) * qJD(5);
t198 = -mrSges(6,1) * t114 + mrSges(6,2) * t115 - Ifges(6,5) * t137 - Ifges(6,6) * t136 - Ifges(6,3) * qJDD(5) - t151 * t129 + t150 * t130;
t211 = -Ifges(5,4) * t182 - Ifges(5,2) * t184;
t212 = -Ifges(5,1) * t182 - Ifges(5,4) * t184;
t223 = -mrSges(5,1) * t119 + mrSges(5,2) * t120 - pkin(4) * t103 + (t182 * t211 - t184 * t212) * t191 + t198;
t221 = mrSges(2,1) + mrSges(3,1);
t220 = mrSges(5,1) * t184;
t219 = mrSges(5,2) * t182;
t210 = -Ifges(5,5) * t182 - Ifges(5,6) * t184;
t216 = t191 * t210;
t126 = -t183 * t141 + t146 * t185;
t204 = qJDD(1) * pkin(3) + qJDD(4) - t126;
t122 = -qJ(4) * t191 + t204;
t168 = t182 ^ 2;
t118 = pkin(4) * t214 + (-qJ(4) + (-t168 - t169) * pkin(6)) * t191 + t204;
t200 = m(6) * t118 - t136 * mrSges(6,1) + t137 * mrSges(6,2) - t150 * t142 + t151 * t143;
t195 = -m(5) * t122 + qJDD(1) * t219 - t200 + (t168 * t191 + t218) * mrSges(5,3);
t106 = m(4) * t126 - mrSges(4,2) * t191 + (-mrSges(4,1) - t220) * qJDD(1) + t195;
t206 = mrSges(5,3) * qJDD(1) + t191 * (-t219 + t220);
t101 = m(5) * t119 + t206 * t182 + t103;
t213 = -t186 * t111 + t188 * t112;
t102 = m(5) * t120 - t206 * t184 + t213;
t99 = -t101 * t182 + t184 * t102;
t95 = m(4) * t127 - mrSges(4,1) * t191 + qJDD(1) * mrSges(4,2) + t99;
t92 = -t106 * t183 + t185 * t95;
t91 = t106 * t185 + t183 * t95;
t98 = t101 * t184 + t102 * t182;
t147 = -pkin(1) * t191 + t203;
t205 = m(3) * t147 + qJDD(1) * mrSges(3,3) + t92;
t97 = -m(4) * t178 - t98;
t149 = -qJDD(1) * pkin(1) + t202;
t201 = -m(3) * t149 + qJDD(1) * mrSges(3,1) + t191 * mrSges(3,3) - t91;
t128 = Ifges(6,5) * t151 + Ifges(6,6) * t150 + Ifges(6,3) * qJD(5);
t104 = -mrSges(6,1) * t118 + mrSges(6,3) * t115 + Ifges(6,4) * t137 + Ifges(6,2) * t136 + Ifges(6,6) * qJDD(5) + qJD(5) * t130 - t128 * t151;
t105 = mrSges(6,2) * t118 - mrSges(6,3) * t114 + Ifges(6,1) * t137 + Ifges(6,4) * t136 + Ifges(6,5) * qJDD(5) - qJD(5) * t129 + t128 * t150;
t90 = -mrSges(5,1) * t122 + mrSges(5,3) * t120 - pkin(4) * t200 + pkin(6) * t213 + t211 * qJDD(1) + t188 * t104 + t186 * t105 + t182 * t216;
t93 = mrSges(5,2) * t122 - mrSges(5,3) * t119 - pkin(6) * t103 + t212 * qJDD(1) - t104 * t186 + t105 * t188 - t184 * t216;
t84 = mrSges(4,2) * t178 - mrSges(4,3) * t126 - Ifges(4,5) * qJDD(1) - Ifges(4,6) * t191 - qJ(4) * t98 - t182 * t90 + t184 * t93;
t85 = (-Ifges(4,6) - t210) * qJDD(1) + t191 * Ifges(4,5) - mrSges(4,1) * t178 + mrSges(4,3) * t127 - pkin(3) * t98 + t223;
t199 = mrSges(3,2) * t149 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t191 * Ifges(3,6) - qJ(3) * t91 - t183 * t85 + t185 * t84;
t197 = mrSges(3,2) * t147 - pkin(2) * t97 - qJ(3) * t92 - t183 * t84 - t185 * t85;
t196 = mrSges(4,1) * t126 + pkin(3) * (-mrSges(5,1) * t214 + t195) + qJ(4) * t99 + t182 * t93 + t184 * t90 - mrSges(4,2) * t127 - Ifges(4,3) * qJDD(1);
t193 = -mrSges(3,1) * t149 + mrSges(3,3) * t147 + Ifges(3,2) * qJDD(1) - pkin(2) * t91 - t196;
t192 = -mrSges(2,2) * t158 + mrSges(2,1) * t157 + Ifges(2,3) * qJDD(1) + t193 + qJ(2) * (-mrSges(3,1) * t191 + t205) + pkin(1) * t201;
t96 = -m(3) * g(3) + t97;
t87 = m(2) * t157 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t191 + t201;
t86 = m(2) * t158 - qJDD(1) * mrSges(2,2) - t221 * t191 + t205;
t82 = -mrSges(2,2) * g(3) - mrSges(2,3) * t157 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t191 - qJ(2) * t96 + t199;
t81 = mrSges(2,3) * t158 - pkin(1) * t96 + (Ifges(3,4) + Ifges(2,5)) * t191 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t221 * g(3) + t197;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t189 * t82 - t187 * t81 - pkin(5) * (t187 * t86 + t189 * t87), t82, t199, t84, t93, t105; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t187 * t82 + t189 * t81 + pkin(5) * (-t187 * t87 + t189 * t86), t81, t193, t85, t90, t104; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t192, t192, -mrSges(3,1) * g(3) - Ifges(3,4) * t191 + Ifges(3,6) * qJDD(1) - t197, t196, t210 * qJDD(1) - t223, -t198;];
m_new = t1;
