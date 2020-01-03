% Calculate vector of cutting torques with Newton-Euler for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:39
% EndTime: 2019-12-31 16:40:40
% DurationCPUTime: 1.16s
% Computational Cost: add. (7186->193), mult. (16609->248), div. (0->0), fcn. (9324->6), ass. (0->87)
t186 = sin(qJ(1));
t188 = cos(qJ(1));
t163 = -t188 * g(1) - t186 * g(2);
t189 = qJD(1) ^ 2;
t227 = -t189 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t163;
t183 = sin(pkin(6));
t179 = t183 ^ 2;
t184 = cos(pkin(6));
t180 = t184 ^ 2;
t215 = t180 * t189;
t226 = t179 * t189 + t215;
t138 = -t184 * g(3) - t227 * t183;
t139 = -t183 * g(3) + t227 * t184;
t153 = (-mrSges(4,1) * t184 - mrSges(4,3) * t183) * qJD(1);
t202 = Ifges(4,5) * t183 - Ifges(4,3) * t184;
t155 = t202 * qJD(1);
t159 = (Ifges(4,1) * t183 - Ifges(4,5) * t184) * qJD(1);
t217 = qJ(3) * t183;
t201 = -pkin(2) * t184 - t217;
t152 = t201 * qJD(1);
t210 = qJD(1) * t183;
t122 = t152 * t210 + qJDD(3) - t138;
t209 = qJD(1) * t184;
t124 = t152 * t209 + t139;
t114 = (-pkin(3) * t184 * t189 - pkin(5) * qJDD(1)) * t183 + t122;
t206 = qJDD(1) * t184;
t115 = -pkin(3) * t215 - pkin(5) * t206 + t124;
t185 = sin(qJ(4));
t187 = cos(qJ(4));
t112 = t187 * t114 - t185 * t115;
t113 = t185 * t114 + t187 * t115;
t199 = -t183 * t185 - t184 * t187;
t149 = t199 * qJD(1);
t200 = t183 * t187 - t184 * t185;
t150 = t200 * qJD(1);
t126 = Ifges(5,4) * t150 + Ifges(5,2) * t149 + Ifges(5,6) * qJD(4);
t127 = Ifges(5,1) * t150 + Ifges(5,4) * t149 + Ifges(5,5) * qJD(4);
t136 = -t150 * qJD(4) + t199 * qJDD(1);
t137 = t149 * qJD(4) + t200 * qJDD(1);
t198 = mrSges(5,1) * t112 - mrSges(5,2) * t113 + Ifges(5,5) * t137 + Ifges(5,6) * t136 + Ifges(5,3) * qJDD(4) + t150 * t126 - t149 * t127;
t207 = qJDD(1) * t183;
t130 = -t149 * mrSges(5,1) + t150 * mrSges(5,2);
t140 = -qJD(4) * mrSges(5,2) + t149 * mrSges(5,3);
t107 = m(5) * t112 + qJDD(4) * mrSges(5,1) - t137 * mrSges(5,3) + qJD(4) * t140 - t150 * t130;
t141 = qJD(4) * mrSges(5,1) - t150 * mrSges(5,3);
t108 = m(5) * t113 - qJDD(4) * mrSges(5,2) + t136 * mrSges(5,3) - qJD(4) * t141 + t149 * t130;
t98 = t187 * t107 + t185 * t108;
t193 = -mrSges(4,1) * t122 + mrSges(4,3) * t124 + Ifges(4,4) * t207 - pkin(3) * t98 - t198;
t196 = -m(4) * t122 - t98;
t219 = Ifges(3,1) * t183;
t99 = -t185 * t107 + t187 * t108;
t97 = m(4) * t124 + mrSges(4,2) * t206 + t153 * t209 + t99;
t225 = ((t155 - (Ifges(3,4) * t183 + Ifges(3,2) * t184) * qJD(1)) * t183 + (t159 + (Ifges(3,4) * t184 + t219) * qJD(1)) * t184) * qJD(1) - mrSges(3,1) * t138 + mrSges(3,2) * t139 - pkin(2) * ((-qJDD(1) * mrSges(4,2) - qJD(1) * t153) * t183 + t196) - qJ(3) * t97 - t193;
t218 = Ifges(3,5) * t183;
t224 = t218 + (Ifges(3,6) - Ifges(4,6)) * t184;
t222 = Ifges(3,4) - Ifges(4,5);
t220 = mrSges(3,2) * t183;
t214 = t189 * qJ(2);
t162 = t186 * g(1) - t188 * g(2);
t204 = -qJDD(2) + t162;
t154 = (-mrSges(3,1) * t184 + t220) * qJD(1);
t94 = m(3) * t138 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t153 - t154) * qJD(1)) * t183 + t196;
t95 = m(3) * t139 + (qJDD(1) * mrSges(3,3) + qJD(1) * t154) * t184 + t97;
t203 = -t183 * t94 + t184 * t95;
t197 = -0.2e1 * qJD(3) * t210 - t204;
t119 = (qJ(2) + (-t179 - t180) * pkin(5)) * t189 + (t217 + pkin(1) + (pkin(2) + pkin(3)) * t184) * qJDD(1) - t197;
t109 = -m(5) * t119 + t136 * mrSges(5,1) - t137 * mrSges(5,2) + t149 * t140 - t150 * t141;
t135 = -t214 + (-pkin(1) + t201) * qJDD(1) + t197;
t105 = m(4) * t135 - mrSges(4,1) * t206 - mrSges(4,2) * t226 - mrSges(4,3) * t207 + t109;
t148 = -qJDD(1) * pkin(1) - t204 - t214;
t191 = -m(3) * t148 + mrSges(3,1) * t206 + mrSges(3,3) * t226 - t105;
t156 = (Ifges(3,6) * t184 + t218) * qJD(1);
t157 = (Ifges(4,4) * t183 - Ifges(4,6) * t184) * qJD(1);
t125 = Ifges(5,5) * t150 + Ifges(5,6) * t149 + Ifges(5,3) * qJD(4);
t101 = -mrSges(5,1) * t119 + mrSges(5,3) * t113 + Ifges(5,4) * t137 + Ifges(5,2) * t136 + Ifges(5,6) * qJDD(4) + qJD(4) * t127 - t150 * t125;
t102 = mrSges(5,2) * t119 - mrSges(5,3) * t112 + Ifges(5,1) * t137 + Ifges(5,4) * t136 + Ifges(5,5) * qJDD(4) - qJD(4) * t126 + t149 * t125;
t192 = mrSges(4,1) * t135 - mrSges(4,2) * t124 + pkin(3) * t109 + pkin(5) * t99 + t187 * t101 + t185 * t102;
t87 = -mrSges(3,1) * t148 + mrSges(3,3) * t139 - pkin(2) * t105 + (-t156 - t157) * t210 + ((Ifges(3,2) + Ifges(4,3)) * t184 + t222 * t183) * qJDD(1) - t192;
t194 = mrSges(4,2) * t122 - mrSges(4,3) * t135 + Ifges(4,1) * t207 - pkin(5) * t98 - t185 * t101 + t187 * t102 + t157 * t209;
t89 = t156 * t209 + mrSges(3,2) * t148 - mrSges(3,3) * t138 - qJ(3) * t105 + (t222 * t184 + t219) * qJDD(1) + t194;
t195 = -mrSges(2,2) * t163 + qJ(2) * t203 + t183 * t89 + t184 * t87 + pkin(1) * (-mrSges(3,2) * t207 + t191) + mrSges(2,1) * t162 + Ifges(2,3) * qJDD(1);
t103 = (mrSges(2,1) - t220) * qJDD(1) + t191 - t189 * mrSges(2,2) + m(2) * t162;
t92 = t183 * t95 + t184 * t94;
t90 = m(2) * t163 - t189 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t203;
t85 = mrSges(2,1) * g(3) + (Ifges(2,6) - t224) * qJDD(1) + t189 * Ifges(2,5) + mrSges(2,3) * t163 - pkin(1) * t92 + t225;
t84 = -mrSges(2,2) * g(3) - mrSges(2,3) * t162 + Ifges(2,5) * qJDD(1) - t189 * Ifges(2,6) - qJ(2) * t92 - t183 * t87 + t184 * t89;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t188 * t84 - t186 * t85 - pkin(4) * (t188 * t103 + t186 * t90), t84, t89, -Ifges(4,5) * t206 + t194, t102; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t186 * t84 + t188 * t85 + pkin(4) * (-t186 * t103 + t188 * t90), t85, t87, -Ifges(4,6) * t206 + (-t183 * t155 - t184 * t159) * qJD(1) + t193, t101; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t195, t195, t224 * qJDD(1) - t225, t202 * qJDD(1) + t157 * t210 + t192, t198;];
m_new = t1;
