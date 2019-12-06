% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:20
% EndTime: 2019-12-05 17:51:25
% DurationCPUTime: 3.42s
% Computational Cost: add. (52562->203), mult. (73380->265), div. (0->0), fcn. (40823->10), ass. (0->101)
t232 = 2 * qJD(4);
t196 = sin(qJ(1));
t199 = cos(qJ(1));
t174 = t199 * g(2) + t196 * g(3);
t167 = qJDD(1) * pkin(1) + t174;
t173 = t196 * g(2) - g(3) * t199;
t200 = qJD(1) ^ 2;
t168 = -pkin(1) * t200 + t173;
t191 = sin(pkin(8));
t193 = cos(pkin(8));
t151 = t193 * t167 - t168 * t191;
t145 = qJDD(1) * pkin(2) + t151;
t152 = t191 * t167 + t193 * t168;
t146 = -pkin(2) * t200 + t152;
t195 = sin(qJ(3));
t198 = cos(qJ(3));
t141 = t145 * t195 + t146 * t198;
t186 = (qJD(1) + qJD(3));
t184 = t186 ^ 2;
t185 = qJDD(1) + qJDD(3);
t138 = -pkin(3) * t184 + qJ(4) * t185 + t141;
t231 = (t186 * t232) + t138;
t190 = sin(pkin(9));
t226 = t186 * t190;
t192 = cos(pkin(9));
t224 = t192 * t186;
t189 = -g(1) + qJDD(2);
t134 = t190 * t189 + t192 * t231;
t214 = -pkin(4) * t192 - pkin(7) * t190;
t163 = t214 * t186;
t132 = t163 * t224 + t134;
t140 = t198 * t145 - t146 * t195;
t207 = -t184 * qJ(4) + qJDD(4) - t140;
t135 = (-pkin(3) + t214) * t185 + t207;
t194 = sin(qJ(5));
t197 = cos(qJ(5));
t128 = -t132 * t194 + t135 * t197;
t171 = qJD(5) - t224;
t220 = t194 * t226;
t154 = -mrSges(6,2) * t171 - mrSges(6,3) * t220;
t156 = (mrSges(6,1) * t194 + mrSges(6,2) * t197) * t226;
t221 = qJD(5) * t186;
t158 = (t185 * t197 - t194 * t221) * t190;
t225 = t192 * t185;
t170 = qJDD(5) - t225;
t219 = t197 * t226;
t126 = m(6) * t128 + mrSges(6,1) * t170 - mrSges(6,3) * t158 + t154 * t171 - t156 * t219;
t129 = t132 * t197 + t135 * t194;
t155 = mrSges(6,1) * t171 - mrSges(6,3) * t219;
t157 = (-t185 * t194 - t197 * t221) * t190;
t127 = m(6) * t129 - mrSges(6,2) * t170 + mrSges(6,3) * t157 - t155 * t171 - t156 * t220;
t120 = -t194 * t126 + t127 * t197;
t223 = t192 * t189;
t131 = -t223 + (t138 + (t232 + t163) * t186) * t190;
t148 = Ifges(6,3) * t171 + (Ifges(6,5) * t197 - Ifges(6,6) * t194) * t226;
t150 = Ifges(6,5) * t171 + (Ifges(6,1) * t197 - Ifges(6,4) * t194) * t226;
t121 = -mrSges(6,1) * t131 + mrSges(6,3) * t129 + Ifges(6,4) * t158 + Ifges(6,2) * t157 + Ifges(6,6) * t170 - t148 * t219 + t150 * t171;
t149 = Ifges(6,6) * t171 + (Ifges(6,4) * t197 - Ifges(6,2) * t194) * t226;
t122 = mrSges(6,2) * t131 - mrSges(6,3) * t128 + Ifges(6,1) * t158 + Ifges(6,4) * t157 + Ifges(6,5) * t170 - t148 * t220 - t149 * t171;
t133 = -t190 * t231 + t223;
t208 = -m(6) * t131 + t157 * mrSges(6,1) - t158 * mrSges(6,2);
t210 = -t154 * t194 - t155 * t197;
t213 = Ifges(5,1) * t190 + Ifges(5,4) * t192;
t230 = -((Ifges(5,4) * t190 + Ifges(5,2) * t192) * t226 - t213 * t224) * t186 - mrSges(5,1) * t133 + mrSges(5,2) * t134 - pkin(4) * (t210 * t226 + t208) - pkin(7) * t120 - t197 * t121 - t194 * t122;
t229 = mrSges(5,2) * t190;
t159 = (-mrSges(5,1) * t192 + t229) * t186;
t227 = mrSges(5,3) * t185;
t117 = m(5) * t134 + (t159 * t186 + t227) * t192 + t120;
t124 = m(5) * t133 + (-t227 + (-t159 + t210) * t186) * t190 + t208;
t215 = t117 * t192 - t124 * t190;
t110 = m(4) * t141 - mrSges(4,1) * t184 - mrSges(4,2) * t185 + t215;
t119 = t197 * t126 + t194 * t127;
t137 = -pkin(3) * t185 + t207;
t205 = -m(5) * t137 + mrSges(5,1) * t225 - t119 + (t190 ^ 2 + t192 ^ 2) * mrSges(5,3) * t184;
t114 = m(4) * t140 - t184 * mrSges(4,2) + (mrSges(4,1) - t229) * t185 + t205;
t103 = t110 * t195 + t114 * t198;
t100 = m(3) * t151 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t200 + t103;
t216 = t110 * t198 - t114 * t195;
t101 = m(3) * t152 - mrSges(3,1) * t200 - qJDD(1) * mrSges(3,2) + t216;
t95 = t100 * t193 + t101 * t191;
t112 = t117 * t190 + t124 * t192;
t218 = m(4) * t189 + t112;
t217 = -t100 * t191 + t101 * t193;
t212 = Ifges(5,5) * t190 + Ifges(5,6) * t192;
t211 = t149 * t197 + t150 * t194;
t160 = t212 * t186;
t105 = mrSges(5,2) * t137 - mrSges(5,3) * t133 - pkin(7) * t119 - t194 * t121 + t197 * t122 + t160 * t224 + t185 * t213;
t204 = mrSges(6,1) * t128 - mrSges(6,2) * t129 + Ifges(6,5) * t158 + Ifges(6,6) * t157 + Ifges(6,3) * t170;
t107 = Ifges(5,2) * t225 - mrSges(5,1) * t137 + mrSges(5,3) * t134 - pkin(4) * t119 + (Ifges(5,4) * t185 + (-t160 - t211) * t186) * t190 - t204;
t206 = -mrSges(4,2) * t141 + qJ(4) * t215 + t190 * t105 + t192 * t107 + pkin(3) * (-t185 * t229 + t205) + mrSges(4,1) * t140 + Ifges(4,3) * t185;
t203 = mrSges(3,1) * t151 - mrSges(3,2) * t152 + Ifges(3,3) * qJDD(1) + pkin(2) * t103 + t206;
t201 = mrSges(2,1) * t174 - mrSges(2,2) * t173 + Ifges(2,3) * qJDD(1) + pkin(1) * t95 + t203;
t96 = -mrSges(4,1) * t189 + mrSges(4,3) * t141 + t184 * Ifges(4,5) - pkin(3) * t112 + (Ifges(4,6) - t212) * t185 + t230;
t93 = m(2) * t174 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t200 + t95;
t92 = m(2) * t173 - mrSges(2,1) * t200 - qJDD(1) * mrSges(2,2) + t217;
t91 = mrSges(4,2) * t189 - mrSges(4,3) * t140 + Ifges(4,5) * t185 - Ifges(4,6) * t184 - qJ(4) * t112 + t105 * t192 - t107 * t190;
t90 = mrSges(3,2) * t189 - mrSges(3,3) * t151 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t200 - pkin(6) * t103 - t195 * t96 + t198 * t91;
t89 = -mrSges(3,1) * t189 + mrSges(3,3) * t152 + t200 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t218 + pkin(6) * t216 + t195 * t91 + t198 * t96;
t88 = -mrSges(2,2) * g(1) - mrSges(2,3) * t174 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t200 - qJ(2) * t95 - t191 * t89 + t193 * t90;
t87 = Ifges(2,6) * qJDD(1) + t200 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t173 + t191 * t90 + t193 * t89 - pkin(1) * (m(3) * t189 + t218) + qJ(2) * t217;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t201, t88, t90, t91, t105, t122; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t196 * t88 - t199 * t87 - pkin(5) * (-t196 * t93 + t199 * t92), t87, t89, t96, t107, t121; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t199 * t88 - t196 * t87 + pkin(5) * (-t196 * t92 - t199 * t93), t201, t203, t206, t185 * t212 - t230, t211 * t226 + t204;];
m_new = t1;
