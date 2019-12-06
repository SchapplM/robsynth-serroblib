% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:40
% EndTime: 2019-12-05 16:15:43
% DurationCPUTime: 3.52s
% Computational Cost: add. (55090->193), mult. (78756->248), div. (0->0), fcn. (51705->10), ass. (0->95)
t187 = qJD(2) + qJD(3);
t183 = t187 ^ 2;
t192 = sin(pkin(8));
t194 = cos(pkin(8));
t172 = t192 * g(1) - t194 * g(2);
t173 = -t194 * g(1) - t192 * g(2);
t197 = sin(qJ(2));
t200 = cos(qJ(2));
t158 = t200 * t172 - t197 * t173;
t155 = qJDD(2) * pkin(2) + t158;
t159 = t197 * t172 + t200 * t173;
t201 = qJD(2) ^ 2;
t156 = -t201 * pkin(2) + t159;
t196 = sin(qJ(3));
t199 = cos(qJ(3));
t141 = t196 * t155 + t199 * t156;
t184 = qJDD(2) + qJDD(3);
t137 = -t183 * pkin(3) + t184 * qJ(4) + t141;
t191 = sin(pkin(9));
t190 = -g(3) + qJDD(1);
t193 = cos(pkin(9));
t222 = qJD(4) * t187;
t223 = t193 * t190 - 0.2e1 * t191 * t222;
t228 = pkin(4) * t193;
t130 = (-pkin(7) * t184 + t183 * t228 - t137) * t191 + t223;
t134 = t191 * t190 + (t137 + 0.2e1 * t222) * t193;
t225 = t184 * t193;
t186 = t193 ^ 2;
t226 = t183 * t186;
t131 = -pkin(4) * t226 + pkin(7) * t225 + t134;
t195 = sin(qJ(5));
t198 = cos(qJ(5));
t128 = t198 * t130 - t195 * t131;
t209 = -t191 * t195 + t193 * t198;
t162 = t209 * t187;
t210 = t191 * t198 + t193 * t195;
t163 = t210 * t187;
t149 = -t162 * mrSges(6,1) + t163 * mrSges(6,2);
t151 = t162 * qJD(5) + t210 * t184;
t160 = -qJD(5) * mrSges(6,2) + t162 * mrSges(6,3);
t125 = m(6) * t128 + qJDD(5) * mrSges(6,1) - t151 * mrSges(6,3) + qJD(5) * t160 - t163 * t149;
t129 = t195 * t130 + t198 * t131;
t150 = -t163 * qJD(5) + t209 * t184;
t161 = qJD(5) * mrSges(6,1) - t163 * mrSges(6,3);
t126 = m(6) * t129 - qJDD(5) * mrSges(6,2) + t150 * mrSges(6,3) - qJD(5) * t161 + t162 * t149;
t116 = t198 * t125 + t195 * t126;
t133 = -t191 * t137 + t223;
t143 = Ifges(6,4) * t163 + Ifges(6,2) * t162 + Ifges(6,6) * qJD(5);
t144 = Ifges(6,1) * t163 + Ifges(6,4) * t162 + Ifges(6,5) * qJD(5);
t206 = -mrSges(6,1) * t128 + mrSges(6,2) * t129 - Ifges(6,5) * t151 - Ifges(6,6) * t150 - Ifges(6,3) * qJDD(5) - t163 * t143 + t162 * t144;
t215 = Ifges(5,4) * t191 + Ifges(5,2) * t193;
t216 = Ifges(5,1) * t191 + Ifges(5,4) * t193;
t229 = -mrSges(5,1) * t133 + mrSges(5,2) * t134 - pkin(4) * t116 - (t191 * t215 - t193 * t216) * t183 + t206;
t227 = mrSges(5,2) * t191;
t212 = mrSges(5,3) * t184 + (-mrSges(5,1) * t193 + t227) * t183;
t114 = m(5) * t133 - t212 * t191 + t116;
t217 = -t195 * t125 + t198 * t126;
t115 = m(5) * t134 + t212 * t193 + t217;
t218 = -t191 * t114 + t193 * t115;
t108 = m(4) * t141 - t183 * mrSges(4,1) - t184 * mrSges(4,2) + t218;
t140 = t199 * t155 - t196 * t156;
t213 = qJDD(4) - t140;
t136 = -t184 * pkin(3) - t183 * qJ(4) + t213;
t185 = t191 ^ 2;
t132 = (-pkin(3) - t228) * t184 + (-qJ(4) + (-t185 - t186) * pkin(7)) * t183 + t213;
t207 = m(6) * t132 - t150 * mrSges(6,1) + t151 * mrSges(6,2) - t162 * t160 + t163 * t161;
t204 = -m(5) * t136 + mrSges(5,1) * t225 - t207 + (t183 * t185 + t226) * mrSges(5,3);
t120 = m(4) * t140 - t183 * mrSges(4,2) + (mrSges(4,1) - t227) * t184 + t204;
t103 = t196 * t108 + t199 * t120;
t98 = m(3) * t158 + qJDD(2) * mrSges(3,1) - t201 * mrSges(3,2) + t103;
t219 = t199 * t108 - t196 * t120;
t99 = m(3) * t159 - t201 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t219;
t93 = t197 * t99 + t200 * t98;
t214 = Ifges(5,5) * t191 + Ifges(5,6) * t193;
t224 = t183 * t214;
t110 = t193 * t114 + t191 * t115;
t221 = m(4) * t190 + t110;
t220 = -t197 * t98 + t200 * t99;
t142 = Ifges(6,5) * t163 + Ifges(6,6) * t162 + Ifges(6,3) * qJD(5);
t117 = -mrSges(6,1) * t132 + mrSges(6,3) * t129 + Ifges(6,4) * t151 + Ifges(6,2) * t150 + Ifges(6,6) * qJDD(5) + qJD(5) * t144 - t163 * t142;
t118 = mrSges(6,2) * t132 - mrSges(6,3) * t128 + Ifges(6,1) * t151 + Ifges(6,4) * t150 + Ifges(6,5) * qJDD(5) - qJD(5) * t143 + t162 * t142;
t101 = -mrSges(5,1) * t136 + mrSges(5,3) * t134 - pkin(4) * t207 + pkin(7) * t217 + t198 * t117 + t195 * t118 + t215 * t184 - t191 * t224;
t105 = mrSges(5,2) * t136 - mrSges(5,3) * t133 - pkin(7) * t116 - t195 * t117 + t198 * t118 + t216 * t184 + t193 * t224;
t208 = -mrSges(4,2) * t141 + qJ(4) * t218 + t193 * t101 + t191 * t105 + pkin(3) * (-t184 * t227 + t204) + mrSges(4,1) * t140 + Ifges(4,3) * t184;
t205 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + Ifges(3,3) * qJDD(2) + pkin(2) * t103 + t208;
t202 = mrSges(2,1) * t172 - mrSges(2,2) * t173 + pkin(1) * t93 + t205;
t94 = (Ifges(4,6) - t214) * t184 - mrSges(4,1) * t190 + t183 * Ifges(4,5) + mrSges(4,3) * t141 - pkin(3) * t110 + t229;
t91 = m(2) * t173 + t220;
t90 = m(2) * t172 + t93;
t89 = mrSges(4,2) * t190 - mrSges(4,3) * t140 + Ifges(4,5) * t184 - t183 * Ifges(4,6) - qJ(4) * t110 - t191 * t101 + t193 * t105;
t88 = mrSges(3,2) * t190 - mrSges(3,3) * t158 + Ifges(3,5) * qJDD(2) - t201 * Ifges(3,6) - pkin(6) * t103 - t196 * t94 + t199 * t89;
t87 = -mrSges(3,1) * t190 + mrSges(3,3) * t159 + t201 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t221 + pkin(6) * t219 + t196 * t89 + t199 * t94;
t86 = mrSges(2,2) * t190 - mrSges(2,3) * t172 - pkin(5) * t93 - t197 * t87 + t200 * t88;
t85 = -mrSges(2,1) * t190 + mrSges(2,3) * t173 + t197 * t88 + t200 * t87 - pkin(1) * (m(3) * t190 + t221) + pkin(5) * t220;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t194 * t86 - t192 * t85 - qJ(1) * (t192 * t91 + t194 * t90), t86, t88, t89, t105, t118; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t192 * t86 + t194 * t85 + qJ(1) * (-t192 * t90 + t194 * t91), t85, t87, t94, t101, t117; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t202, t202, t205, t208, t214 * t184 - t229, -t206;];
m_new = t1;
