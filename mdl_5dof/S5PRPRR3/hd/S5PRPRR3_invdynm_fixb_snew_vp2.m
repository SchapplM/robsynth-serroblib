% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:44
% EndTime: 2019-12-05 15:46:50
% DurationCPUTime: 3.19s
% Computational Cost: add. (42389->214), mult. (76028->276), div. (0->0), fcn. (46790->10), ass. (0->94)
t193 = sin(pkin(8));
t195 = cos(pkin(8));
t179 = -t195 * g(1) - t193 * g(2);
t191 = -g(3) + qJDD(1);
t198 = sin(qJ(2));
t201 = cos(qJ(2));
t161 = -t198 * t179 + t201 * t191;
t157 = qJDD(2) * pkin(2) + t161;
t162 = t201 * t179 + t198 * t191;
t202 = qJD(2) ^ 2;
t158 = -t202 * pkin(2) + t162;
t192 = sin(pkin(9));
t194 = cos(pkin(9));
t143 = t192 * t157 + t194 * t158;
t139 = -t202 * pkin(3) + qJDD(2) * pkin(6) + t143;
t178 = t193 * g(1) - t195 * g(2);
t177 = qJDD(3) - t178;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t135 = -t197 * t139 + t200 * t177;
t218 = qJD(2) * qJD(4);
t217 = t200 * t218;
t174 = t197 * qJDD(2) + t217;
t132 = (-t174 + t217) * pkin(7) + (t197 * t200 * t202 + qJDD(4)) * pkin(4) + t135;
t136 = t200 * t139 + t197 * t177;
t175 = t200 * qJDD(2) - t197 * t218;
t220 = qJD(2) * t197;
t182 = qJD(4) * pkin(4) - pkin(7) * t220;
t190 = t200 ^ 2;
t133 = -t190 * t202 * pkin(4) + t175 * pkin(7) - qJD(4) * t182 + t136;
t196 = sin(qJ(5));
t199 = cos(qJ(5));
t130 = t199 * t132 - t196 * t133;
t166 = (-t196 * t197 + t199 * t200) * qJD(2);
t148 = t166 * qJD(5) + t199 * t174 + t196 * t175;
t167 = (t196 * t200 + t197 * t199) * qJD(2);
t153 = -t166 * mrSges(6,1) + t167 * mrSges(6,2);
t188 = qJD(4) + qJD(5);
t159 = -t188 * mrSges(6,2) + t166 * mrSges(6,3);
t187 = qJDD(4) + qJDD(5);
t127 = m(6) * t130 + t187 * mrSges(6,1) - t148 * mrSges(6,3) - t167 * t153 + t188 * t159;
t131 = t196 * t132 + t199 * t133;
t147 = -t167 * qJD(5) - t196 * t174 + t199 * t175;
t160 = t188 * mrSges(6,1) - t167 * mrSges(6,3);
t128 = m(6) * t131 - t187 * mrSges(6,2) + t147 * mrSges(6,3) + t166 * t153 - t188 * t160;
t117 = t199 * t127 + t196 * t128;
t164 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t197 + Ifges(5,2) * t200) * qJD(2);
t165 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t197 + Ifges(5,4) * t200) * qJD(2);
t150 = Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * t188;
t151 = Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * t188;
t206 = -mrSges(6,1) * t130 + mrSges(6,2) * t131 - Ifges(6,5) * t148 - Ifges(6,6) * t147 - Ifges(6,3) * t187 - t167 * t150 + t166 * t151;
t221 = mrSges(5,1) * t135 - mrSges(5,2) * t136 + Ifges(5,5) * t174 + Ifges(5,6) * t175 + Ifges(5,3) * qJDD(4) + pkin(4) * t117 + (t197 * t164 - t200 * t165) * qJD(2) - t206;
t173 = (-mrSges(5,1) * t200 + mrSges(5,2) * t197) * qJD(2);
t219 = qJD(2) * t200;
t181 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t219;
t115 = m(5) * t135 + qJDD(4) * mrSges(5,1) - t174 * mrSges(5,3) + qJD(4) * t181 - t173 * t220 + t117;
t180 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t220;
t213 = -t196 * t127 + t199 * t128;
t116 = m(5) * t136 - qJDD(4) * mrSges(5,2) + t175 * mrSges(5,3) - qJD(4) * t180 + t173 * t219 + t213;
t214 = -t197 * t115 + t200 * t116;
t106 = m(4) * t143 - t202 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t214;
t142 = t194 * t157 - t192 * t158;
t210 = -qJDD(2) * pkin(3) - t142;
t138 = -t202 * pkin(6) + t210;
t134 = t182 * t220 - t175 * pkin(4) + (-pkin(7) * t190 - pkin(6)) * t202 + t210;
t207 = m(6) * t134 - t147 * mrSges(6,1) + t148 * mrSges(6,2) - t166 * t159 + t167 * t160;
t204 = -m(5) * t138 + t175 * mrSges(5,1) - t174 * mrSges(5,2) - t180 * t220 + t181 * t219 - t207;
t121 = m(4) * t142 + qJDD(2) * mrSges(4,1) - t202 * mrSges(4,2) + t204;
t101 = t192 * t106 + t194 * t121;
t110 = t200 * t115 + t197 * t116;
t215 = t194 * t106 - t192 * t121;
t100 = m(3) * t162 - t202 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t215;
t99 = m(3) * t161 + qJDD(2) * mrSges(3,1) - t202 * mrSges(3,2) + t101;
t216 = t201 * t100 - t198 * t99;
t212 = m(4) * t177 + t110;
t149 = Ifges(6,5) * t167 + Ifges(6,6) * t166 + Ifges(6,3) * t188;
t118 = -mrSges(6,1) * t134 + mrSges(6,3) * t131 + Ifges(6,4) * t148 + Ifges(6,2) * t147 + Ifges(6,6) * t187 - t167 * t149 + t188 * t151;
t119 = mrSges(6,2) * t134 - mrSges(6,3) * t130 + Ifges(6,1) * t148 + Ifges(6,4) * t147 + Ifges(6,5) * t187 + t166 * t149 - t188 * t150;
t163 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t197 + Ifges(5,6) * t200) * qJD(2);
t103 = mrSges(5,2) * t138 - mrSges(5,3) * t135 + Ifges(5,1) * t174 + Ifges(5,4) * t175 + Ifges(5,5) * qJDD(4) - pkin(7) * t117 - qJD(4) * t164 - t196 * t118 + t199 * t119 + t163 * t219;
t97 = -mrSges(5,1) * t138 + mrSges(5,3) * t136 + Ifges(5,4) * t174 + Ifges(5,2) * t175 + Ifges(5,6) * qJDD(4) - pkin(4) * t207 + pkin(7) * t213 + qJD(4) * t165 + t199 * t118 + t196 * t119 - t163 * t220;
t91 = mrSges(4,2) * t177 - mrSges(4,3) * t142 + Ifges(4,5) * qJDD(2) - t202 * Ifges(4,6) - pkin(6) * t110 + t200 * t103 - t197 * t97;
t95 = -mrSges(4,1) * t177 + mrSges(4,3) * t143 + t202 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t110 - t221;
t87 = mrSges(3,1) * t178 + mrSges(3,3) * t162 + t202 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t212 + qJ(3) * t215 + t192 * t91 + t194 * t95;
t90 = -mrSges(3,2) * t178 - mrSges(3,3) * t161 + Ifges(3,5) * qJDD(2) - t202 * Ifges(3,6) - qJ(3) * t101 - t192 * t95 + t194 * t91;
t209 = -mrSges(2,2) * t179 + pkin(5) * t216 + t198 * t90 + t201 * t87 + pkin(1) * (m(3) * t178 - t212) + mrSges(2,1) * t178;
t208 = mrSges(4,1) * t142 - mrSges(4,2) * t143 + Ifges(4,3) * qJDD(2) + pkin(3) * t204 + pkin(6) * t214 + t197 * t103 + t200 * t97;
t205 = mrSges(3,1) * t161 - mrSges(3,2) * t162 + Ifges(3,3) * qJDD(2) + pkin(2) * t101 + t208;
t107 = (m(2) + m(3)) * t178 - t212;
t94 = t198 * t100 + t201 * t99;
t92 = m(2) * t179 + t216;
t88 = -mrSges(2,1) * t191 + mrSges(2,3) * t179 - pkin(1) * t94 - t205;
t85 = mrSges(2,2) * t191 - mrSges(2,3) * t178 - pkin(5) * t94 - t198 * t87 + t201 * t90;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t195 * t85 - t193 * t88 - qJ(1) * (t195 * t107 + t193 * t92), t85, t90, t91, t103, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t193 * t85 + t195 * t88 + qJ(1) * (-t193 * t107 + t195 * t92), t88, t87, t95, t97, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209, t209, t205, t208, t221, -t206;];
m_new = t1;
