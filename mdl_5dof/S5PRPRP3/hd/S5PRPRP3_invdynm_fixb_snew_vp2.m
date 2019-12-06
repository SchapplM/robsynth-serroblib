% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:43
% EndTime: 2019-12-05 15:32:49
% DurationCPUTime: 1.73s
% Computational Cost: add. (18026->210), mult. (31286->260), div. (0->0), fcn. (16415->8), ass. (0->87)
t198 = sin(qJ(4));
t200 = cos(qJ(4));
t170 = (-mrSges(6,1) * t200 + mrSges(6,2) * t198) * qJD(2);
t221 = qJD(2) * qJD(4);
t217 = t200 * t221;
t172 = qJDD(2) * t198 + t217;
t195 = sin(pkin(7));
t197 = cos(pkin(7));
t178 = -g(1) * t197 - g(2) * t195;
t193 = -g(3) + qJDD(1);
t199 = sin(qJ(2));
t201 = cos(qJ(2));
t144 = -t178 * t199 + t201 * t193;
t142 = qJDD(2) * pkin(2) + t144;
t145 = t201 * t178 + t199 * t193;
t202 = qJD(2) ^ 2;
t143 = -pkin(2) * t202 + t145;
t194 = sin(pkin(8));
t196 = cos(pkin(8));
t138 = t194 * t142 + t196 * t143;
t135 = -pkin(3) * t202 + qJDD(2) * pkin(6) + t138;
t177 = g(1) * t195 - t197 * g(2);
t176 = qJDD(3) - t177;
t157 = t200 * t176;
t220 = qJD(2) * qJD(5);
t228 = pkin(4) * t202;
t127 = qJDD(4) * pkin(4) + t157 + (-t172 + t217) * qJ(5) + (t200 * t228 - t135 - 0.2e1 * t220) * t198;
t222 = qJD(2) * t200;
t182 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t222;
t219 = m(6) * t127 + qJDD(4) * mrSges(6,1) + qJD(4) * t182;
t223 = qJD(2) * t198;
t122 = -t172 * mrSges(6,3) - t170 * t223 + t219;
t131 = -t198 * t135 + t157;
t132 = t200 * t135 + t198 * t176;
t153 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t198 + Ifges(5,2) * t200) * qJD(2);
t154 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t198 + Ifges(6,4) * t200) * qJD(2);
t155 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t198 + Ifges(5,4) * t200) * qJD(2);
t173 = qJDD(2) * t200 - t198 * t221;
t179 = qJD(4) * pkin(4) - qJ(5) * t223;
t192 = t200 ^ 2;
t128 = qJ(5) * t173 - qJD(4) * t179 - t192 * t228 + 0.2e1 * t200 * t220 + t132;
t152 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t198 + Ifges(6,2) * t200) * qJD(2);
t208 = -mrSges(6,1) * t127 + mrSges(6,2) * t128 - Ifges(6,5) * t172 - Ifges(6,6) * t173 - Ifges(6,3) * qJDD(4) - t152 * t223;
t230 = mrSges(5,1) * t131 - mrSges(5,2) * t132 + Ifges(5,5) * t172 + Ifges(5,6) * t173 + Ifges(5,3) * qJDD(4) + pkin(4) * t122 - (-t198 * t153 + (t154 + t155) * t200) * qJD(2) - t208;
t227 = -mrSges(5,2) - mrSges(6,2);
t171 = (-mrSges(5,1) * t200 + mrSges(5,2) * t198) * qJD(2);
t183 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t222;
t118 = m(5) * t131 + qJDD(4) * mrSges(5,1) + qJD(4) * t183 + (-mrSges(5,3) - mrSges(6,3)) * t172 + (-t170 - t171) * t223 + t219;
t218 = m(6) * t128 + t173 * mrSges(6,3) + t170 * t222;
t180 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t223;
t224 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t223 - t180;
t119 = m(5) * t132 + t173 * mrSges(5,3) + qJD(4) * t224 + qJDD(4) * t227 + t171 * t222 + t218;
t214 = -t118 * t198 + t119 * t200;
t107 = m(4) * t138 - mrSges(4,1) * t202 - qJDD(2) * mrSges(4,2) + t214;
t137 = t196 * t142 - t194 * t143;
t211 = -qJDD(2) * pkin(3) - t137;
t134 = -pkin(6) * t202 + t211;
t130 = t179 * t223 - t173 * pkin(4) + qJDD(5) + (-qJ(5) * t192 - pkin(6)) * t202 + t211;
t212 = -m(6) * t130 + t173 * mrSges(6,1) + t182 * t222;
t205 = -m(5) * t134 + t173 * mrSges(5,1) + t172 * t227 + t183 * t222 + t223 * t224 + t212;
t114 = m(4) * t137 + qJDD(2) * mrSges(4,1) - t202 * mrSges(4,2) + t205;
t100 = t107 * t194 + t114 * t196;
t111 = t118 * t200 + t119 * t198;
t98 = m(3) * t144 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t202 + t100;
t215 = t107 * t196 - t114 * t194;
t99 = m(3) * t145 - mrSges(3,1) * t202 - qJDD(2) * mrSges(3,2) + t215;
t216 = -t199 * t98 + t201 * t99;
t213 = m(4) * t176 + t111;
t150 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t198 + Ifges(6,6) * t200) * qJD(2);
t151 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t198 + Ifges(5,6) * t200) * qJD(2);
t209 = -mrSges(6,1) * t130 + mrSges(6,3) * t128 + Ifges(6,4) * t172 + Ifges(6,2) * t173 + Ifges(6,6) * qJDD(4) + qJD(4) * t154;
t102 = Ifges(5,4) * t172 + Ifges(5,2) * t173 + Ifges(5,6) * qJDD(4) + qJD(4) * t155 - mrSges(5,1) * t134 + mrSges(5,3) * t132 - pkin(4) * (t172 * mrSges(6,2) - t212) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t180 + t218) + (-pkin(4) * t180 - t150 - t151) * t223 + t209;
t207 = mrSges(6,2) * t130 - mrSges(6,3) * t127 + Ifges(6,1) * t172 + Ifges(6,4) * t173 + Ifges(6,5) * qJDD(4) + t150 * t222;
t104 = t151 * t222 + mrSges(5,2) * t134 - mrSges(5,3) * t131 + Ifges(5,1) * t172 + Ifges(5,4) * t173 + Ifges(5,5) * qJDD(4) - qJ(5) * t122 + (-t152 - t153) * qJD(4) + t207;
t92 = mrSges(4,2) * t176 - mrSges(4,3) * t137 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t202 - pkin(6) * t111 - t102 * t198 + t104 * t200;
t96 = -mrSges(4,1) * t176 + mrSges(4,3) * t138 + t202 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t111 - t230;
t89 = mrSges(3,1) * t177 + mrSges(3,3) * t145 + t202 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t213 + qJ(3) * t215 + t194 * t92 + t196 * t96;
t91 = -mrSges(3,2) * t177 - mrSges(3,3) * t144 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t202 - qJ(3) * t100 - t194 * t96 + t196 * t92;
t210 = -mrSges(2,2) * t178 + pkin(5) * t216 + t199 * t91 + t201 * t89 + pkin(1) * (m(3) * t177 - t213) + mrSges(2,1) * t177;
t206 = mrSges(4,1) * t137 - mrSges(4,2) * t138 + Ifges(4,3) * qJDD(2) + pkin(3) * t205 + pkin(6) * t214 + t102 * t200 + t104 * t198;
t204 = mrSges(3,1) * t144 - mrSges(3,2) * t145 + Ifges(3,3) * qJDD(2) + pkin(2) * t100 + t206;
t108 = (m(2) + m(3)) * t177 - t213;
t95 = t199 * t99 + t201 * t98;
t93 = m(2) * t178 + t216;
t87 = -mrSges(2,1) * t193 + mrSges(2,3) * t178 - pkin(1) * t95 - t204;
t86 = mrSges(2,2) * t193 - mrSges(2,3) * t177 - pkin(5) * t95 - t199 * t89 + t201 * t91;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t197 * t86 - t195 * t87 - qJ(1) * (t108 * t197 + t195 * t93), t86, t91, t92, t104, -qJD(4) * t152 + t207; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t195 * t86 + t197 * t87 + qJ(1) * (-t108 * t195 + t197 * t93), t87, t89, t96, t102, -t150 * t223 + t209; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t210, t210, t204, t206, t230, -t154 * t222 - t208;];
m_new = t1;
