% Calculate vector of cutting torques with Newton-Euler for
% S4RRRR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:20
% EndTime: 2019-12-31 17:24:25
% DurationCPUTime: 3.24s
% Computational Cost: add. (38401->238), mult. (83185->306), div. (0->0), fcn. (54825->8), ass. (0->97)
t201 = cos(qJ(2));
t198 = sin(qJ(1));
t202 = cos(qJ(1));
t185 = -t202 * g(1) - t198 * g(2);
t203 = qJD(1) ^ 2;
t173 = -t203 * pkin(1) + qJDD(1) * pkin(5) + t185;
t197 = sin(qJ(2));
t219 = t197 * t173;
t160 = -t201 * g(3) - t219;
t161 = -t197 * g(3) + t201 * t173;
t168 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t197 + Ifges(3,2) * t201) * qJD(1);
t169 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t197 + Ifges(3,4) * t201) * qJD(1);
t216 = qJD(1) * qJD(2);
t178 = t197 * qJDD(1) + t201 * t216;
t179 = t201 * qJDD(1) - t197 * t216;
t220 = pkin(2) * t203;
t143 = qJDD(2) * pkin(2) - t178 * pkin(6) - t219 + (pkin(6) * t216 + t197 * t220 - g(3)) * t201;
t218 = qJD(1) * t197;
t183 = qJD(2) * pkin(2) - pkin(6) * t218;
t194 = t201 ^ 2;
t144 = t179 * pkin(6) - qJD(2) * t183 - t194 * t220 + t161;
t196 = sin(qJ(3));
t200 = cos(qJ(3));
t131 = t200 * t143 - t196 * t144;
t170 = (-t196 * t197 + t200 * t201) * qJD(1);
t149 = t170 * qJD(3) + t200 * t178 + t196 * t179;
t171 = (t196 * t201 + t197 * t200) * qJD(1);
t191 = qJDD(2) + qJDD(3);
t192 = qJD(2) + qJD(3);
t119 = (t170 * t192 - t149) * pkin(7) + (t170 * t171 + t191) * pkin(3) + t131;
t132 = t196 * t143 + t200 * t144;
t148 = -t171 * qJD(3) - t196 * t178 + t200 * t179;
t164 = t192 * pkin(3) - t171 * pkin(7);
t166 = t170 ^ 2;
t120 = -t166 * pkin(3) + t148 * pkin(7) - t192 * t164 + t132;
t195 = sin(qJ(4));
t199 = cos(qJ(4));
t117 = t199 * t119 - t195 * t120;
t157 = t199 * t170 - t195 * t171;
t128 = t157 * qJD(4) + t195 * t148 + t199 * t149;
t158 = t195 * t170 + t199 * t171;
t138 = -t157 * mrSges(5,1) + t158 * mrSges(5,2);
t189 = qJD(4) + t192;
t151 = -t189 * mrSges(5,2) + t157 * mrSges(5,3);
t188 = qJDD(4) + t191;
t114 = m(5) * t117 + t188 * mrSges(5,1) - t128 * mrSges(5,3) - t158 * t138 + t189 * t151;
t118 = t195 * t119 + t199 * t120;
t127 = -t158 * qJD(4) + t199 * t148 - t195 * t149;
t152 = t189 * mrSges(5,1) - t158 * mrSges(5,3);
t115 = m(5) * t118 - t188 * mrSges(5,2) + t127 * mrSges(5,3) + t157 * t138 - t189 * t152;
t106 = t199 * t114 + t195 * t115;
t154 = Ifges(4,4) * t171 + Ifges(4,2) * t170 + Ifges(4,6) * t192;
t155 = Ifges(4,1) * t171 + Ifges(4,4) * t170 + Ifges(4,5) * t192;
t134 = Ifges(5,4) * t158 + Ifges(5,2) * t157 + Ifges(5,6) * t189;
t135 = Ifges(5,1) * t158 + Ifges(5,4) * t157 + Ifges(5,5) * t189;
t208 = -mrSges(5,1) * t117 + mrSges(5,2) * t118 - Ifges(5,5) * t128 - Ifges(5,6) * t127 - Ifges(5,3) * t188 - t158 * t134 + t157 * t135;
t206 = -mrSges(4,1) * t131 + mrSges(4,2) * t132 - Ifges(4,5) * t149 - Ifges(4,6) * t148 - Ifges(4,3) * t191 - pkin(3) * t106 - t171 * t154 + t170 * t155 + t208;
t159 = -t170 * mrSges(4,1) + t171 * mrSges(4,2);
t162 = -t192 * mrSges(4,2) + t170 * mrSges(4,3);
t103 = m(4) * t131 + t191 * mrSges(4,1) - t149 * mrSges(4,3) - t171 * t159 + t192 * t162 + t106;
t163 = t192 * mrSges(4,1) - t171 * mrSges(4,3);
t213 = -t195 * t114 + t199 * t115;
t104 = m(4) * t132 - t191 * mrSges(4,2) + t148 * mrSges(4,3) + t170 * t159 - t192 * t163 + t213;
t99 = t200 * t103 + t196 * t104;
t221 = mrSges(3,1) * t160 - mrSges(3,2) * t161 + Ifges(3,5) * t178 + Ifges(3,6) * t179 + Ifges(3,3) * qJDD(2) + pkin(2) * t99 + (t197 * t168 - t201 * t169) * qJD(1) - t206;
t217 = qJD(1) * t201;
t184 = t198 * g(1) - t202 * g(2);
t177 = (-mrSges(3,1) * t201 + mrSges(3,2) * t197) * qJD(1);
t182 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t217;
t97 = m(3) * t160 + qJDD(2) * mrSges(3,1) - t178 * mrSges(3,3) + qJD(2) * t182 - t177 * t218 + t99;
t181 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t218;
t214 = -t196 * t103 + t200 * t104;
t98 = m(3) * t161 - qJDD(2) * mrSges(3,2) + t179 * mrSges(3,3) - qJD(2) * t181 + t177 * t217 + t214;
t215 = -t197 * t97 + t201 * t98;
t210 = -qJDD(1) * pkin(1) - t184;
t150 = -t179 * pkin(2) + t183 * t218 + (-pkin(6) * t194 - pkin(5)) * t203 + t210;
t122 = -t148 * pkin(3) - t166 * pkin(7) + t171 * t164 + t150;
t212 = m(5) * t122 - t127 * mrSges(5,1) + t128 * mrSges(5,2) - t157 * t151 + t158 * t152;
t172 = -t203 * pkin(5) + t210;
t207 = m(4) * t150 - t148 * mrSges(4,1) + t149 * mrSges(4,2) - t170 * t162 + t171 * t163 + t212;
t205 = -m(3) * t172 + t179 * mrSges(3,1) - t178 * mrSges(3,2) - t181 * t218 + t182 * t217 - t207;
t167 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t197 + Ifges(3,6) * t201) * qJD(1);
t133 = Ifges(5,5) * t158 + Ifges(5,6) * t157 + Ifges(5,3) * t189;
t107 = -mrSges(5,1) * t122 + mrSges(5,3) * t118 + Ifges(5,4) * t128 + Ifges(5,2) * t127 + Ifges(5,6) * t188 - t158 * t133 + t189 * t135;
t108 = mrSges(5,2) * t122 - mrSges(5,3) * t117 + Ifges(5,1) * t128 + Ifges(5,4) * t127 + Ifges(5,5) * t188 + t157 * t133 - t189 * t134;
t153 = Ifges(4,5) * t171 + Ifges(4,6) * t170 + Ifges(4,3) * t192;
t94 = -mrSges(4,1) * t150 + mrSges(4,3) * t132 + Ifges(4,4) * t149 + Ifges(4,2) * t148 + Ifges(4,6) * t191 - pkin(3) * t212 + pkin(7) * t213 + t199 * t107 + t195 * t108 - t171 * t153 + t192 * t155;
t95 = mrSges(4,2) * t150 - mrSges(4,3) * t131 + Ifges(4,1) * t149 + Ifges(4,4) * t148 + Ifges(4,5) * t191 - pkin(7) * t106 - t195 * t107 + t199 * t108 + t170 * t153 - t192 * t154;
t88 = -mrSges(3,1) * t172 + mrSges(3,3) * t161 + Ifges(3,4) * t178 + Ifges(3,2) * t179 + Ifges(3,6) * qJDD(2) - pkin(2) * t207 + pkin(6) * t214 + qJD(2) * t169 - t167 * t218 + t196 * t95 + t200 * t94;
t90 = mrSges(3,2) * t172 - mrSges(3,3) * t160 + Ifges(3,1) * t178 + Ifges(3,4) * t179 + Ifges(3,5) * qJDD(2) - pkin(6) * t99 - qJD(2) * t168 + t167 * t217 - t196 * t94 + t200 * t95;
t209 = mrSges(2,1) * t184 - mrSges(2,2) * t185 + Ifges(2,3) * qJDD(1) + pkin(1) * t205 + pkin(5) * t215 + t197 * t90 + t201 * t88;
t109 = m(2) * t184 + qJDD(1) * mrSges(2,1) - t203 * mrSges(2,2) + t205;
t93 = t197 * t98 + t201 * t97;
t91 = m(2) * t185 - t203 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t215;
t86 = mrSges(2,1) * g(3) + mrSges(2,3) * t185 + t203 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t93 - t221;
t85 = -mrSges(2,2) * g(3) - mrSges(2,3) * t184 + Ifges(2,5) * qJDD(1) - t203 * Ifges(2,6) - pkin(5) * t93 - t197 * t88 + t201 * t90;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t202 * t85 - t198 * t86 - pkin(4) * (t202 * t109 + t198 * t91), t85, t90, t95, t108; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t198 * t85 + t202 * t86 + pkin(4) * (-t198 * t109 + t202 * t91), t86, t88, t94, t107; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209, t209, t221, -t206, -t208;];
m_new = t1;
