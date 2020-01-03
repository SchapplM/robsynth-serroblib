% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR6
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:26
% EndTime: 2019-12-31 17:04:31
% DurationCPUTime: 3.07s
% Computational Cost: add. (34362->237), mult. (79502->308), div. (0->0), fcn. (51357->8), ass. (0->95)
t201 = cos(qJ(2));
t199 = sin(qJ(1));
t202 = cos(qJ(1));
t187 = -t202 * g(1) - t199 * g(2);
t203 = qJD(1) ^ 2;
t175 = -t203 * pkin(1) + qJDD(1) * pkin(5) + t187;
t198 = sin(qJ(2));
t219 = t198 * t175;
t160 = -t201 * g(3) - t219;
t161 = -t198 * g(3) + t201 * t175;
t172 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t198 + Ifges(3,2) * t201) * qJD(1);
t173 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t198 + Ifges(3,4) * t201) * qJD(1);
t216 = qJD(1) * qJD(2);
t180 = t198 * qJDD(1) + t201 * t216;
t181 = t201 * qJDD(1) - t198 * t216;
t220 = pkin(2) * t203;
t143 = qJDD(2) * pkin(2) - t180 * qJ(3) - t219 + (qJ(3) * t216 + t198 * t220 - g(3)) * t201;
t218 = qJD(1) * t198;
t183 = qJD(2) * pkin(2) - qJ(3) * t218;
t194 = t201 ^ 2;
t144 = t181 * qJ(3) - qJD(2) * t183 - t194 * t220 + t161;
t195 = sin(pkin(7));
t196 = cos(pkin(7));
t170 = (t195 * t201 + t196 * t198) * qJD(1);
t124 = -0.2e1 * qJD(3) * t170 + t196 * t143 - t195 * t144;
t159 = t196 * t180 + t195 * t181;
t169 = (-t195 * t198 + t196 * t201) * qJD(1);
t119 = (qJD(2) * t169 - t159) * pkin(6) + (t169 * t170 + qJDD(2)) * pkin(3) + t124;
t125 = 0.2e1 * qJD(3) * t169 + t195 * t143 + t196 * t144;
t158 = -t195 * t180 + t196 * t181;
t164 = qJD(2) * pkin(3) - t170 * pkin(6);
t168 = t169 ^ 2;
t120 = -t168 * pkin(3) + t158 * pkin(6) - qJD(2) * t164 + t125;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t117 = t200 * t119 - t197 * t120;
t151 = t200 * t169 - t197 * t170;
t132 = t151 * qJD(4) + t197 * t158 + t200 * t159;
t152 = t197 * t169 + t200 * t170;
t138 = -t151 * mrSges(5,1) + t152 * mrSges(5,2);
t191 = qJD(2) + qJD(4);
t146 = -t191 * mrSges(5,2) + t151 * mrSges(5,3);
t190 = qJDD(2) + qJDD(4);
t114 = m(5) * t117 + t190 * mrSges(5,1) - t132 * mrSges(5,3) - t152 * t138 + t191 * t146;
t118 = t197 * t119 + t200 * t120;
t131 = -t152 * qJD(4) + t200 * t158 - t197 * t159;
t147 = t191 * mrSges(5,1) - t152 * mrSges(5,3);
t115 = m(5) * t118 - t190 * mrSges(5,2) + t131 * mrSges(5,3) + t151 * t138 - t191 * t147;
t106 = t200 * t114 + t197 * t115;
t149 = Ifges(4,4) * t170 + Ifges(4,2) * t169 + Ifges(4,6) * qJD(2);
t150 = Ifges(4,1) * t170 + Ifges(4,4) * t169 + Ifges(4,5) * qJD(2);
t134 = Ifges(5,4) * t152 + Ifges(5,2) * t151 + Ifges(5,6) * t191;
t135 = Ifges(5,1) * t152 + Ifges(5,4) * t151 + Ifges(5,5) * t191;
t208 = -mrSges(5,1) * t117 + mrSges(5,2) * t118 - Ifges(5,5) * t132 - Ifges(5,6) * t131 - Ifges(5,3) * t190 - t152 * t134 + t151 * t135;
t206 = -mrSges(4,1) * t124 + mrSges(4,2) * t125 - Ifges(4,5) * t159 - Ifges(4,6) * t158 - Ifges(4,3) * qJDD(2) - pkin(3) * t106 - t170 * t149 + t169 * t150 + t208;
t154 = -t169 * mrSges(4,1) + t170 * mrSges(4,2);
t162 = -qJD(2) * mrSges(4,2) + t169 * mrSges(4,3);
t103 = m(4) * t124 + qJDD(2) * mrSges(4,1) - t159 * mrSges(4,3) + qJD(2) * t162 - t170 * t154 + t106;
t163 = qJD(2) * mrSges(4,1) - t170 * mrSges(4,3);
t213 = -t197 * t114 + t200 * t115;
t104 = m(4) * t125 - qJDD(2) * mrSges(4,2) + t158 * mrSges(4,3) - qJD(2) * t163 + t169 * t154 + t213;
t99 = t196 * t103 + t195 * t104;
t221 = mrSges(3,1) * t160 - mrSges(3,2) * t161 + Ifges(3,5) * t180 + Ifges(3,6) * t181 + Ifges(3,3) * qJDD(2) + pkin(2) * t99 + (t198 * t172 - t201 * t173) * qJD(1) - t206;
t217 = qJD(1) * t201;
t186 = t199 * g(1) - t202 * g(2);
t179 = (-mrSges(3,1) * t201 + mrSges(3,2) * t198) * qJD(1);
t185 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t217;
t97 = m(3) * t160 + qJDD(2) * mrSges(3,1) - t180 * mrSges(3,3) + qJD(2) * t185 - t179 * t218 + t99;
t184 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t218;
t214 = -t195 * t103 + t196 * t104;
t98 = m(3) * t161 - qJDD(2) * mrSges(3,2) + t181 * mrSges(3,3) - qJD(2) * t184 + t179 * t217 + t214;
t215 = -t198 * t97 + t201 * t98;
t210 = -qJDD(1) * pkin(1) - t186;
t145 = -t181 * pkin(2) + qJDD(3) + t183 * t218 + (-qJ(3) * t194 - pkin(5)) * t203 + t210;
t122 = -t158 * pkin(3) - t168 * pkin(6) + t170 * t164 + t145;
t212 = m(5) * t122 - t131 * mrSges(5,1) + t132 * mrSges(5,2) - t151 * t146 + t152 * t147;
t174 = -t203 * pkin(5) + t210;
t207 = m(4) * t145 - t158 * mrSges(4,1) + t159 * mrSges(4,2) - t169 * t162 + t170 * t163 + t212;
t205 = -m(3) * t174 + t181 * mrSges(3,1) - t180 * mrSges(3,2) - t184 * t218 + t185 * t217 - t207;
t171 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t198 + Ifges(3,6) * t201) * qJD(1);
t133 = Ifges(5,5) * t152 + Ifges(5,6) * t151 + Ifges(5,3) * t191;
t107 = -mrSges(5,1) * t122 + mrSges(5,3) * t118 + Ifges(5,4) * t132 + Ifges(5,2) * t131 + Ifges(5,6) * t190 - t152 * t133 + t191 * t135;
t108 = mrSges(5,2) * t122 - mrSges(5,3) * t117 + Ifges(5,1) * t132 + Ifges(5,4) * t131 + Ifges(5,5) * t190 + t151 * t133 - t191 * t134;
t148 = Ifges(4,5) * t170 + Ifges(4,6) * t169 + Ifges(4,3) * qJD(2);
t94 = -mrSges(4,1) * t145 + mrSges(4,3) * t125 + Ifges(4,4) * t159 + Ifges(4,2) * t158 + Ifges(4,6) * qJDD(2) - pkin(3) * t212 + pkin(6) * t213 + qJD(2) * t150 + t200 * t107 + t197 * t108 - t170 * t148;
t95 = mrSges(4,2) * t145 - mrSges(4,3) * t124 + Ifges(4,1) * t159 + Ifges(4,4) * t158 + Ifges(4,5) * qJDD(2) - pkin(6) * t106 - qJD(2) * t149 - t197 * t107 + t200 * t108 + t169 * t148;
t88 = -mrSges(3,1) * t174 + mrSges(3,3) * t161 + Ifges(3,4) * t180 + Ifges(3,2) * t181 + Ifges(3,6) * qJDD(2) - pkin(2) * t207 + qJ(3) * t214 + qJD(2) * t173 - t171 * t218 + t195 * t95 + t196 * t94;
t90 = mrSges(3,2) * t174 - mrSges(3,3) * t160 + Ifges(3,1) * t180 + Ifges(3,4) * t181 + Ifges(3,5) * qJDD(2) - qJ(3) * t99 - qJD(2) * t172 + t171 * t217 - t195 * t94 + t196 * t95;
t209 = mrSges(2,1) * t186 - mrSges(2,2) * t187 + Ifges(2,3) * qJDD(1) + pkin(1) * t205 + pkin(5) * t215 + t198 * t90 + t201 * t88;
t109 = m(2) * t186 + qJDD(1) * mrSges(2,1) - t203 * mrSges(2,2) + t205;
t93 = t198 * t98 + t201 * t97;
t91 = m(2) * t187 - t203 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t215;
t86 = mrSges(2,1) * g(3) + mrSges(2,3) * t187 + t203 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t93 - t221;
t85 = -mrSges(2,2) * g(3) - mrSges(2,3) * t186 + Ifges(2,5) * qJDD(1) - t203 * Ifges(2,6) - pkin(5) * t93 - t198 * t88 + t201 * t90;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t202 * t85 - t199 * t86 - pkin(4) * (t202 * t109 + t199 * t91), t85, t90, t95, t108; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t199 * t85 + t202 * t86 + pkin(4) * (-t199 * t109 + t202 * t91), t86, t88, t94, t107; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209, t209, t221, -t206, -t208;];
m_new = t1;
