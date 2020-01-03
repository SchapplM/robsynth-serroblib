% Calculate vector of cutting torques with Newton-Euler for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:37
% EndTime: 2019-12-31 16:41:38
% DurationCPUTime: 0.95s
% Computational Cost: add. (7810->171), mult. (16489->211), div. (0->0), fcn. (9176->6), ass. (0->79)
t164 = qJD(1) ^ 2;
t158 = sin(pkin(6));
t149 = t158 ^ 2;
t159 = cos(pkin(6));
t191 = t159 ^ 2 + t149;
t183 = t191 * mrSges(4,3);
t199 = t164 * t183;
t161 = sin(qJ(1));
t163 = cos(qJ(1));
t139 = t161 * g(1) - t163 * g(2);
t176 = -t164 * qJ(2) + qJDD(2) - t139;
t193 = -pkin(1) - qJ(3);
t198 = -(2 * qJD(1) * qJD(3)) + t193 * qJDD(1) + t176;
t140 = -t163 * g(1) - t161 * g(2);
t197 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t140;
t196 = pkin(3) * t164;
t195 = mrSges(2,1) - mrSges(3,2);
t194 = -Ifges(2,6) + Ifges(3,5);
t160 = sin(qJ(4));
t162 = cos(qJ(4));
t116 = t158 * g(3) + t198 * t159;
t103 = (-pkin(5) * qJDD(1) - t158 * t196) * t159 + t116;
t117 = -t159 * g(3) + t198 * t158;
t188 = qJDD(1) * t158;
t104 = -pkin(5) * t188 - t149 * t196 + t117;
t100 = t162 * t103 - t160 * t104;
t180 = -t158 * t162 - t159 * t160;
t132 = t180 * qJD(1);
t179 = -t158 * t160 + t159 * t162;
t133 = t179 * qJD(1);
t112 = -t132 * mrSges(5,1) + t133 * mrSges(5,2);
t119 = t132 * qJD(4) + t179 * qJDD(1);
t124 = -qJD(4) * mrSges(5,2) + t132 * mrSges(5,3);
t97 = m(5) * t100 + qJDD(4) * mrSges(5,1) - t119 * mrSges(5,3) + qJD(4) * t124 - t133 * t112;
t101 = t160 * t103 + t162 * t104;
t118 = -t133 * qJD(4) + t180 * qJDD(1);
t125 = qJD(4) * mrSges(5,1) - t133 * mrSges(5,3);
t98 = m(5) * t101 - qJDD(4) * mrSges(5,2) + t118 * mrSges(5,3) - qJD(4) * t125 + t132 * t112;
t88 = t160 * t98 + t162 * t97;
t192 = Ifges(4,6) * t158;
t190 = t164 * (Ifges(4,5) * t159 - t192);
t187 = qJDD(1) * t159;
t185 = Ifges(3,4) + t192;
t178 = -qJDD(1) * mrSges(4,3) - t164 * (mrSges(4,1) * t158 + mrSges(4,2) * t159);
t85 = m(4) * t116 + t178 * t159 + t88;
t184 = -t160 * t97 + t162 * t98;
t86 = m(4) * t117 + t178 * t158 + t184;
t83 = -t158 * t85 + t159 * t86;
t182 = Ifges(4,1) * t159 - Ifges(4,4) * t158;
t181 = Ifges(4,4) * t159 - Ifges(4,2) * t158;
t82 = t158 * t86 + t159 * t85;
t175 = qJDD(3) + t197;
t106 = pkin(3) * t188 + (-t191 * pkin(5) + t193) * t164 + t175;
t174 = -m(5) * t106 + t118 * mrSges(5,1) - t119 * mrSges(5,2) + t132 * t124 - t133 * t125;
t131 = -qJDD(1) * pkin(1) + t176;
t173 = -m(3) * t131 + t164 * mrSges(3,3) - t82;
t108 = Ifges(5,4) * t133 + Ifges(5,2) * t132 + Ifges(5,6) * qJD(4);
t109 = Ifges(5,1) * t133 + Ifges(5,4) * t132 + Ifges(5,5) * qJD(4);
t172 = -mrSges(5,1) * t100 + mrSges(5,2) * t101 - Ifges(5,5) * t119 - Ifges(5,6) * t118 - Ifges(5,3) * qJDD(4) - t133 * t108 + t132 * t109;
t127 = t164 * pkin(1) - t197;
t123 = t193 * t164 + t175;
t107 = Ifges(5,5) * t133 + Ifges(5,6) * t132 + Ifges(5,3) * qJD(4);
t89 = -mrSges(5,1) * t106 + mrSges(5,3) * t101 + Ifges(5,4) * t119 + Ifges(5,2) * t118 + Ifges(5,6) * qJDD(4) + qJD(4) * t109 - t133 * t107;
t90 = mrSges(5,2) * t106 - mrSges(5,3) * t100 + Ifges(5,1) * t119 + Ifges(5,4) * t118 + Ifges(5,5) * qJDD(4) - qJD(4) * t108 + t132 * t107;
t76 = -mrSges(4,1) * t123 + mrSges(4,3) * t117 + pkin(3) * t174 + pkin(5) * t184 + t181 * qJDD(1) - t159 * t190 + t160 * t90 + t162 * t89;
t78 = mrSges(4,2) * t123 - mrSges(4,3) * t116 - pkin(5) * t88 + t182 * qJDD(1) - t158 * t190 - t160 * t89 + t162 * t90;
t171 = mrSges(3,2) * t131 - mrSges(3,3) * t127 + Ifges(3,1) * qJDD(1) - qJ(3) * t82 - t158 * t76 + t159 * t78;
t170 = -m(4) * t123 - mrSges(4,1) * t188 - mrSges(4,2) * t187 + t174;
t169 = -mrSges(3,1) * t127 - pkin(2) * (t170 + t199) - qJ(3) * t83 - t158 * t78 - t159 * t76;
t168 = -m(3) * t127 + t164 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t170;
t167 = -mrSges(2,2) * t140 + mrSges(2,1) * t139 + Ifges(2,3) * qJDD(1) + t171 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t173) + qJ(2) * (t168 - t199);
t166 = -mrSges(4,1) * t116 + mrSges(4,2) * t117 - Ifges(4,5) * t187 - pkin(3) * t88 + t172 + (-t158 * t182 - t159 * t181) * t164;
t165 = -mrSges(3,1) * t131 - pkin(2) * t82 + t166;
t91 = t168 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t183) * t164 + m(2) * t140;
t81 = -m(3) * g(3) + t83;
t79 = m(2) * t139 - t164 * mrSges(2,2) + t195 * qJDD(1) + t173;
t75 = -t165 + t194 * t164 + (Ifges(2,5) - t185) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t139 - qJ(2) * t81;
t74 = mrSges(2,3) * t140 - pkin(1) * t81 + (-Ifges(3,4) + Ifges(2,5)) * t164 - t194 * qJDD(1) + t195 * g(3) + t169;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t163 * t75 - t161 * t74 - pkin(4) * (t161 * t91 + t163 * t79), t75, t171, t78, t90; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t161 * t75 + t163 * t74 + pkin(4) * (-t161 * t79 + t163 * t91), t74, -mrSges(3,3) * g(3) - t164 * Ifges(3,5) + t185 * qJDD(1) + t165, t76, t89; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t167, t167, mrSges(3,2) * g(3) + t164 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t169, -Ifges(4,6) * t188 - t166, -t172;];
m_new = t1;
