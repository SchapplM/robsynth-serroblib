% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:34
% EndTime: 2019-12-05 16:21:40
% DurationCPUTime: 2.05s
% Computational Cost: add. (4731->182), mult. (10659->268), div. (0->0), fcn. (11751->8), ass. (0->104)
t131 = sin(pkin(9));
t133 = sin(qJ(3));
t136 = cos(qJ(3));
t167 = cos(pkin(9));
t149 = t167 * t136;
t109 = -t131 * t133 + t149;
t110 = -t131 * t136 - t133 * t167;
t132 = sin(qJ(5));
t135 = cos(qJ(5));
t147 = t109 * t135 + t110 * t132;
t176 = -qJ(4) - pkin(6);
t116 = t176 * t133;
t117 = t176 * t136;
t216 = t116 * t167 + t117 * t131;
t225 = pkin(7) * t110 + t216;
t92 = t116 * t131 - t117 * t167;
t66 = -pkin(7) * t109 - t92;
t231 = t132 * t66 + t135 * t225;
t37 = t132 * t225 - t135 * t66;
t87 = t109 * t132 - t110 * t135;
t6 = -t37 * mrSges(6,1) - t231 * mrSges(6,2) + Ifges(6,5) * t147 - Ifges(6,6) * t87;
t245 = t6 * qJD(5);
t212 = m(6) / 0.2e1;
t194 = t133 * pkin(3);
t98 = -pkin(4) * t110 + t194;
t243 = t98 * t212;
t226 = t87 * mrSges(6,3);
t238 = -mrSges(5,1) * t110 + t109 * mrSges(5,2);
t201 = Ifges(6,4) * t87;
t204 = t87 / 0.2e1;
t205 = -t87 / 0.2e1;
t222 = t147 / 0.2e1;
t227 = t87 * mrSges(6,1);
t232 = mrSges(6,2) * t147 + t227;
t125 = -pkin(3) * t136 - pkin(2);
t97 = -pkin(4) * t109 + t125;
t4 = t97 * t232 + (0.2e1 * Ifges(6,4) * t147 + (Ifges(6,1) - Ifges(6,2)) * t87) * t222 + (Ifges(6,2) * t147 + t201) * t205 + (Ifges(6,1) * t147 - t201) * t204;
t237 = t110 ^ 2;
t156 = t227 / 0.2e1;
t134 = sin(qJ(2));
t100 = t110 * t134;
t166 = t133 * t134;
t99 = t131 * t166 - t134 * t149;
t151 = t100 * t135 + t132 * t99;
t64 = t100 * t132 - t135 * t99;
t13 = -t64 * mrSges(6,1) - mrSges(6,2) * t151;
t234 = t13 * qJD(5);
t218 = t147 * mrSges(6,3);
t215 = t232 + t238;
t169 = t136 * mrSges(4,1);
t214 = mrSges(4,2) * t133 - t169;
t137 = cos(qJ(2));
t101 = t110 * t137;
t102 = t137 * t109;
t62 = t101 * t135 - t102 * t132;
t65 = t101 * t132 + t102 * t135;
t175 = t62 * mrSges(6,1) / 0.2e1 - t65 * mrSges(6,2) / 0.2e1;
t213 = m(5) / 0.2e1;
t211 = m(5) * pkin(3);
t202 = -t137 / 0.2e1;
t196 = pkin(3) * t131;
t172 = t109 * mrSges(5,3);
t171 = t110 * mrSges(5,3);
t168 = t136 * mrSges(4,2);
t165 = t133 * t137;
t164 = t137 * t134;
t158 = t133 ^ 2 + t136 ^ 2;
t15 = m(6) * (t151 * t62 + t64 * t65 - t164) + m(5) * (t100 * t101 - t102 * t99 - t164) + m(4) * (-0.1e1 + t158) * t164;
t163 = t15 * qJD(1);
t154 = t167 * pkin(3);
t124 = t154 + pkin(4);
t103 = t124 * t135 - t132 * t196;
t104 = t124 * t132 + t135 * t196;
t157 = t211 / 0.2e1;
t140 = (-t103 * t87 + t104 * t147) * t212 + (t109 * t131 + t110 * t167) * t157;
t143 = t133 * t157 + t243;
t16 = -t140 + t143 + t215;
t162 = t16 * qJD(2);
t25 = 0.2e1 * mrSges(6,2) * t222 + 0.2e1 * t156;
t160 = t25 * qJD(2);
t28 = -mrSges(6,1) * t104 - mrSges(6,2) * t103;
t159 = t28 * qJD(5);
t155 = t232 * t202;
t47 = -mrSges(6,1) * t147 + mrSges(6,2) * t87;
t88 = -mrSges(5,1) * t109 - mrSges(5,2) * t110;
t1 = (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t133 + pkin(3) * t88) * t133 - Ifges(5,4) * t237 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t136 + (Ifges(4,1) - Ifges(4,2)) * t133) * t136 + (Ifges(5,4) * t109 + (-Ifges(5,1) + Ifges(5,2)) * t110) * t109 + (m(5) * t194 + t238) * t125 + t4 + (m(6) * t97 + t47) * t98;
t138 = -pkin(3) * t165 * t213 - t137 * t243 + (mrSges(4,1) * t133 + t168 + t215) * t202;
t139 = (t103 * t62 + t104 * t65) * t212 + t101 * mrSges(5,1) / 0.2e1 - t102 * mrSges(5,2) / 0.2e1 + (t101 * t167 + t102 * t131) * t157 - mrSges(4,1) * t165 / 0.2e1 + t168 * t202 + t175;
t2 = -t138 + t139;
t146 = -t2 * qJD(1) + t1 * qJD(2);
t7 = t155 - t175;
t145 = qJD(1) * t7 + qJD(2) * t4;
t11 = (t147 ^ 2 + t87 ^ 2) * mrSges(6,3) + (t109 ^ 2 + t237) * mrSges(5,3) + m(6) * (t147 * t37 - t231 * t87) + m(5) * (t109 * t92 + t110 * t216);
t141 = (t147 * t64 - t151 * t87) * t212 + (t100 * t110 - t109 * t99) * t213;
t20 = (-m(6) / 0.2e1 - m(5) / 0.2e1) * t134 + t141;
t144 = -qJD(1) * t20 - qJD(2) * t11;
t5 = (t205 + t204) * Ifges(6,6) + (t222 - t147 / 0.2e1) * Ifges(6,5);
t142 = qJD(2) * t5 + qJD(3) * t28;
t26 = t156 - t227 / 0.2e1;
t22 = t140 + t143;
t19 = t141 + (m(5) + m(6)) * t134 / 0.2e1;
t8 = t155 + t175;
t3 = t138 + t139;
t9 = [t15 * qJD(2), t3 * qJD(3) + t19 * qJD(4) + t8 * qJD(5) + t163 + (t101 * t171 + t102 * t172 - t62 * t226 + t65 * t218 + (mrSges(4,3) * t158 - mrSges(3,2)) * t137 + (-mrSges(3,1) + t47 + t88 + t214) * t134 + 0.2e1 * (t134 * t97 + t231 * t62 + t37 * t65) * t212 + 0.2e1 * (t101 * t216 + t102 * t92 + t125 * t134) * t213 + m(4) * (pkin(6) * t137 * t158 - t134 * pkin(2))) * qJD(2), t3 * qJD(2) + (mrSges(4,2) * t166 - t134 * t169 + m(6) * (-t103 * t64 + t104 * t151) + (t100 * t131 + t167 * t99) * t211 - t100 * mrSges(5,2) + t99 * mrSges(5,1) + t13) * qJD(3) + t234, t19 * qJD(2), qJD(2) * t8 + qJD(3) * t13 + t234; -qJD(3) * t2 + qJD(4) * t20 + qJD(5) * t7 - t163, qJD(3) * t1 + qJD(4) * t11 + qJD(5) * t4, (-t104 * t226 - t103 * t218 + t171 * t196 - t154 * t172 + m(6) * (-t103 * t37 + t104 * t231) - t92 * mrSges(5,1) - t216 * mrSges(5,2) + (t131 * t216 - t167 * t92) * t211 + Ifges(5,5) * t109 + Ifges(5,6) * t110 + Ifges(4,5) * t136 - Ifges(4,6) * t133 + t214 * pkin(6) + t6) * qJD(3) + t22 * qJD(4) + t245 + t146, qJD(3) * t22 + qJD(5) * t26 - t144, t6 * qJD(3) + t26 * qJD(4) + t145 + t245; qJD(2) * t2, -qJD(4) * t16 + qJD(5) * t5 - t146, t159, -t162, t142 + t159; -t20 * qJD(2), qJD(3) * t16 + qJD(5) * t25 + t144, t162, 0, t160; -t7 * qJD(2), -qJD(3) * t5 - qJD(4) * t25 - t145, -t142, -t160, 0;];
Cq = t9;
