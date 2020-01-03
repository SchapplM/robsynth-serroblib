% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:23
% EndTime: 2019-12-31 20:03:27
% DurationCPUTime: 1.59s
% Computational Cost: add. (3163->199), mult. (5887->246), div. (0->0), fcn. (5356->4), ass. (0->103)
t124 = sin(qJ(2));
t180 = pkin(6) * t124;
t109 = -t124 * pkin(7) + t180;
t187 = cos(qJ(2));
t153 = t187 * pkin(6);
t110 = -t187 * pkin(7) + t153;
t125 = cos(qJ(4));
t186 = sin(qJ(4));
t134 = t186 * t109 + t125 * t110;
t98 = -t124 * t186 - t187 * t125;
t131 = t98 * qJ(5) + t134;
t63 = -t125 * t109 + t186 * t110;
t99 = t124 * t125 - t187 * t186;
t141 = -qJ(5) * t99 - t63;
t224 = -t134 * mrSges(5,1) - t131 * mrSges(6,1) + t63 * mrSges(5,2) - t141 * mrSges(6,2);
t117 = t124 * qJ(3);
t200 = -t187 * pkin(2) - t117;
t106 = -pkin(1) + t200;
t87 = t187 * pkin(3) - t106;
t49 = -t98 * pkin(4) + t87;
t211 = m(6) * t49 - mrSges(6,1) * t98 + mrSges(6,2) * t99;
t223 = m(5) * t87 - mrSges(5,1) * t98 + mrSges(5,2) * t99;
t181 = Ifges(6,6) * t99;
t182 = Ifges(5,6) * t99;
t221 = -t181 - t182 + t224;
t216 = t125 * t131 - t186 * t141;
t215 = m(5) * (t125 * t134 + t186 * t63);
t214 = mrSges(5,1) + mrSges(6,1);
t126 = -pkin(2) - pkin(3);
t104 = -t186 * qJ(3) + t125 * t126;
t103 = -pkin(4) + t104;
t210 = -t103 + t104;
t206 = -Ifges(4,5) + Ifges(3,4);
t89 = t98 * mrSges(6,2);
t203 = t99 * mrSges(6,1) + t89;
t107 = -t187 * mrSges(4,1) - t124 * mrSges(4,3);
t105 = t125 * qJ(3) + t186 * t126;
t159 = t105 * t125;
t198 = -t104 * t186 + t159;
t197 = m(4) * t106 + t107;
t170 = mrSges(5,2) + mrSges(6,2);
t196 = t170 * t104 + t214 * t105;
t130 = t170 * t125 + t214 * t186;
t194 = m(6) / 0.2e1;
t193 = m(6) * pkin(4);
t192 = t98 / 0.2e1;
t189 = pkin(4) * t131;
t188 = pkin(4) * t99;
t185 = mrSges(6,3) * t98;
t94 = Ifges(5,4) * t99;
t92 = Ifges(6,4) * t99;
t184 = Ifges(5,5) * t98;
t183 = Ifges(6,5) * t98;
t118 = t187 * qJ(3);
t45 = Ifges(6,2) * t98 + t92;
t46 = Ifges(5,2) * t98 + t94;
t91 = Ifges(6,4) * t98;
t47 = Ifges(6,1) * t99 + t91;
t93 = Ifges(5,4) * t98;
t48 = Ifges(5,1) * t99 + t93;
t95 = t126 * t124 + t118;
t54 = t95 - t188;
t1 = (t45 / 0.2e1 + t46 / 0.2e1 - t49 * mrSges(6,1) - t87 * mrSges(5,1) + t92 / 0.2e1 + t94 / 0.2e1) * t99 + (-t47 / 0.2e1 - t48 / 0.2e1 - t49 * mrSges(6,2) - t87 * mrSges(5,2) + (-Ifges(6,4) / 0.2e1 - Ifges(5,4) / 0.2e1) * t98 + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t99) * t98 + (-pkin(1) * mrSges(3,1) + t106 * mrSges(4,1) - t124 * t206) * t124 + t197 * (pkin(2) * t124 - t118) + (-pkin(1) * mrSges(3,2) - mrSges(4,3) * t106 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t124 + t206 * t187) * t187 + t223 * t95 + t211 * t54;
t167 = t1 * qJD(1);
t166 = t103 * t99;
t165 = t105 * t98;
t164 = t105 * t99;
t162 = t125 * t99;
t12 = m(6) * (t131 * t98 - t141 * t99) + (t98 ^ 2 + t99 ^ 2) * mrSges(6,3);
t161 = qJD(1) * t12;
t16 = (-t197 + t211 + t223) * t124;
t160 = qJD(1) * t16;
t15 = 0.2e1 * (-t166 / 0.4e1 + t165 / 0.4e1 - t54 / 0.4e1) * m(6) + t203;
t157 = t15 * qJD(1);
t27 = -t89 + (-mrSges(6,1) - t193) * t99;
t156 = t27 * qJD(1);
t149 = t186 * t98;
t34 = (t149 / 0.2e1 - t162 / 0.2e1 - t124 / 0.2e1) * m(6);
t155 = t34 * qJD(1);
t151 = t187 * mrSges(4,2);
t146 = -t103 / 0.2e1 + t104 / 0.2e1;
t144 = t210 * t131;
t142 = t186 * t193;
t140 = -t103 * t186 + t159;
t2 = t87 * (t99 * mrSges(5,1) + t98 * mrSges(5,2)) + t49 * t203 - (t45 + t46) * t99 / 0.2e1 + (-t92 - t94 + (Ifges(5,1) + Ifges(6,1)) * t98) * t99 / 0.2e1 + (t47 + t48 + t91 + t93 + (-Ifges(5,2) - Ifges(6,2)) * t99) * t192 + t211 * t188;
t138 = t2 * qJD(1);
t13 = m(6) * t210 * t105 + t196;
t3 = (t144 / 0.2e1 - t189 / 0.2e1) * m(6) + (-pkin(4) / 0.2e1 + t146) * t185;
t137 = t3 * qJD(1) + t13 * qJD(2);
t129 = (t140 - t198) * t194;
t17 = t142 / 0.2e1 + t129 + t130;
t136 = t17 * qJD(2);
t21 = m(4) * qJ(3) + m(5) * t198 + m(6) * t140 + mrSges(4,3) + t130;
t127 = t216 * t194 + t215 / 0.2e1;
t128 = -t215 / 0.2e1 - m(6) * t216 / 0.2e1;
t6 = t127 + t128;
t135 = -qJD(1) * t6 - qJD(2) * t21;
t33 = (t149 - t162 + t124) * t194;
t18 = -t142 / 0.2e1 + t129;
t14 = (t165 - t166 + t54) * t194;
t5 = m(4) * t153 + t127 - t128 + t151 + (mrSges(5,3) + mrSges(6,3)) * (0.2e1 * t125 * t192 + t186 * t99);
t4 = -t184 / 0.2e1 - t183 / 0.2e1 + t181 / 0.2e1 + t182 / 0.2e1 + pkin(4) * t185 / 0.2e1 + (Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t99 + (-Ifges(5,5) / 0.2e1 - Ifges(6,5) / 0.2e1 + t146 * mrSges(6,3)) * t98 + (t189 + t144) * t194 - t224;
t7 = [qJD(2) * t1 + qJD(3) * t16 + qJD(4) * t2 + qJD(5) * t12, t5 * qJD(3) + t4 * qJD(4) + t14 * qJD(5) + t167 + (t103 * t185 + mrSges(6,3) * t164 - mrSges(3,1) * t153 - pkin(2) * t151 + m(6) * (t103 * t131 - t105 * t141) + m(5) * (t104 * t134 + t105 * t63) - mrSges(4,2) * t117 + mrSges(3,2) * t180 + t184 + t183 + (m(4) * t200 + t107) * pkin(6) + (Ifges(4,4) + Ifges(3,5)) * t187 + (Ifges(4,6) - Ifges(3,6)) * t124 + (t104 * t98 + t164) * mrSges(5,3) + t221) * qJD(2), qJD(2) * t5 + qJD(5) * t33 + t160, t4 * qJD(2) + t138 + (-t131 * t193 + (-mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t98 + t221) * qJD(4), qJD(2) * t14 + qJD(3) * t33 + t161; qJD(3) * t6 + qJD(4) * t3 + qJD(5) * t15 - t167, qJD(3) * t21 + qJD(4) * t13, qJD(4) * t18 - t135, t18 * qJD(3) + (-t105 * t193 - t196) * qJD(4) + t137, t157; -qJD(2) * t6 + qJD(5) * t34 - t160, qJD(4) * t17 + t135, 0, (-t130 - t142) * qJD(4) + t136, t155; -qJD(2) * t3 + qJD(5) * t27 - t138, -qJD(3) * t17 - t137, -t136, 0, t156; -qJD(2) * t15 - qJD(3) * t34 - qJD(4) * t27 - t161, -t157, -t155, -t156, 0;];
Cq = t7;
