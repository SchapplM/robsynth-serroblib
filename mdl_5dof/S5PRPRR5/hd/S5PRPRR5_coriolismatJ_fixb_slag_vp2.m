% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:31
% EndTime: 2019-12-05 15:53:36
% DurationCPUTime: 1.38s
% Computational Cost: add. (4141->161), mult. (9460->240), div. (0->0), fcn. (10592->8), ass. (0->93)
t117 = sin(pkin(9));
t118 = cos(pkin(9));
t173 = sin(qJ(4));
t174 = cos(qJ(4));
t102 = -t173 * t117 + t174 * t118;
t103 = -t174 * t117 - t173 * t118;
t119 = sin(qJ(5));
t121 = cos(qJ(5));
t133 = t121 * t102 + t103 * t119;
t153 = pkin(6) + qJ(3);
t106 = t153 * t117;
t107 = t153 * t118;
t86 = -t174 * t106 - t173 * t107;
t127 = t103 * pkin(7) + t86;
t87 = -t173 * t106 + t174 * t107;
t65 = t102 * pkin(7) + t87;
t187 = -t119 * t65 + t121 * t127;
t35 = t119 * t127 + t121 * t65;
t84 = t102 * t119 - t103 * t121;
t7 = -t35 * mrSges(6,1) - t187 * mrSges(6,2) + Ifges(6,5) * t133 - Ifges(6,6) * t84;
t207 = t7 * qJD(5);
t172 = Ifges(6,4) * t84;
t177 = t84 / 0.2e1;
t193 = t133 / 0.2e1;
t200 = -t84 / 0.2e1;
t197 = t84 * mrSges(6,1);
t203 = t133 * mrSges(6,2) + t197;
t111 = -pkin(3) * t118 - pkin(2);
t91 = -pkin(4) * t102 + t111;
t4 = t91 * t203 + (Ifges(6,2) * t133 + t172) * t200 + (Ifges(6,1) * t133 - t172) * t177 + (0.2e1 * Ifges(6,4) * t133 + (Ifges(6,1) - Ifges(6,2)) * t84) * t193;
t139 = t197 / 0.2e1;
t120 = sin(qJ(2));
t92 = t103 * t120;
t94 = t102 * t120;
t135 = -t119 * t94 + t121 * t92;
t61 = t119 * t92 + t121 * t94;
t14 = -t61 * mrSges(6,1) - t135 * mrSges(6,2);
t205 = t14 * qJD(5);
t183 = m(6) * pkin(4);
t141 = -t183 / 0.2e1;
t122 = cos(qJ(2));
t93 = t103 * t122;
t95 = t102 * t122;
t59 = -t119 * t95 + t121 * t93;
t62 = t119 * t93 + t121 * t95;
t152 = t59 * mrSges(6,1) / 0.2e1 - t62 * mrSges(6,2) / 0.2e1;
t204 = -t93 * mrSges(5,1) / 0.2e1 + t95 * mrSges(5,2) / 0.2e1 + (t119 * t62 + t121 * t59) * t141 - t152;
t195 = t193 - t133 / 0.2e1;
t85 = -t103 * mrSges(5,1) + t102 * mrSges(5,2);
t192 = t203 + t85;
t186 = t103 ^ 2;
t185 = m(5) / 0.2e1;
t184 = m(6) / 0.2e1;
t176 = t120 / 0.2e1;
t175 = -t122 / 0.2e1;
t167 = t103 * pkin(4);
t150 = t119 * t133;
t148 = t121 * t84;
t115 = t117 ^ 2;
t116 = t118 ^ 2;
t142 = t115 + t116;
t146 = t120 * t122;
t12 = m(6) * (t135 * t59 + t61 * t62 - t146) + m(5) * (t92 * t93 + t94 * t95 - t146) + m(4) * (-0.1e1 + t142) * t146;
t147 = t12 * qJD(1);
t18 = (-t103 / 0.2e1 - t150 / 0.2e1 + t148 / 0.2e1) * t183 + t192;
t144 = t18 * qJD(2);
t23 = 0.2e1 * t193 * mrSges(6,2) + 0.2e1 * t139;
t143 = t23 * qJD(2);
t140 = m(4) * t176;
t138 = t203 * t175;
t134 = t142 * mrSges(4,3);
t132 = t142 * qJ(3);
t47 = -mrSges(6,1) * t133 + mrSges(6,2) * t84;
t1 = t111 * t85 - Ifges(5,4) * t186 + (Ifges(5,4) * t102 + (-Ifges(5,1) + Ifges(5,2)) * t103) * t102 + (-m(6) * t91 - t47) * t167 + t4;
t123 = t195 * mrSges(6,3) * t135 + t122 * t167 * t184;
t2 = (-t203 / 0.2e1 - t85 / 0.2e1) * t122 + t123 + t204;
t131 = t2 * qJD(1) + t1 * qJD(2);
t5 = t138 - t152;
t130 = t5 * qJD(1) + t4 * qJD(2);
t11 = (t133 ^ 2 + t84 ^ 2) * mrSges(6,3) + (t102 ^ 2 + t186) * mrSges(5,3) + t134 + m(6) * (t133 * t35 - t187 * t84) + m(5) * (t102 * t87 + t103 * t86) + m(4) * t132;
t124 = (t133 * t61 - t135 * t84) * t184 + (t102 * t94 + t103 * t92) * t185;
t16 = (-m(6) / 0.2e1 - m(5) / 0.2e1 + (t115 / 0.2e1 + t116 / 0.2e1 - 0.1e1 / 0.2e1) * m(4)) * t120 + t124;
t129 = -qJD(1) * t16 - qJD(2) * t11;
t105 = (mrSges(6,1) * t119 + mrSges(6,2) * t121) * pkin(4);
t8 = (t200 + t177) * Ifges(6,6) + t195 * Ifges(6,5);
t126 = -qJD(2) * t8 + qJD(4) * t105;
t104 = t105 * qJD(5);
t27 = t103 * t141 + (-t148 + t150) * t183 / 0.2e1;
t24 = t139 - t197 / 0.2e1;
t15 = t142 * t140 + t124 + t140 + (m(5) + m(6)) * t176;
t6 = t138 + t152;
t3 = t192 * t175 + t123 - t204;
t9 = [t12 * qJD(2), t15 * qJD(3) + t3 * qJD(4) + t6 * qJD(5) + t147 + ((-mrSges(3,2) + t134) * t122 + (-t118 * mrSges(4,1) - t102 * mrSges(5,1) + t117 * mrSges(4,2) - t103 * mrSges(5,2) - mrSges(3,1) + t47) * t120 + 0.2e1 * (t120 * t91 + t187 * t59 + t35 * t62) * t184 + 0.2e1 * (t111 * t120 + t86 * t93 + t87 * t95) * t185 + m(4) * (-t120 * pkin(2) + t122 * t132) + (t133 * t62 - t59 * t84) * mrSges(6,3) + (t102 * t95 + t103 * t93) * mrSges(5,3)) * qJD(2), t15 * qJD(2), t3 * qJD(2) + (-t94 * mrSges(5,1) - t92 * mrSges(5,2) + (t119 * t135 - t121 * t61) * t183 + t14) * qJD(4) + t205, t6 * qJD(2) + t14 * qJD(4) + t205; qJD(3) * t16 + qJD(4) * t2 + qJD(5) * t5 - t147, qJD(3) * t11 + qJD(4) * t1 + qJD(5) * t4, qJD(4) * t27 + qJD(5) * t24 - t129, t27 * qJD(3) + t207 + t131 + (-t87 * mrSges(5,1) - t86 * mrSges(5,2) + Ifges(5,5) * t102 + Ifges(5,6) * t103 + (m(6) * (t119 * t187 - t121 * t35) + (-t119 * t84 - t121 * t133) * mrSges(6,3)) * pkin(4) + t7) * qJD(4), t24 * qJD(3) + t7 * qJD(4) + t130 + t207; -t16 * qJD(2), qJD(4) * t18 + qJD(5) * t23 + t129, 0, t144, t143; -qJD(2) * t2, -qJD(3) * t18 + qJD(5) * t8 - t131, -t144, -t104, -t104 - t126; -t5 * qJD(2), -qJD(3) * t23 - qJD(4) * t8 - t130, -t143, t126, 0;];
Cq = t9;
