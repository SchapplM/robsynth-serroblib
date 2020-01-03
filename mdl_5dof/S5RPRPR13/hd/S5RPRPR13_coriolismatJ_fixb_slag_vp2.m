% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:59
% EndTime: 2019-12-31 18:32:02
% DurationCPUTime: 1.29s
% Computational Cost: add. (4092->237), mult. (8150->325), div. (0->0), fcn. (8556->6), ass. (0->133)
t202 = -mrSges(4,1) + mrSges(5,2);
t125 = sin(qJ(5));
t121 = t125 ^ 2;
t127 = cos(qJ(5));
t122 = t127 ^ 2;
t158 = t121 + t122;
t192 = -t127 / 0.2e1;
t194 = -t125 / 0.2e1;
t201 = -Ifges(6,5) * t194 - Ifges(6,6) * t192 - Ifges(4,4) - Ifges(5,6);
t123 = sin(pkin(8));
t124 = cos(pkin(8));
t126 = sin(qJ(3));
t190 = cos(qJ(3));
t109 = t190 * t123 + t126 * t124;
t200 = 0.2e1 * t109;
t199 = m(6) / 0.2e1;
t198 = pkin(3) + pkin(7);
t107 = t123 * t126 - t190 * t124;
t167 = qJ(4) * t107;
t72 = pkin(3) * t109 + t167;
t197 = m(5) * t72;
t196 = -t107 / 0.2e1;
t195 = t109 / 0.2e1;
t193 = t125 / 0.2e1;
t191 = t127 / 0.2e1;
t189 = mrSges(4,3) + mrSges(5,1);
t188 = pkin(6) + qJ(2);
t187 = mrSges(6,3) * t107;
t186 = Ifges(6,4) * t125;
t185 = Ifges(6,4) * t127;
t184 = Ifges(6,5) * t109;
t183 = Ifges(6,6) * t109;
t165 = t107 * t125;
t68 = t109 * mrSges(6,1) - mrSges(6,3) * t165;
t171 = t127 * t68;
t164 = t107 * t127;
t70 = -t109 * mrSges(6,2) + mrSges(6,3) * t164;
t175 = t125 * t70;
t156 = -pkin(2) * t124 - pkin(1);
t142 = -qJ(4) * t109 + t156;
t35 = t198 * t107 + t142;
t110 = t188 * t124;
t152 = t188 * t123;
t74 = t110 * t126 + t190 * t152;
t42 = pkin(4) * t109 + t74;
t24 = -t125 * t35 + t127 * t42;
t25 = t125 * t42 + t127 * t35;
t75 = t190 * t110 - t126 * t152;
t44 = -t107 * pkin(4) + t75;
t174 = t127 * mrSges(6,1);
t178 = t125 * mrSges(6,2);
t111 = -t174 + t178;
t61 = t111 * t107;
t5 = (t189 * t107 - t61) * t107 + (t189 * t109 + t171 + t175) * t109 + m(6) * (-t107 * t44 + (t125 * t25 + t127 * t24) * t109) + (m(3) * qJ(2) + mrSges(3,3)) * (t123 ^ 2 + t124 ^ 2) + (m(5) + m(4)) * (-t107 * t75 + t109 * t74);
t182 = qJD(1) * t5;
t170 = t127 * t70;
t176 = t125 * t68;
t60 = pkin(3) * t107 + t142;
t73 = -mrSges(5,2) * t107 - mrSges(5,3) * t109;
t8 = (-t170 + t176 + m(6) * (t125 * t24 - t127 * t25) - t73 - m(5) * t60) * t109;
t181 = qJD(1) * t8;
t102 = t107 * mrSges(5,3);
t103 = t107 * mrSges(4,2);
t41 = t198 * t109 + t167;
t28 = -t125 * t41 + t127 * t44;
t29 = t125 * t44 + t127 * t41;
t148 = Ifges(6,2) * t127 + t186;
t36 = t148 * t107 + t183;
t37 = -Ifges(6,6) * t107 + t148 * t109;
t99 = Ifges(6,4) * t164;
t38 = Ifges(6,1) * t165 + t184 + t99;
t149 = Ifges(6,1) * t125 + t185;
t39 = -Ifges(6,5) * t107 + t149 * t109;
t62 = t111 * t109;
t163 = t109 * t125;
t69 = -t107 * mrSges(6,1) - mrSges(6,3) * t163;
t162 = t109 * t127;
t71 = t107 * mrSges(6,2) + mrSges(6,3) * t162;
t1 = -t42 * t61 + t44 * t62 + t28 * t68 + t24 * t69 + t29 * t70 + t25 * t71 + t72 * t73 - t156 * t103 + m(6) * (t24 * t28 + t25 * t29 - t42 * t44) + (t156 * mrSges(4,1) + t201 * t109 + t36 * t191 + t38 * t193) * t109 + (-t109 * mrSges(5,2) + t102 + t197) * t60 + (t37 * t191 + t39 * t193 + (Ifges(5,3) - Ifges(4,1) - Ifges(5,2) + Ifges(4,2) - Ifges(6,3)) * t109 - t201 * t107) * t107;
t180 = t1 * qJD(1);
t179 = t125 * mrSges(6,1);
t177 = t125 * t29;
t173 = t127 * mrSges(6,2);
t172 = t127 * t28;
t112 = t173 + t179;
t63 = t107 * t112;
t64 = -Ifges(6,2) * t165 + t99;
t116 = Ifges(6,1) * t127 - t186;
t65 = t107 * t116;
t98 = Ifges(6,5) * t164;
t2 = t44 * t63 + t98 * t195 + t24 * t70 - t25 * t68 + ((t38 / 0.2e1 + t64 / 0.2e1 - t24 * mrSges(6,3)) * t127 + (t65 / 0.2e1 - t36 / 0.2e1 - t183 / 0.2e1 - t25 * mrSges(6,3)) * t125) * t107;
t169 = t2 * qJD(1);
t151 = t158 * t109;
t129 = -t197 / 0.2e1 + (-t151 * t198 - t167) * t199 + t112 * t196;
t130 = t197 / 0.2e1 + (-t125 * t28 + t127 * t29) * t199 + t69 * t194 + t71 * t191;
t153 = t121 / 0.2e1 + t122 / 0.2e1;
t150 = t153 * mrSges(6,3);
t6 = t103 - t102 + (-t150 + t202) * t109 + t129 - t130;
t168 = t6 * qJD(1);
t132 = (t179 / 0.2e1 + t173 / 0.2e1) * t109;
t135 = t176 / 0.2e1 - t170 / 0.2e1;
t10 = t153 * t187 + t132 + t135;
t166 = t10 * qJD(1);
t136 = t178 / 0.2e1 - t174 / 0.2e1;
t131 = t136 * t109;
t134 = -t175 / 0.2e1 - t171 / 0.2e1;
t13 = -t131 - t134;
t161 = t13 * qJD(1);
t21 = (-m(5) / 0.2e1 - t153 * m(6)) * t200;
t160 = t21 * qJD(1);
t114 = -Ifges(6,2) * t125 + t185;
t155 = -t149 / 0.4e1 - t114 / 0.4e1;
t154 = t116 / 0.4e1 - t148 / 0.4e1;
t146 = t172 + t177;
t30 = -qJ(4) * t111 + (-t149 / 0.2e1 - t114 / 0.2e1) * t127 + (-t116 / 0.2e1 + t148 / 0.2e1) * t125;
t137 = qJ(4) * t63 / 0.2e1 - t44 * t111 / 0.2e1;
t139 = -t38 / 0.4e1 - t64 / 0.4e1 + t198 * t68 / 0.2e1;
t140 = -t36 / 0.4e1 + t65 / 0.4e1 - t198 * t70 / 0.2e1;
t141 = -t28 * mrSges(6,1) / 0.2e1 + t29 * mrSges(6,2) / 0.2e1;
t143 = t198 * t150;
t4 = (Ifges(6,3) / 0.2e1 + t143) * t107 + (-0.3e1 / 0.4e1 * t183 + t154 * t107 + t140) * t127 + (-0.3e1 / 0.4e1 * t184 + t155 * t107 + t139) * t125 + t137 + t141;
t145 = -t4 * qJD(1) - t30 * qJD(3);
t101 = mrSges(5,3) + (m(5) + m(6)) * qJ(4) + t112;
t133 = t69 * t192 + t71 * t194;
t12 = t136 * t107 + 0.2e1 * (t44 / 0.4e1 - t177 / 0.4e1 - t172 / 0.4e1) * m(6) + t133;
t144 = qJD(1) * t12 + qJD(3) * t101;
t20 = m(5) * t195 + t151 * t199 + (-m(6) * t158 / 0.4e1 - m(5) / 0.4e1) * t200;
t14 = -t131 + t134;
t11 = t132 - t135 - t158 * t187 / 0.2e1;
t9 = m(5) * t75 + (-mrSges(5,1) + t136) * t107 - t133 + (t146 + t44) * t199;
t7 = -t109 * t150 + t129 + t130;
t3 = Ifges(6,5) * t163 / 0.2e1 + Ifges(6,6) * t162 / 0.2e1 + Ifges(6,3) * t196 + (-t183 / 0.4e1 + t140) * t127 + (-t184 / 0.4e1 + t139) * t125 + (t155 * t125 + t154 * t127 + t143) * t107 + t137 - t141;
t15 = [qJD(2) * t5 + qJD(3) * t1 + qJD(4) * t8 + qJD(5) * t2, qJD(3) * t7 + qJD(4) * t20 + qJD(5) * t14 + t182, t7 * qJD(2) + t9 * qJD(4) + t3 * qJD(5) + t180 + (qJ(4) * t62 - t42 * t112 + (-t198 * t69 - t28 * mrSges(6,3) + t39 / 0.2e1) * t127 + (-t198 * t71 - t29 * mrSges(6,3) - t37 / 0.2e1) * t125 + 0.2e1 * (-qJ(4) * t42 - t146 * t198) * t199 + (pkin(3) * mrSges(5,1) + Ifges(6,5) * t192 + Ifges(6,6) * t193 + Ifges(5,4) - Ifges(4,5)) * t107 + (-qJ(4) * mrSges(5,1) + t114 * t191 + t116 * t193 + Ifges(5,5) - Ifges(4,6)) * t109 + (-m(5) * pkin(3) + t202) * t75 + (-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t74) * qJD(3), qJD(2) * t20 + qJD(3) * t9 + qJD(5) * t11 + t181, t169 + t14 * qJD(2) + t3 * qJD(3) + t11 * qJD(4) + (-mrSges(6,1) * t25 - mrSges(6,2) * t24 - Ifges(6,6) * t165 + t98) * qJD(5); -qJD(3) * t6 + qJD(4) * t21 - qJD(5) * t13 - t182, 0, -t168, t160, qJD(5) * t111 - t161; qJD(2) * t6 + qJD(4) * t12 + qJD(5) * t4 - t180, t168, qJD(4) * t101 + qJD(5) * t30, t144, ((mrSges(6,2) * t198 - Ifges(6,6)) * t127 + (mrSges(6,1) * t198 - Ifges(6,5)) * t125) * qJD(5) - t145; -qJD(2) * t21 - qJD(3) * t12 - qJD(5) * t10 - t181, -t160, -t144, 0, -qJD(5) * t112 - t166; qJD(2) * t13 - qJD(3) * t4 + qJD(4) * t10 - t169, t161, t145, t166, 0;];
Cq = t15;
