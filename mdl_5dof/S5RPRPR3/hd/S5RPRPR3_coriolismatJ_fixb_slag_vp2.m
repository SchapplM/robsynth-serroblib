% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:38
% EndTime: 2022-01-23 09:20:40
% DurationCPUTime: 0.90s
% Computational Cost: add. (2862->154), mult. (5978->221), div. (0->0), fcn. (4842->8), ass. (0->111)
t104 = sin(qJ(5));
t179 = t104 / 0.2e1;
t178 = qJD(1) + qJD(3);
t103 = cos(pkin(9));
t100 = t103 ^ 2;
t106 = cos(qJ(5));
t147 = t106 * mrSges(6,2);
t150 = t104 * mrSges(6,1);
t127 = t147 + t150;
t101 = sin(pkin(9));
t99 = t101 ^ 2;
t177 = t99 * t127 + (t100 + t99) * mrSges(5,3);
t143 = t101 * t104;
t83 = (-Ifges(6,5) * t104 - Ifges(6,6) * t106) * t101;
t176 = -t103 * t83 / 0.2e1 - (-Ifges(6,5) * t103 + (-0.2e1 * Ifges(6,4) * t104 + (Ifges(6,1) - Ifges(6,2)) * t106) * t101) * t143 / 0.2e1;
t96 = t99 * qJ(4);
t128 = mrSges(6,1) * t106 - mrSges(6,2) * t104;
t175 = t101 * t128;
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t159 = pkin(1) * sin(pkin(8));
t95 = cos(pkin(8)) * pkin(1) + pkin(2);
t81 = t105 * t95 + t107 * t159;
t77 = qJ(4) + t81;
t65 = t99 * t77;
t155 = t65 + t96;
t141 = t103 * t104;
t80 = -t105 * t159 + t107 * t95;
t88 = -t103 * pkin(4) - t101 * pkin(7) - pkin(3);
t55 = -t80 + t88;
t33 = t106 * t55 - t77 * t141;
t140 = t103 * t106;
t34 = t104 * t55 + t77 * t140;
t61 = -qJ(4) * t141 + t106 * t88;
t62 = qJ(4) * t140 + t104 * t88;
t174 = ((t34 + t62) * t106 + (-t33 - t61) * t104) * t103 + t155;
t37 = t106 * t81 - t80 * t141;
t38 = t104 * t81 + t80 * t140;
t137 = mrSges(6,3) * t143;
t86 = t103 * mrSges(6,2) - t137;
t142 = t101 * t106;
t136 = mrSges(6,3) * t142;
t87 = -t103 * mrSges(6,1) - t136;
t173 = (-mrSges(4,2) + t177) * t80 + t37 * t87 + t38 * t86;
t165 = t106 ^ 2;
t166 = t104 ^ 2;
t109 = (t165 / 0.2e1 + t166 / 0.2e1) * t99 * mrSges(6,3) + (t103 * t128 / 0.2e1 + t86 * t179 + t106 * t87 / 0.2e1) * t101;
t172 = qJD(2) * t109;
t116 = (-t147 / 0.2e1 - t150 / 0.2e1) * t103;
t169 = -t106 * t86 / 0.2e1 + t87 * t179;
t111 = t116 + t169;
t171 = qJD(4) * t111;
t113 = t116 - t169;
t170 = t113 * qJD(4);
t56 = (t165 - 0.1e1 + t166) * t103 * t101;
t160 = m(6) * t56;
t138 = t160 / 0.2e1;
t156 = qJD(4) * t138 - qJD(5) * t109;
t163 = -m(6) / 0.2e1;
t135 = qJD(2) * t163;
t168 = -qJD(5) * t111 + t56 * t135;
t167 = qJD(2) * t138 + qJD(5) * t113;
t164 = m(5) / 0.2e1;
t162 = m(6) / 0.2e1;
t19 = (-t103 * t80 - t104 * t37 + t106 * t38) * t101;
t161 = m(6) * t19;
t158 = t33 * t86;
t157 = t34 * t87;
t153 = Ifges(6,4) * t106;
t152 = t100 * t77;
t149 = t104 * t33;
t119 = t86 * t140 - t87 * t141 + t177;
t126 = -t106 * t34 + t149;
t10 = m(6) * (-t126 * t103 + t65) + m(5) * (t65 + t152) + t119;
t145 = qJD(1) * t10;
t144 = t100 * qJ(4);
t139 = t161 / 0.2e1;
t132 = -t142 / 0.2e1;
t131 = t142 / 0.2e1;
t73 = -Ifges(6,6) * t103 + (-Ifges(6,2) * t104 + t153) * t101;
t85 = (-Ifges(6,1) * t104 - t153) * t101;
t130 = t85 * t131 + t73 * t132 + t176;
t118 = t77 * t175;
t4 = t158 - t157 + (t126 * mrSges(6,3) + t118) * t101 + t130;
t125 = t4 * qJD(1) - t172;
t15 = m(6) * (t96 + (-t104 * t61 + t106 * t62) * t103) + m(5) * (t96 + t144) + t119;
t110 = ((qJ(4) + t77) * t100 + t155) * t164 + t119;
t112 = (t104 * t38 + t106 * t37) * t162 + t81 * t164;
t5 = t174 * t163 - t110 + t112;
t124 = qJD(1) * t5 - qJD(3) * t15;
t123 = t178 * t109;
t122 = t178 * t111;
t121 = t37 * mrSges(6,1) / 0.2e1 - t38 * mrSges(6,2) / 0.2e1;
t39 = t80 * t65;
t89 = -mrSges(5,1) * t103 + mrSges(5,2) * t101;
t3 = (-mrSges(4,1) + t89) * t81 + m(6) * (t33 * t37 + t34 * t38 + t39) + m(5) * (t80 * t152 + t39 + (-pkin(3) - t80) * t81) + t173;
t120 = -t3 * qJD(1) + t19 * t135;
t43 = t61 * t86;
t44 = t62 * t87;
t50 = t61 * t137;
t72 = t128 * t96;
t108 = (t118 / 0.2e1 + (t149 / 0.2e1 + (-t34 / 0.2e1 - t62 / 0.2e1) * t106) * mrSges(6,3)) * t101 + t43 / 0.2e1 - t44 / 0.2e1 + t50 / 0.2e1 + t72 / 0.2e1 + t158 / 0.2e1 - t157 / 0.2e1 + t130;
t1 = t108 - t121;
t8 = t73 * t131 + t85 * t132 + t62 * t136 - t176 - t43 + t44 - t50 - t72;
t117 = t1 * qJD(1) - t8 * qJD(3) - t172;
t64 = t80 * t96;
t29 = 0.2e1 * (qJD(1) / 0.4e1 + qJD(3) / 0.4e1) * t160;
t7 = qJD(3) * t139 + t156;
t6 = t174 * t162 + t110 + t112;
t2 = t108 + t121;
t9 = [qJD(3) * t3 + qJD(4) * t10 + qJD(5) * t4, t7, t6 * qJD(4) + t2 * qJD(5) - t120 + (-t81 * mrSges(4,1) + t81 * t89 + 0.2e1 * (t61 * t37 + t62 * t38 + t64) * t162 + 0.2e1 * (-pkin(3) * t81 + t80 * t144 + t64) * t164 + t173) * qJD(3), qJD(3) * t6 + t145 + t167, t2 * qJD(3) + t170 + (-mrSges(6,1) * t34 - mrSges(6,2) * t33 + t83) * qJD(5) + t125; t7, 0, qJD(1) * t139 + t156, t29, -qJD(5) * t175 - t123; -t5 * qJD(4) + t1 * qJD(5) + t120, -qJD(1) * t161 / 0.2e1 + t156, qJD(4) * t15 - qJD(5) * t8, -t124 + t167, t170 + (-mrSges(6,1) * t62 - mrSges(6,2) * t61 + t83) * qJD(5) + t117; qJD(3) * t5 - t145 + t168, -t29, t124 + t168, 0, -t127 * qJD(5) - t122; -qJD(3) * t1 - t125 + t171, t123, -t117 + t171, t122, 0;];
Cq = t9;
