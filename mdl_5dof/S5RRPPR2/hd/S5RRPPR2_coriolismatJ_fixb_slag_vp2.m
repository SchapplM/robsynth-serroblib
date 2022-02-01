% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:30
% EndTime: 2022-01-20 10:05:33
% DurationCPUTime: 1.00s
% Computational Cost: add. (2988->162), mult. (6286->233), div. (0->0), fcn. (5136->8), ass. (0->121)
t108 = sin(qJ(5));
t188 = t108 / 0.2e1;
t187 = qJD(1) + qJD(2);
t104 = sin(pkin(9));
t106 = cos(pkin(9));
t110 = cos(qJ(5));
t150 = t104 * t108;
t83 = (-Ifges(6,5) * t108 - Ifges(6,6) * t110) * t104;
t186 = -t106 * t83 / 0.2e1 - (-Ifges(6,5) * t106 + (-0.2e1 * Ifges(6,4) * t108 + (Ifges(6,1) - Ifges(6,2)) * t110) * t104) * t150 / 0.2e1;
t102 = t104 ^ 2;
t105 = sin(pkin(8));
t100 = pkin(2) * t105 + qJ(4);
t92 = t102 * t100;
t132 = mrSges(6,1) * t110 - mrSges(6,2) * t108;
t185 = t104 * t132;
t176 = t110 ^ 2;
t177 = t108 ^ 2;
t143 = mrSges(6,3) * t150;
t89 = t106 * mrSges(6,2) - t143;
t149 = t104 * t110;
t142 = mrSges(6,3) * t149;
t90 = -t106 * mrSges(6,1) - t142;
t114 = (t176 / 0.2e1 + t177 / 0.2e1) * t102 * mrSges(6,3) + (t106 * t132 / 0.2e1 + t89 * t188 + t110 * t90 / 0.2e1) * t104;
t184 = qJD(3) * t114;
t155 = t110 * mrSges(6,2);
t158 = t108 * mrSges(6,1);
t120 = (-t155 / 0.2e1 - t158 / 0.2e1) * t106;
t181 = -t110 * t89 / 0.2e1 + t90 * t188;
t115 = t120 + t181;
t183 = qJD(4) * t115;
t118 = t120 - t181;
t182 = t118 * qJD(4);
t61 = (t176 - 0.1e1 + t177) * t106 * t104;
t171 = m(6) * t61;
t144 = t171 / 0.2e1;
t164 = qJD(4) * t144 - qJD(5) * t114;
t174 = -m(6) / 0.2e1;
t141 = qJD(3) * t174;
t180 = -qJD(5) * t115 + t61 * t141;
t179 = qJD(3) * t144 + qJD(5) * t118;
t109 = sin(qJ(2));
t111 = cos(qJ(2));
t103 = t106 ^ 2;
t131 = t155 + t158;
t123 = t102 * t131 + (t102 + t103) * mrSges(5,3);
t148 = t106 * t108;
t107 = cos(pkin(8));
t146 = t107 * t109;
t87 = (t105 * t111 + t146) * pkin(1);
t170 = pkin(1) * t111;
t99 = t105 * t109 * pkin(1);
t88 = t107 * t170 - t99;
t44 = t110 * t87 - t88 * t148;
t147 = t106 * t110;
t45 = t108 * t87 + t88 * t147;
t178 = (-mrSges(4,2) + t123) * t88 + t44 * t90 + t45 * t89 + (-t109 * mrSges(3,1) - t111 * mrSges(3,2)) * pkin(1);
t175 = m(5) / 0.2e1;
t173 = m(6) / 0.2e1;
t23 = (-t106 * t88 - t108 * t44 + t110 * t45) * t104;
t172 = m(6) * t23;
t169 = pkin(2) * t107;
t126 = -pkin(4) * t106 - pkin(7) * t104 - pkin(3);
t101 = pkin(2) + t170;
t136 = t101 * t107 - t99;
t59 = t126 - t136;
t135 = pkin(1) * t146 + t105 * t101;
t79 = qJ(4) + t135;
t34 = t110 * t59 - t79 * t148;
t168 = t34 * t89;
t35 = t108 * t59 + t79 * t147;
t167 = t35 * t90;
t68 = t102 * t79;
t166 = t35 * t147 + t68;
t86 = t126 - t169;
t54 = t100 * t147 + t108 * t86;
t165 = t54 * t147 + t92;
t69 = t103 * t79;
t163 = t69 + t68;
t93 = t103 * t100;
t162 = t93 + t92;
t160 = Ifges(6,4) * t110;
t152 = t34 * t108;
t116 = t89 * t147 - t90 * t148 + t123;
t10 = m(6) * (-t34 * t148 + t166) + m(5) * t163 + t116;
t151 = qJD(1) * t10;
t145 = t172 / 0.2e1;
t138 = -t149 / 0.2e1;
t137 = t149 / 0.2e1;
t75 = -Ifges(6,6) * t106 + (-Ifges(6,2) * t108 + t160) * t104;
t85 = (-Ifges(6,1) * t108 - t160) * t104;
t134 = t85 * t137 + t75 * t138 + t186;
t53 = -t100 * t148 + t110 * t86;
t13 = m(6) * (-t53 * t148 + t165) + m(5) * t162 + t116;
t113 = (t162 + t163) * t175 + ((-t34 - t53) * t148 + t165 + t166) * t173 + t116;
t117 = -m(5) * t87 / 0.2e1 + (t108 * t45 + t110 * t44) * t174;
t6 = t113 + t117;
t130 = -qJD(1) * t6 - qJD(2) * t13;
t122 = t79 * t185;
t4 = t168 - t167 + (t122 + (-t35 * t110 + t152) * mrSges(6,3)) * t104 + t134;
t129 = t4 * qJD(1) - t184;
t128 = t187 * t114;
t127 = t187 * t115;
t125 = t44 * mrSges(6,1) / 0.2e1 - t45 * mrSges(6,2) / 0.2e1;
t49 = t88 * t68;
t91 = -mrSges(5,1) * t106 + mrSges(5,2) * t104;
t3 = (-mrSges(4,1) + t91) * t87 + m(6) * (t34 * t44 + t35 * t45 + t49) + m(5) * (t88 * t69 + t49 + (-pkin(3) - t136) * t87) + m(4) * (t135 * t88 - t136 * t87) + t178;
t124 = -t3 * qJD(1) + t23 * t141;
t38 = t53 * t89;
t39 = t54 * t90;
t46 = t53 * t143;
t62 = t132 * t92;
t112 = (t122 / 0.2e1 + (t152 / 0.2e1 + (-t35 / 0.2e1 - t54 / 0.2e1) * t110) * mrSges(6,3)) * t104 + t38 / 0.2e1 - t39 / 0.2e1 + t46 / 0.2e1 + t62 / 0.2e1 + t168 / 0.2e1 - t167 / 0.2e1 + t134;
t1 = t112 - t125;
t7 = t75 * t137 + t85 * t138 + t54 * t142 - t186 - t38 + t39 - t46 - t62;
t121 = t1 * qJD(1) - t7 * qJD(2) - t184;
t63 = t88 * t92;
t32 = 0.2e1 * (qJD(1) / 0.4e1 + qJD(2) / 0.4e1) * t171;
t8 = qJD(2) * t145 + t164;
t5 = t113 - t117;
t2 = t112 + t125;
t9 = [qJD(2) * t3 + qJD(4) * t10 + qJD(5) * t4, t5 * qJD(4) + t2 * qJD(5) - t124 + (-t87 * mrSges(4,1) + t87 * t91 + 0.2e1 * (t53 * t44 + t54 * t45 + t63) * t173 + 0.2e1 * (t88 * t93 + t63 + (-pkin(3) - t169) * t87) * t175 + m(4) * (t105 * t88 - t107 * t87) * pkin(2) + t178) * qJD(2), t8, qJD(2) * t5 + t151 + t179, t2 * qJD(2) + t182 + (-mrSges(6,1) * t35 - mrSges(6,2) * t34 + t83) * qJD(5) + t129; t6 * qJD(4) + t1 * qJD(5) + t124, qJD(4) * t13 - qJD(5) * t7, -qJD(1) * t172 / 0.2e1 + t164, -t130 + t179, t182 + (-mrSges(6,1) * t54 - mrSges(6,2) * t53 + t83) * qJD(5) + t121; t8, qJD(1) * t145 + t164, 0, t32, -qJD(5) * t185 - t128; -qJD(2) * t6 - t151 + t180, t130 + t180, -t32, 0, -t131 * qJD(5) - t127; -qJD(2) * t1 - t129 + t183, -t121 + t183, t128, t127, 0;];
Cq = t9;
