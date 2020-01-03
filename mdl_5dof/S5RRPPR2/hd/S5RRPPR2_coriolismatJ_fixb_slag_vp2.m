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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:57:23
% EndTime: 2020-01-03 11:57:26
% DurationCPUTime: 0.96s
% Computational Cost: add. (2988->165), mult. (6286->236), div. (0->0), fcn. (5136->8), ass. (0->120)
t107 = sin(qJ(5));
t189 = t107 / 0.2e1;
t188 = qJD(1) + qJD(2);
t103 = sin(pkin(9));
t101 = t103 ^ 2;
t105 = cos(pkin(9));
t102 = t105 ^ 2;
t109 = cos(qJ(5));
t154 = t109 * mrSges(6,2);
t158 = t107 * mrSges(6,1);
t131 = t154 + t158;
t187 = t101 * t131 + (t101 + t102) * mrSges(5,3);
t150 = t103 * t107;
t81 = (-Ifges(6,5) * t107 - Ifges(6,6) * t109) * t103;
t186 = -t105 * t81 / 0.2e1 - (-Ifges(6,5) * t105 + (-0.2e1 * Ifges(6,4) * t107 + (Ifges(6,1) - Ifges(6,2)) * t109) * t103) * t150 / 0.2e1;
t104 = sin(pkin(8));
t97 = t104 * pkin(2) + qJ(4);
t90 = t101 * t97;
t132 = mrSges(6,1) * t109 - mrSges(6,2) * t107;
t185 = t103 * t132;
t106 = cos(pkin(8));
t108 = sin(qJ(2));
t146 = t106 * t108;
t110 = cos(qJ(2));
t168 = t110 * pkin(1);
t98 = pkin(2) + t168;
t135 = pkin(1) * t146 + t104 * t98;
t77 = qJ(4) + t135;
t66 = t101 * t77;
t164 = t66 + t90;
t148 = t105 * t107;
t125 = -t105 * pkin(4) - t103 * pkin(7) - pkin(3);
t96 = t104 * t108 * pkin(1);
t136 = t106 * t98 - t96;
t57 = t125 - t136;
t33 = t109 * t57 - t77 * t148;
t147 = t105 * t109;
t34 = t107 * t57 + t77 * t147;
t169 = t106 * pkin(2);
t84 = t125 - t169;
t51 = t109 * t84 - t97 * t148;
t52 = t107 * t84 + t97 * t147;
t184 = ((t34 + t52) * t109 + (-t33 - t51) * t107) * t105 + t164;
t175 = t109 ^ 2;
t176 = t107 ^ 2;
t143 = mrSges(6,3) * t150;
t87 = t105 * mrSges(6,2) - t143;
t149 = t103 * t109;
t142 = mrSges(6,3) * t149;
t88 = -t105 * mrSges(6,1) - t142;
t112 = (t175 / 0.2e1 + t176 / 0.2e1) * t101 * mrSges(6,3) + (t105 * t132 / 0.2e1 + t87 * t189 + t109 * t88 / 0.2e1) * t103;
t183 = qJD(3) * t112;
t119 = (-t154 / 0.2e1 - t158 / 0.2e1) * t105;
t180 = -t109 * t87 / 0.2e1 + t88 * t189;
t114 = t119 + t180;
t182 = qJD(4) * t114;
t116 = t119 - t180;
t181 = t116 * qJD(4);
t59 = (t175 - 0.1e1 + t176) * t105 * t103;
t170 = m(6) * t59;
t144 = t170 / 0.2e1;
t165 = qJD(4) * t144 - qJD(5) * t112;
t173 = -m(6) / 0.2e1;
t141 = qJD(3) * t173;
t179 = -qJD(5) * t114 + t59 * t141;
t178 = qJD(3) * t144 + qJD(5) * t116;
t85 = (t104 * t110 + t146) * pkin(1);
t86 = t106 * t168 - t96;
t43 = t109 * t85 - t86 * t148;
t44 = t107 * t85 + t86 * t147;
t177 = (-mrSges(4,2) + t187) * t86 + t43 * t88 + t44 * t87 + (-t108 * mrSges(3,1) - t110 * mrSges(3,2)) * pkin(1);
t174 = m(5) / 0.2e1;
t172 = m(6) / 0.2e1;
t23 = (-t105 * t86 - t107 * t43 + t109 * t44) * t103;
t171 = m(6) * t23;
t167 = t33 * t87;
t166 = t34 * t88;
t162 = Ifges(6,4) * t109;
t161 = t102 * t77;
t160 = t102 * t97;
t157 = t107 * t33;
t122 = t87 * t147 - t88 * t148 + t187;
t130 = -t109 * t34 + t157;
t10 = m(5) * (t66 + t161) + m(6) * (-t130 * t105 + t66) + t122;
t151 = t10 * qJD(1);
t145 = t171 / 0.2e1;
t138 = -t149 / 0.2e1;
t137 = t149 / 0.2e1;
t72 = -Ifges(6,6) * t105 + (-Ifges(6,2) * t107 + t162) * t103;
t83 = (-Ifges(6,1) * t107 - t162) * t103;
t134 = t83 * t137 + t72 * t138 + t186;
t13 = m(6) * (t90 + (-t107 * t51 + t109 * t52) * t105) + m(5) * (t90 + t160) + t122;
t113 = ((t77 + t97) * t102 + t164) * t174 + t122;
t115 = t85 * t174 + (t107 * t44 + t109 * t43) * t172;
t5 = t184 * t173 - t113 + t115;
t129 = t5 * qJD(1) - t13 * qJD(2);
t121 = t77 * t185;
t4 = t167 - t166 + (t130 * mrSges(6,3) + t121) * t103 + t134;
t128 = t4 * qJD(1) - t183;
t127 = t188 * t112;
t126 = t188 * t114;
t124 = t43 * mrSges(6,1) / 0.2e1 - t44 * mrSges(6,2) / 0.2e1;
t47 = t86 * t66;
t89 = -t105 * mrSges(5,1) + t103 * mrSges(5,2);
t3 = (-mrSges(4,1) + t89) * t85 + m(5) * (t86 * t161 + t47 + (-pkin(3) - t136) * t85) + m(6) * (t33 * t43 + t34 * t44 + t47) + m(4) * (t135 * t86 - t136 * t85) + t177;
t123 = -t3 * qJD(1) + t23 * t141;
t37 = t51 * t87;
t38 = t52 * t88;
t45 = t51 * t143;
t60 = t132 * t90;
t111 = (t121 / 0.2e1 + (t157 / 0.2e1 + (-t34 / 0.2e1 - t52 / 0.2e1) * t109) * mrSges(6,3)) * t103 + t37 / 0.2e1 - t38 / 0.2e1 + t45 / 0.2e1 + t60 / 0.2e1 + t167 / 0.2e1 - t166 / 0.2e1 + t134;
t1 = t111 - t124;
t7 = t72 * t137 + t83 * t138 + t52 * t142 - t186 - t37 + t38 - t45 - t60;
t120 = t1 * qJD(1) - t7 * qJD(2) - t183;
t61 = t86 * t90;
t32 = 0.2e1 * (qJD(1) / 0.4e1 + qJD(2) / 0.4e1) * t170;
t8 = qJD(2) * t145 + t165;
t6 = t184 * t172 + t113 + t115;
t2 = t111 + t124;
t9 = [qJD(2) * t3 + qJD(4) * t10 + qJD(5) * t4, t6 * qJD(4) + t2 * qJD(5) - t123 + (-t85 * mrSges(4,1) + t85 * t89 + 0.2e1 * (t86 * t160 + t61 + (-pkin(3) - t169) * t85) * t174 + 0.2e1 * (t51 * t43 + t52 * t44 + t61) * t172 + m(4) * (t104 * t86 - t106 * t85) * pkin(2) + t177) * qJD(2), t8, t6 * qJD(2) + t151 + t178, t2 * qJD(2) + t181 + (-t34 * mrSges(6,1) - t33 * mrSges(6,2) + t81) * qJD(5) + t128; -t5 * qJD(4) + t1 * qJD(5) + t123, qJD(4) * t13 - qJD(5) * t7, -qJD(1) * t171 / 0.2e1 + t165, -t129 + t178, t181 + (-t52 * mrSges(6,1) - t51 * mrSges(6,2) + t81) * qJD(5) + t120; t8, qJD(1) * t145 + t165, 0, t32, -qJD(5) * t185 - t127; t5 * qJD(2) - t151 + t179, t129 + t179, -t32, 0, -t131 * qJD(5) - t126; -qJD(2) * t1 - t128 + t182, -t120 + t182, t127, t126, 0;];
Cq = t9;
