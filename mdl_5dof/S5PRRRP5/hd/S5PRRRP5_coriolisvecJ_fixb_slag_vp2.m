% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:46
% EndTime: 2019-12-05 16:47:55
% DurationCPUTime: 3.04s
% Computational Cost: add. (1875->268), mult. (4899->375), div. (0->0), fcn. (3131->6), ass. (0->136)
t163 = Ifges(5,4) + Ifges(6,4);
t186 = Ifges(5,1) + Ifges(6,1);
t192 = Ifges(5,5) + Ifges(6,5);
t185 = Ifges(5,2) + Ifges(6,2);
t191 = Ifges(5,6) + Ifges(6,6);
t122 = sin(qJ(4));
t123 = sin(qJ(3));
t125 = cos(qJ(4));
t126 = cos(qJ(3));
t131 = t122 * t123 - t125 * t126;
t90 = t131 * qJD(2);
t194 = t163 * t90;
t100 = t122 * t126 + t123 * t125;
t91 = t100 * qJD(2);
t193 = t163 * t91;
t175 = -pkin(7) - pkin(6);
t141 = qJD(3) * t175;
t105 = t123 * t141;
t106 = t126 * t141;
t109 = t175 * t123;
t110 = t175 * t126;
t145 = qJD(4) * t125;
t146 = qJD(4) * t122;
t127 = cos(qJ(2));
t152 = qJD(1) * t127;
t184 = t125 * t105 + t122 * t106 + t109 * t145 + t110 * t146 + t131 * t152;
t130 = t100 * t127;
t65 = t122 * t109 - t125 * t110;
t183 = qJD(1) * t130 - t65 * qJD(4) - t105 * t122 + t125 * t106;
t119 = qJD(3) + qJD(4);
t190 = t191 * t119 - t185 * t90 + t193;
t189 = t192 * t119 + t186 * t91 - t194;
t60 = t119 * t100;
t188 = -qJ(5) * t60 - qJD(5) * t131 + t184;
t59 = t119 * t131;
t187 = qJ(5) * t59 - qJD(5) * t100 + t183;
t124 = sin(qJ(2));
t153 = qJD(1) * t124;
t112 = qJD(2) * pkin(6) + t153;
t138 = pkin(7) * qJD(2) + t112;
t82 = t138 * t126;
t74 = t122 * t82;
t81 = t138 * t123;
t77 = qJD(3) * pkin(3) - t81;
t31 = t125 * t77 - t74;
t85 = t91 * qJ(5);
t13 = t31 - t85;
t53 = mrSges(6,1) * t90 + mrSges(6,2) * t91;
t54 = mrSges(5,1) * t90 + mrSges(5,2) * t91;
t182 = -t53 - t54;
t84 = t131 * t124;
t120 = t123 ^ 2;
t121 = t126 ^ 2;
t179 = t90 / 0.2e1;
t176 = t91 / 0.2e1;
t174 = pkin(4) * t91;
t171 = -t123 / 0.2e1;
t170 = t123 / 0.2e1;
t168 = mrSges(5,3) * t90;
t167 = mrSges(6,3) * t90;
t166 = t91 * mrSges(5,3);
t68 = -mrSges(6,2) * t119 - t167;
t69 = -mrSges(5,2) * t119 - t168;
t162 = t68 + t69;
t70 = mrSges(6,1) * t119 - t91 * mrSges(6,3);
t71 = mrSges(5,1) * t119 - t166;
t161 = t70 + t71;
t34 = -t125 * t81 - t74;
t160 = Ifges(4,4) * t123;
t159 = qJ(5) * t90;
t50 = t60 * qJD(2);
t158 = t122 * t50;
t76 = t125 * t82;
t157 = t126 * Ifges(4,2);
t156 = Ifges(4,5) * qJD(3);
t155 = Ifges(4,6) * qJD(3);
t113 = -qJD(2) * pkin(2) - t152;
t154 = t113 * t124;
t144 = qJD(1) * qJD(2);
t115 = t124 * t144;
t148 = qJD(3) * t123;
t142 = pkin(3) * t148;
t97 = qJD(2) * t142 + t115;
t151 = qJD(2) * t123;
t150 = qJD(2) * t124;
t149 = qJD(2) * t126;
t147 = qJD(3) * t126;
t143 = qJD(2) * qJD(3);
t117 = -pkin(3) * t126 - pkin(2);
t140 = t112 * t147;
t49 = t59 * qJD(2);
t10 = t50 * mrSges(6,1) - t49 * mrSges(6,2);
t139 = t127 * t144;
t33 = t122 * t81 - t76;
t137 = (t120 + t121) * t112;
t64 = t125 * t109 + t110 * t122;
t135 = mrSges(4,1) * t123 + mrSges(4,2) * t126;
t32 = t122 * t77 + t76;
t111 = t126 * t139;
t72 = -t112 * t148 + t111;
t73 = -t123 * t139 - t140;
t134 = -t123 * t73 + t126 * t72;
t107 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t151;
t108 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t149;
t133 = t126 * t107 + t123 * t108;
t132 = t123 * t107 - t126 * t108;
t61 = -qJD(3) * t81 + t111;
t62 = -t140 + (-pkin(7) * t147 - t123 * t152) * qJD(2);
t5 = t122 * t62 + t125 * t61 + t77 * t145 - t146 * t82;
t96 = qJD(2) * t117 - t152;
t6 = -qJD(4) * t32 - t122 * t61 + t125 * t62;
t12 = pkin(4) * t119 + t13;
t2 = -qJ(5) * t50 - qJD(5) * t90 + t5;
t3 = qJ(5) * t49 - qJD(5) * t91 + t6;
t58 = pkin(4) * t90 + qJD(5) + t96;
t129 = t6 * mrSges(5,1) + t3 * mrSges(6,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2) - t12 * t167 - t31 * t168 - t58 * (mrSges(6,1) * t91 - mrSges(6,2) * t90) - t96 * (mrSges(5,1) * t91 - mrSges(5,2) * t90) - t191 * t50 - t192 * t49 - (-t186 * t90 - t193) * t91 / 0.2e1 + t190 * t176 - (-t191 * t91 - t192 * t90) * t119 / 0.2e1 + (-t185 * t91 + t189 - t194) * t179;
t128 = qJD(2) ^ 2;
t118 = Ifges(4,4) * t149;
t116 = pkin(3) * t125 + pkin(4);
t95 = t135 * t143;
t89 = Ifges(4,1) * t151 + t118 + t156;
t88 = t155 + (t157 + t160) * qJD(2);
t83 = t100 * t124;
t80 = pkin(4) * t131 + t117;
t67 = pkin(3) * t151 + t174;
t43 = pkin(4) * t60 + t142;
t42 = -qJ(5) * t131 + t65;
t41 = -qJ(5) * t100 + t64;
t30 = pkin(4) * t50 + t97;
t19 = -qJD(2) * t130 + t119 * t84;
t18 = -t124 * t60 - t127 * t90;
t17 = -t85 + t34;
t16 = t33 + t159;
t14 = t32 - t159;
t11 = mrSges(5,1) * t50 - mrSges(5,2) * t49;
t1 = [t161 * t19 + t162 * t18 + (-t128 * mrSges(3,2) - qJD(2) * t132 - t10 - t11 - t95) * t127 + (-t128 * mrSges(3,1) - t133 * qJD(3) + (qJD(2) * (-mrSges(4,1) * t126 + mrSges(4,2) * t123) - t182) * qJD(2)) * t124 + m(4) * (t134 * t124 + (t154 + (t137 - t153) * t127) * qJD(2)) + m(5) * (-t127 * t97 + t150 * t96 + t18 * t32 + t19 * t31 - t5 * t84 - t6 * t83) + m(6) * (t12 * t19 - t127 * t30 + t14 * t18 + t150 * t58 - t2 * t84 - t3 * t83) + (mrSges(6,3) + mrSges(5,3)) * (-t49 * t83 + t50 * t84); t187 * t70 + (t2 * t42 + t3 * t41 + t30 * t80 + (-t153 + t43) * t58 + t188 * t14 + t187 * t12) * m(6) + t188 * t68 + t183 * t71 + (t117 * t97 + t5 * t65 + t6 * t64 + (t142 - t153) * t96 + t184 * t32 + t183 * t31) * m(5) + t184 * t69 + (t182 * t124 + t132 * t127) * qJD(1) + (-(t127 * t137 + t154) * qJD(1) - pkin(2) * t115 + pkin(6) * t134) * m(4) - (-t163 * t59 - t185 * t60) * t90 / 0.2e1 + (-t163 * t60 - t186 * t59) * t176 + t96 * (mrSges(5,1) * t60 - mrSges(5,2) * t59) + t58 * (mrSges(6,1) * t60 - mrSges(6,2) * t59) + t117 * t11 - pkin(2) * t95 + (0.3e1 / 0.2e1 * t121 - 0.3e1 / 0.2e1 * t120) * Ifges(4,4) * t143 + t80 * t10 + t43 * t53 + (mrSges(5,2) * t97 + mrSges(6,2) * t30 - mrSges(5,3) * t6 - mrSges(6,3) * t3 - t163 * t50 - t186 * t49) * t100 - (-mrSges(5,1) * t97 - mrSges(6,1) * t30 + mrSges(5,3) * t5 + mrSges(6,3) * t2 - t163 * t49 - t185 * t50) * t131 + (t12 * t59 - t14 * t60 + t41 * t49 - t42 * t50) * mrSges(6,3) + (t31 * t59 - t32 * t60 + t64 * t49 - t65 * t50) * mrSges(5,3) + t134 * mrSges(4,3) + ((-pkin(6) * t107 + t113 * mrSges(4,2) + t89 / 0.2e1 + t156 / 0.2e1) * t126 + (-pkin(6) * t108 + pkin(3) * t54 + t113 * mrSges(4,1) - t88 / 0.2e1 - t155 / 0.2e1 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t149) * t123) * qJD(3) - t189 * t59 / 0.2e1 - t190 * t60 / 0.2e1 + (-t191 * t60 - t192 * t59) * t119 / 0.2e1; t129 + (-t113 * t135 + t88 * t170 + (t157 * t170 + (Ifges(4,1) * t126 - t160) * t171) * qJD(2) + (Ifges(4,5) * t126 / 0.2e1 + Ifges(4,6) * t171) * qJD(3) - (t118 + t89) * t126 / 0.2e1) * qJD(2) + t133 * t112 - t34 * t69 - t16 * t70 - t33 * t71 - t72 * mrSges(4,2) + t73 * mrSges(4,1) - t67 * t53 - t17 * t68 + t32 * t166 - m(5) * (t31 * t33 + t32 * t34) + (-mrSges(6,3) * t158 - t54 * t151 + (t125 * t49 - t158) * mrSges(5,3) + (-t122 * t161 + t125 * t162) * qJD(4) + (t122 * t5 + t125 * t6 + t145 * t32 - t146 * t31 - t151 * t96) * m(5)) * pkin(3) + (t116 * t49 + t14 * t91) * mrSges(6,3) + (t116 * t3 + (-t12 * t146 + t122 * t2 + t14 * t145) * pkin(3) - t12 * t16 - t14 * t17 - t58 * t67) * m(6); t129 - t31 * t69 + t14 * t70 + t32 * t71 - t13 * t68 + pkin(4) * t49 * mrSges(6,3) + (-t58 * t174 - (-t12 + t13) * t14 + t3 * pkin(4)) * m(6) + (t32 * mrSges(5,3) + t14 * mrSges(6,3) - pkin(4) * t53) * t91; t90 * t68 + t91 * t70 + 0.2e1 * (t30 / 0.2e1 + t12 * t176 + t14 * t179) * m(6) + t10;];
tauc = t1(:);
