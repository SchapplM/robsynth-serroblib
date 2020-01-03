% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:09
% EndTime: 2019-12-31 17:20:16
% DurationCPUTime: 2.80s
% Computational Cost: add. (1339->309), mult. (3568->413), div. (0->0), fcn. (1900->4), ass. (0->160)
t195 = qJD(2) / 0.2e1;
t194 = Ifges(4,1) + Ifges(5,1);
t191 = Ifges(4,5) + Ifges(5,4);
t90 = cos(qJ(2));
t142 = qJD(1) * t90;
t81 = qJD(3) - t142;
t125 = Ifges(3,5) * t195;
t193 = Ifges(4,6) - Ifges(5,6);
t134 = qJD(1) * qJD(2);
t88 = sin(qJ(2));
t126 = t88 * t134;
t133 = qJD(2) * qJD(3);
t89 = cos(qJ(3));
t135 = t89 * qJD(2);
t87 = sin(qJ(3));
t138 = qJD(3) * t87;
t44 = t89 * t133 + (t135 * t90 - t138 * t88) * qJD(1);
t137 = qJD(3) * t89;
t139 = qJD(2) * t90;
t45 = t87 * t133 + (t137 * t88 + t139 * t87) * qJD(1);
t192 = (-Ifges(4,4) + Ifges(5,5)) * t45 + t194 * t44 + t191 * t126;
t103 = pkin(3) * t87 - qJ(4) * t89;
t85 = pkin(5) * t142;
t190 = -qJD(4) * t87 + t81 * t103 - t85;
t155 = Ifges(5,5) * t87;
t160 = Ifges(4,4) * t87;
t189 = t194 * t89 + t155 - t160;
t154 = Ifges(5,5) * t89;
t159 = Ifges(4,4) * t89;
t188 = t194 * t87 - t154 + t159;
t143 = qJD(1) * t88;
t83 = Ifges(3,4) * t142;
t74 = -pkin(2) * t90 - t88 * pkin(6) - pkin(1);
t63 = t74 * qJD(1);
t78 = qJD(2) * pkin(6) + t85;
t30 = t63 * t89 - t78 * t87;
t31 = t63 * t87 + t78 * t89;
t101 = t30 * t89 + t31 * t87;
t184 = qJD(4) - t30;
t17 = -pkin(3) * t81 + t184;
t18 = qJ(4) * t81 + t31;
t102 = t17 * t89 - t18 * t87;
t106 = Ifges(5,3) * t87 + t154;
t110 = -Ifges(4,2) * t87 + t159;
t115 = mrSges(5,1) * t87 - mrSges(5,3) * t89;
t117 = mrSges(4,1) * t87 + mrSges(4,2) * t89;
t152 = Ifges(5,6) * t87;
t153 = Ifges(4,6) * t87;
t157 = Ifges(4,5) * t89;
t158 = Ifges(5,4) * t89;
t168 = t89 / 0.2e1;
t170 = t87 / 0.2e1;
t171 = -t87 / 0.2e1;
t131 = t89 * t143;
t70 = qJD(2) * t87 + t131;
t174 = t70 / 0.2e1;
t132 = t87 * t143;
t69 = t132 - t135;
t176 = t69 / 0.2e1;
t177 = -t69 / 0.2e1;
t156 = Ifges(5,5) * t69;
t67 = Ifges(4,4) * t69;
t185 = t191 * t81 + t194 * t70 + t156 - t67;
t186 = t81 / 0.2e1;
t77 = -qJD(2) * pkin(2) + pkin(5) * t143;
t19 = pkin(3) * t69 - qJ(4) * t70 + t77;
t66 = Ifges(5,5) * t70;
t20 = Ifges(5,6) * t81 + Ifges(5,3) * t69 + t66;
t161 = Ifges(4,4) * t70;
t23 = -Ifges(4,2) * t69 + Ifges(4,6) * t81 + t161;
t91 = t102 * mrSges(5,2) - t101 * mrSges(4,3) + t106 * t176 + t110 * t177 + t19 * t115 + t77 * t117 + t20 * t170 + t23 * t171 + t189 * t174 + (t152 + t158 - t153 + t157) * t186 + t185 * t168;
t187 = t91 + Ifges(3,1) * t143 / 0.2e1 + t125 + t83 / 0.2e1;
t124 = -Ifges(3,6) * qJD(2) / 0.2e1;
t183 = t191 * t87 + t193 * t89;
t136 = qJD(3) * t90;
t140 = qJD(2) * t88;
t121 = pkin(2) * t88 - pkin(6) * t90;
t72 = t121 * qJD(2);
t16 = pkin(5) * (-t136 * t89 + t140 * t87) - t138 * t74 + t89 * t72;
t64 = qJD(1) * t72;
t182 = qJD(3) * t31 - t64 * t89;
t128 = -Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t129 = Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1;
t130 = Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t162 = Ifges(3,4) * t88;
t181 = t128 * t69 - t129 * t81 - t130 * t70 - t18 * mrSges(5,3) - t30 * mrSges(4,1) - Ifges(4,6) * t177 - Ifges(5,6) * t176 - t124 + (t90 * Ifges(3,2) + t162) * qJD(1) / 0.2e1 + t17 * mrSges(5,1) + t31 * mrSges(4,2) - (Ifges(4,3) + Ifges(5,2)) * t186 - t191 * t174;
t180 = t44 / 0.2e1;
t179 = -t45 / 0.2e1;
t178 = t45 / 0.2e1;
t175 = -t70 / 0.2e1;
t173 = -t81 / 0.2e1;
t169 = -t89 / 0.2e1;
t167 = pkin(1) * mrSges(3,1);
t166 = pkin(1) * mrSges(3,2);
t165 = pkin(5) * t90;
t164 = mrSges(4,3) * t69;
t163 = mrSges(4,3) * t70;
t71 = t121 * qJD(1);
t150 = t71 * t89;
t149 = t89 * t74;
t48 = -mrSges(4,2) * t81 - t164;
t51 = -mrSges(5,2) * t69 + mrSges(5,3) * t81;
t147 = -t48 - t51;
t49 = mrSges(4,1) * t81 - t163;
t50 = -mrSges(5,1) * t81 + t70 * mrSges(5,2);
t146 = -t49 + t50;
t54 = t89 * t165 + t87 * t74;
t141 = qJD(2) * mrSges(3,2);
t127 = pkin(5) * t87 + pkin(3);
t123 = pkin(5) * t126;
t122 = m(4) * t77 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t69 + mrSges(4,2) * t70 + mrSges(3,3) * t143;
t100 = t63 * t137 - t138 * t78 + t87 * t64;
t97 = (-pkin(5) * t89 + qJ(4)) * t143;
t1 = qJD(2) * t97 + qJD(4) * t81 + t100;
t98 = t127 * t143;
t3 = -qJD(2) * t98 + t182;
t120 = t1 * t89 + t3 * t87;
t4 = -t123 * t89 + t100;
t5 = t123 * t87 - t182;
t119 = t4 * t89 - t5 * t87;
t118 = mrSges(4,1) * t89 - mrSges(4,2) * t87;
t116 = mrSges(5,1) * t89 + mrSges(5,3) * t87;
t109 = Ifges(4,2) * t89 + t160;
t105 = -Ifges(5,3) * t89 + t155;
t104 = pkin(3) * t89 + qJ(4) * t87;
t99 = pkin(5) + t103;
t28 = -mrSges(5,1) * t126 + t44 * mrSges(5,2);
t94 = -t5 * mrSges(4,1) + t3 * mrSges(5,1) + t4 * mrSges(4,2) - t1 * mrSges(5,3);
t15 = t87 * t72 + t74 * t137 + (-t135 * t88 - t136 * t87) * pkin(5);
t80 = Ifges(5,2) * t126;
t79 = Ifges(4,3) * t126;
t76 = mrSges(3,3) * t142 - t141;
t73 = -pkin(2) - t104;
t61 = t87 * t71;
t56 = t99 * t88;
t53 = -t165 * t87 + t149;
t47 = -pkin(5) * t131 + t61;
t46 = pkin(5) * t132 + t150;
t43 = t127 * t90 - t149;
t42 = -qJ(4) * t90 + t54;
t41 = Ifges(5,4) * t44;
t40 = Ifges(4,5) * t44;
t39 = Ifges(4,6) * t45;
t38 = Ifges(5,6) * t45;
t35 = mrSges(5,1) * t69 - mrSges(5,3) * t70;
t34 = pkin(3) * t70 + qJ(4) * t69;
t33 = -t98 - t150;
t32 = t61 + t97;
t29 = -mrSges(4,2) * t126 - mrSges(4,3) * t45;
t27 = mrSges(4,1) * t126 - mrSges(4,3) * t44;
t26 = -mrSges(5,2) * t45 + mrSges(5,3) * t126;
t14 = (qJD(3) * t104 - qJD(4) * t89) * t88 + t99 * t139;
t13 = -pkin(3) * t140 - t16;
t12 = mrSges(4,1) * t45 + mrSges(4,2) * t44;
t11 = mrSges(5,1) * t45 - mrSges(5,3) * t44;
t10 = qJ(4) * t140 - qJD(4) * t90 + t15;
t7 = t44 * Ifges(4,4) - t45 * Ifges(4,2) + Ifges(4,6) * t126;
t6 = t44 * Ifges(5,5) + Ifges(5,6) * t126 + t45 * Ifges(5,3);
t2 = t45 * pkin(3) - t44 * qJ(4) + qJD(2) * t85 - t70 * qJD(4);
t8 = [t10 * t51 + t56 * t11 + t13 * t50 + t14 * t35 + t15 * t48 + t16 * t49 + t42 * t26 + t53 * t27 + t43 * t28 + t54 * t29 + m(4) * (t31 * t15 + t30 * t16 + t4 * t54 + t5 * t53) + m(5) * (t1 * t42 + t10 * t18 + t13 * t17 + t14 * t19 + t2 * t56 + t3 * t43) + (-t79 / 0.2e1 - t80 / 0.2e1 - t40 / 0.2e1 - t41 / 0.2e1 - t38 / 0.2e1 + t39 / 0.2e1 + t128 * t45 - t130 * t44 + ((0.3e1 / 0.2e1 * Ifges(3,4) * t90 - 0.2e1 * t166) * qJD(1) + t122 * pkin(5) + t125 + t187) * qJD(2) + t94) * t90 + (pkin(5) * t12 + t7 * t171 + t6 * t170 + t2 * t115 + t106 * t178 + t110 * t179 + (-t4 * t87 - t5 * t89) * mrSges(4,3) + (-t1 * t87 + t3 * t89) * mrSges(5,2) + (-pkin(5) * t76 + t124 - t181) * qJD(2) + (t77 * t118 + t109 * t176 + t105 * t177 + t19 * t116 + t23 * t169 + (t30 * t87 - t31 * t89) * mrSges(4,3) + (-t17 * t87 - t18 * t89) * mrSges(5,2) + t188 * t175 + t183 * t173 + t185 * t171) * qJD(3) + ((-0.3e1 / 0.2e1 * Ifges(3,4) + t158 / 0.2e1 + t152 / 0.2e1 + t157 / 0.2e1 - t153 / 0.2e1) * t88 - 0.2e1 * t167 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(5) + t117) * pkin(5) - t129) * t90) * t134 + t189 * t180 + (t20 * qJD(3) + t192) * t168) * t88; t91 * qJD(3) - m(4) * (t30 * t46 + t31 * t47) + t190 * t35 + t73 * t11 - t46 * t49 - t33 * t50 - t32 * t51 - t47 * t48 - pkin(2) * t12 - t2 * t116 + t119 * mrSges(4,3) + t120 * mrSges(5,2) + ((t124 + (t167 + t162 / 0.2e1) * qJD(1) + (t76 + t141) * pkin(5) + t183 * t195 + t181) * t88 + (((-m(4) * pkin(2) - mrSges(3,1) - t118) * qJD(2) - t122) * pkin(5) - t83 / 0.2e1 + (t166 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t88) * qJD(1) + t125 - t187) * t90) * qJD(1) + t7 * t168 + t6 * t169 + t105 * t178 + t109 * t179 + t188 * t180 + t192 * t170 + (-t17 * t33 - t18 * t32 + t190 * t19 + t2 * t73) * m(5) + ((-m(4) * t101 + m(5) * t102 + t146 * t89 + t147 * t87) * qJD(3) + (t26 + t29) * t89 + (-t27 + t28) * t87 + m(5) * t120 + m(4) * t119) * pkin(6); (-t146 + t163) * t31 + (t147 - t164) * t30 + (t17 * t69 + t18 * t70) * mrSges(5,2) - t94 + t79 + t80 - t77 * (mrSges(4,1) * t70 - mrSges(4,2) * t69) - t19 * (mrSges(5,1) * t70 + mrSges(5,3) * t69) + qJD(4) * t51 - t34 * t35 + t40 + t41 + t38 - t39 + qJ(4) * t26 - pkin(3) * t28 + t23 * t174 + (Ifges(5,3) * t70 - t156) * t177 + (-t191 * t69 - t193 * t70) * t173 + (-t3 * pkin(3) + t1 * qJ(4) - t17 * t31 + t184 * t18 - t19 * t34) * m(5) + (-Ifges(4,2) * t70 + t185 - t67) * t176 + (-t194 * t69 - t161 + t20 + t66) * t175; t70 * t35 - t81 * t51 + 0.2e1 * (t3 / 0.2e1 + t18 * t173 + t19 * t174) * m(5) + t28;];
tauc = t8(:);
