% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:00
% EndTime: 2019-12-05 16:12:10
% DurationCPUTime: 3.19s
% Computational Cost: add. (1206->322), mult. (3274->454), div. (0->0), fcn. (1854->6), ass. (0->150)
t195 = Ifges(5,1) + Ifges(6,1);
t190 = Ifges(6,4) + Ifges(5,5);
t101 = cos(qJ(3));
t140 = qJD(2) * t101;
t194 = t140 / 0.2e1;
t193 = qJD(2) / 0.2e1;
t192 = qJD(3) / 0.2e1;
t191 = -Ifges(5,4) + Ifges(6,5);
t97 = sin(pkin(8));
t169 = Ifges(6,5) * t97;
t173 = Ifges(5,4) * t97;
t98 = cos(pkin(8));
t189 = t195 * t98 + t169 - t173;
t99 = sin(qJ(3));
t147 = qJD(2) * t99;
t73 = t97 * qJD(3) + t147 * t98;
t178 = t73 / 0.2e1;
t136 = qJD(2) * qJD(3);
t187 = t136 / 0.2e1;
t185 = t147 / 0.2e1;
t125 = -Ifges(4,6) * qJD(3) / 0.2e1;
t184 = Ifges(4,5) * t192;
t100 = sin(qJ(2));
t143 = qJD(1) * t100;
t88 = qJD(2) * pkin(6) + t143;
t183 = (t101 ^ 2 + t99 ^ 2) * t88;
t102 = cos(qJ(2));
t142 = qJD(1) * t102;
t79 = -pkin(3) * t101 - qJ(4) * t99 - pkin(2);
t53 = qJD(2) * t79 - t142;
t77 = t101 * t88;
t67 = qJD(3) * qJ(4) + t77;
t10 = t97 * t53 + t98 * t67;
t168 = Ifges(6,5) * t98;
t113 = Ifges(6,3) * t97 + t168;
t172 = Ifges(5,4) * t98;
t114 = -Ifges(5,2) * t97 + t172;
t117 = mrSges(5,1) * t97 + mrSges(5,2) * t98;
t175 = t98 / 0.2e1;
t176 = t97 / 0.2e1;
t177 = -t97 / 0.2e1;
t9 = t53 * t98 - t67 * t97;
t6 = pkin(4) * t140 + qJD(5) - t9;
t76 = t99 * t88;
t63 = -qJD(3) * pkin(3) + qJD(4) + t76;
t7 = -qJ(5) * t140 + t10;
t72 = -t98 * qJD(3) + t147 * t97;
t8 = pkin(4) * t72 - qJ(5) * t73 + t63;
t89 = -qJD(2) * pkin(2) - t142;
t93 = Ifges(4,4) * t140;
t182 = (t6 * t98 - t7 * t97) * mrSges(6,2) - (t10 * t97 + t9 * t98) * mrSges(5,3) + t63 * t117 + t8 * (mrSges(6,1) * t97 - mrSges(6,3) * t98) + t89 * mrSges(4,2) + Ifges(4,1) * t185 + t184 + t93 / 0.2e1 + (t73 * Ifges(6,5) - Ifges(6,6) * t140) * t176 + (t73 * Ifges(5,4) - Ifges(5,6) * t140) * t177 + t189 * t178 + (-t190 * t140 + t195 * t73) * t175 + (Ifges(6,3) * t176 - Ifges(5,2) * t177 - t114 / 0.2e1 + t113 / 0.2e1 + t191 * t175) * t72;
t181 = 0.2e1 * (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t72 - (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t73 - t7 * mrSges(6,3) - t89 * mrSges(4,1) - t9 * mrSges(5,1) - t125 + (Ifges(4,4) * t99 + t101 * Ifges(4,2)) * t193 + t10 * mrSges(5,2) + t6 * mrSges(6,1) + (Ifges(5,3) + Ifges(6,2)) * t194 - t190 * t178;
t112 = pkin(4) * t97 - qJ(5) * t98;
t108 = t112 * t101;
t138 = qJD(3) * t101;
t127 = t88 * t138;
t5 = t127 - qJD(5) * t73 + (qJD(3) * t108 + t142 * t99) * qJD(2);
t179 = m(6) * t5;
t139 = qJD(2) * t102;
t124 = qJD(1) * t139;
t43 = t124 * t99 + t127;
t174 = m(5) * t43;
t171 = Ifges(6,4) * t98;
t170 = Ifges(5,5) * t98;
t167 = Ifges(5,6) * t97;
t166 = Ifges(6,6) * t97;
t165 = t43 * t99;
t111 = pkin(3) * t99 - qJ(4) * t101;
t62 = qJD(3) * t111 - qJD(4) * t99;
t164 = t62 * t98;
t75 = t111 * qJD(2);
t163 = t75 * t98;
t162 = t79 * t98;
t86 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t140;
t161 = t99 * t86;
t158 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t72 + mrSges(5,2) * t73 + mrSges(4,3) * t147;
t37 = (t62 + t143) * qJD(2);
t87 = t101 * t124;
t40 = t87 + (qJD(4) - t76) * qJD(3);
t4 = t97 * t37 + t98 * t40;
t46 = -mrSges(6,2) * t72 - mrSges(6,3) * t140;
t47 = mrSges(5,2) * t140 - mrSges(5,3) * t72;
t157 = t46 + t47;
t48 = -mrSges(5,1) * t140 - mrSges(5,3) * t73;
t49 = mrSges(6,1) * t140 + mrSges(6,2) * t73;
t156 = t48 - t49;
t123 = t101 * t136;
t120 = t98 * t123;
t121 = t97 * t123;
t50 = mrSges(6,1) * t121 - mrSges(6,3) * t120;
t51 = mrSges(5,1) * t121 + mrSges(5,2) * t120;
t155 = t50 + t51;
t152 = t101 * t97;
t58 = (-mrSges(5,2) * t99 - mrSges(5,3) * t152) * t136;
t61 = (-mrSges(6,2) * t152 + mrSges(6,3) * t99) * t136;
t154 = t58 + t61;
t151 = t101 * t98;
t59 = (mrSges(5,1) * t99 - mrSges(5,3) * t151) * t136;
t126 = t99 * t136;
t83 = mrSges(6,2) * t120;
t60 = -mrSges(6,1) * t126 + t83;
t153 = -t59 + t60;
t45 = pkin(6) * t151 + t97 * t79;
t146 = qJD(3) * t99;
t145 = t100 * t101;
t144 = t101 * t102;
t141 = qJD(2) * t100;
t137 = qJD(5) * t101;
t24 = mrSges(6,1) * t72 - mrSges(6,3) * t73;
t135 = t24 + t158;
t134 = pkin(6) * t146;
t133 = t97 * t144;
t130 = pkin(6) * t97 + pkin(4);
t129 = t99 * t139;
t128 = t100 * t146;
t122 = m(5) * t63 + t158;
t3 = t37 * t98 - t40 * t97;
t110 = pkin(6) + t112;
t109 = t100 * t97 + t144 * t98;
t1 = (qJ(5) * t146 - t137) * qJD(2) + t4;
t107 = t5 * mrSges(6,1) + (Ifges(6,6) * t99 + t101 * t113) * t187 - (Ifges(5,6) * t99 + t101 * t114) * t136 / 0.2e1 - t1 * mrSges(6,2) - t4 * mrSges(5,3);
t2 = -pkin(4) * t126 - t3;
t106 = t2 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (t189 * t101 + t190 * t99) * t187;
t103 = qJD(2) ^ 2;
t74 = (-mrSges(4,1) * t101 + mrSges(4,2) * t99) * qJD(2);
t70 = (mrSges(4,1) * t99 + mrSges(4,2) * t101) * t136;
t66 = t97 * t75;
t65 = -t102 * t97 + t145 * t98;
t64 = t102 * t98 + t145 * t97;
t57 = t109 * qJD(1);
t56 = qJD(1) * t133 - t143 * t98;
t54 = t110 * t99;
t52 = t97 * t62;
t44 = -pkin(6) * t152 + t162;
t42 = -t146 * t88 + t87;
t39 = t101 * t130 - t162;
t38 = -qJ(5) * t101 + t45;
t32 = qJD(2) * t109 - t128 * t98;
t31 = qJD(2) * t133 - t128 * t97 - t141 * t98;
t29 = -t134 * t98 + t52;
t28 = t134 * t97 + t164;
t27 = qJD(2) * t108 + t77;
t23 = -t76 * t98 + t66;
t22 = t76 * t97 + t163;
t21 = -qJD(5) * t98 * t99 + t110 * t138;
t20 = -t130 * t146 - t164;
t13 = -t163 + (-pkin(4) * qJD(2) - t88 * t97) * t99;
t12 = t66 + (qJ(5) * qJD(2) - t88 * t98) * t99;
t11 = -t137 + t52 + (-pkin(6) * t98 + qJ(5)) * t146;
t14 = [t154 * t65 + t153 * t64 + t157 * t32 - t156 * t31 + (-t103 * mrSges(3,2) - t70 + (t101 * t86 + t135 * t99) * qJD(2)) * t102 + m(6) * (t1 * t65 + t8 * t129 + t2 * t64 + t31 * t6 + t32 * t7) + m(5) * (t10 * t32 + t63 * t129 - t3 * t64 - t31 * t9 + t4 * t65) + m(4) * t139 * t183 + (-t103 * mrSges(3,1) + qJD(2) * t74 + t155 * t99 + (t101 * t135 - t161) * qJD(3) + m(4) * (qJD(2) * t89 + t101 * t42 - t124 + t165) + m(6) * (t138 * t8 + t5 * t99) + m(5) * (t138 * t63 + t165)) * t100; -t74 * t143 - pkin(2) * t70 + t21 * t24 + t38 * t61 + t39 * t60 + t44 * t59 + t45 * t58 + t54 * t50 + (t20 - t56) * t49 + (t28 + t56) * t48 + (t29 - t57) * t47 + (t11 - t57) * t46 - m(6) * (t56 * t6 + t57 * t7) - m(5) * (t10 * t57 - t56 * t9) + m(6) * (t1 * t38 + t11 * t7 + t2 * t39 + t20 * t6 + t21 * t8 + t5 * t54) + m(5) * (t10 * t29 + t28 * t9 + t3 * t44 + t4 * t45) - (pkin(2) * t141 + t100 * t89 + t102 * t183) * m(4) * qJD(1) + (t106 * t98 + t107 * t97 + (mrSges(4,3) + t117) * t43 + (mrSges(4,2) * t141 + (-m(6) * t8 - t122 - t24) * t102) * qJD(1) + ((-0.3e1 / 0.2e1 * qJD(2) * Ifges(4,4) + (t166 + t171 - t167 + t170) * t193) * t99 + t125 - t181) * qJD(3) + (t51 + (m(5) + m(4)) * t43 - t86 * qJD(3)) * pkin(6)) * t99 + (-t3 * mrSges(5,1) + t2 * mrSges(6,1) + t4 * mrSges(5,2) - t1 * mrSges(6,3) + (m(4) * pkin(6) + mrSges(4,3)) * t42 + (-mrSges(4,1) * t141 - t102 * t86) * qJD(1) + (((0.3e1 / 0.2e1 * Ifges(4,4) + (-0.3e1 / 0.2e1 * Ifges(6,4) - 0.3e1 / 0.2e1 * Ifges(5,5)) * t98 + (-0.3e1 / 0.2e1 * Ifges(6,6) + 0.3e1 / 0.2e1 * Ifges(5,6)) * t97) * t101 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t98 ^ 2 + ((Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t97 + t191 * t98) * t97) * t99) * qJD(2) + t122 * pkin(6) + t184 + t182) * qJD(3)) * t101; -t43 * mrSges(4,1) - t42 * mrSges(4,2) - t12 * t46 - t13 * t49 - t22 * t48 - t23 * t47 - t27 * t24 + (-t101 * t158 + t161) * t88 - m(5) * (t10 * t23 + t22 * t9 + t63 * t77) - m(6) * (t12 * t7 + t13 * t6 + t27 * t8) + (-t43 * mrSges(5,1) + t157 * qJD(4) + t154 * qJ(4) + m(6) * (qJ(4) * t1 + qJD(4) * t7) + m(5) * (qJ(4) * t4 + qJD(4) * t10) - t107) * t98 + (t43 * mrSges(5,2) - qJD(5) * t24 - t156 * qJD(4) + t153 * qJ(4) + m(6) * (qJ(4) * t2 + qJD(4) * t6 - qJD(5) * t8) + m(5) * (-qJ(4) * t3 - qJD(4) * t9) + t106) * t97 + ((t125 + Ifges(4,4) * t185 + ((Ifges(5,6) - Ifges(6,6)) * t98 + t190 * t97) * t192 + t181) * t99 + ((Ifges(4,5) / 0.2e1 + (Ifges(5,2) * t98 + t173) * t177 + (-Ifges(6,3) * t98 + t169) * t176 + (t195 * t97 - t168 + t172) * t175) * qJD(3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t147 + (t170 / 0.2e1 - t167 / 0.2e1 + t171 / 0.2e1 + t166 / 0.2e1) * t140 - t93 / 0.2e1 - t182) * t101) * qJD(2) + (t50 + t179) * (-pkin(4) * t98 - qJ(5) * t97 - pkin(3)) + (-t51 - t174) * pkin(3); t156 * t73 + t157 * t72 + t174 + t179 - m(5) * (-t10 * t72 - t9 * t73) - m(6) * (t6 * t73 - t7 * t72) + t155; t73 * t24 + t83 + (-mrSges(6,1) * t146 + t101 * t46) * qJD(2) + 0.2e1 * (t2 / 0.2e1 + t7 * t194 + t8 * t178) * m(6);];
tauc = t14(:);
