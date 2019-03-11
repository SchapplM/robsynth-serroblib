% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:18
% EndTime: 2019-03-09 01:35:22
% DurationCPUTime: 2.58s
% Computational Cost: add. (2260->305), mult. (4216->446), div. (0->0), fcn. (2073->6), ass. (0->155)
t78 = cos(qJ(5));
t146 = qJD(1) * t78;
t124 = mrSges(6,3) * t146;
t77 = cos(qJ(6));
t140 = qJD(5) * t77;
t75 = sin(qJ(6));
t53 = t146 * t75 + t140;
t137 = t75 * qJD(5);
t54 = t146 * t77 - t137;
t151 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t53 - mrSges(7,2) * t54 - t124;
t139 = qJD(5) * t78;
t73 = cos(pkin(9));
t125 = t73 * t139;
t76 = sin(qJ(5));
t157 = t75 * t76;
t72 = sin(pkin(9));
t155 = t76 * t77;
t93 = t155 * t73 - t72 * t75;
t197 = qJD(6) * t93 + t125 * t75 - (-t157 * t72 + t73 * t77) * qJD(1);
t45 = t157 * t73 + t72 * t77;
t196 = qJD(6) * t45 - t125 * t77 - (t155 * t72 + t73 * t75) * qJD(1);
t136 = qJ(2) * qJD(1);
t79 = -pkin(1) - pkin(2);
t61 = qJD(1) * t79 + qJD(2);
t42 = -t72 * t136 + t61 * t73;
t37 = qJD(1) * pkin(3) + qJD(4) - t42;
t29 = qJD(1) * pkin(7) + t37;
t22 = qJD(3) * t78 + t29 * t76;
t20 = qJD(5) * pkin(8) + t22;
t110 = -pkin(5) * t76 + pkin(8) * t78;
t160 = t61 * t72;
t68 = t73 * qJ(2);
t71 = qJD(1) * qJ(4);
t25 = t160 - t71 + (t110 + t68) * qJD(1);
t5 = -t20 * t75 + t25 * t77;
t6 = t20 * t77 + t25 * t75;
t107 = t5 * t77 + t6 * t75;
t147 = qJD(1) * t76;
t63 = qJD(6) - t147;
t35 = -mrSges(7,2) * t63 + mrSges(7,3) * t53;
t36 = mrSges(7,1) * t63 + mrSges(7,3) * t54;
t83 = -m(7) * t107 - t75 * t35 - t77 * t36;
t141 = qJD(5) * t76;
t21 = -qJD(3) * t76 + t29 * t78;
t135 = qJD(1) * qJD(2);
t121 = t72 * t135;
t143 = qJD(5) * t21;
t12 = t121 * t76 + t143;
t142 = qJD(5) * t22;
t13 = -t121 * t78 + t142;
t165 = t13 * t78;
t98 = t12 * t76 - t165;
t195 = t22 * t139 - t21 * t141 + t98;
t133 = qJD(5) * qJD(6);
t138 = qJD(6) * t78;
t88 = t138 * t75 + t140 * t76;
t33 = qJD(1) * t88 + t133 * t77;
t182 = t33 / 0.2e1;
t87 = t137 * t76 - t138 * t77;
t34 = -qJD(1) * t87 - t133 * t75;
t181 = t34 / 0.2e1;
t193 = -mrSges(5,2) + mrSges(4,1);
t192 = m(3) * qJ(2) + mrSges(3,3);
t145 = qJD(2) * t72;
t114 = -t72 * qJ(2) + t73 * t79;
t112 = pkin(3) - t114;
t51 = pkin(7) + t112;
t191 = t51 * t139 + t76 * t145;
t70 = qJD(1) * qJD(4);
t111 = -pkin(5) * t78 - pkin(8) * t76;
t144 = qJD(2) * t73;
t85 = qJD(5) * t111 + t144;
t32 = qJD(1) * t85 - t70;
t1 = qJD(6) * t5 + t12 * t77 + t32 * t75;
t2 = -qJD(6) * t6 - t12 * t75 + t32 * t77;
t108 = t1 * t77 - t2 * t75;
t100 = Ifges(7,5) * t77 - Ifges(7,6) * t75;
t167 = Ifges(7,4) * t77;
t102 = -Ifges(7,2) * t75 + t167;
t168 = Ifges(7,4) * t75;
t104 = Ifges(7,1) * t77 - t168;
t105 = mrSges(7,1) * t75 + mrSges(7,2) * t77;
t169 = Ifges(7,4) * t54;
t17 = Ifges(7,2) * t53 + Ifges(7,6) * t63 - t169;
t122 = -t75 * t17 / 0.2e1;
t50 = Ifges(7,4) * t53;
t18 = -Ifges(7,1) * t54 + Ifges(7,5) * t63 + t50;
t154 = t77 * t18;
t177 = t63 / 0.2e1;
t178 = -t54 / 0.2e1;
t179 = t53 / 0.2e1;
t19 = -qJD(5) * pkin(5) - t21;
t188 = t107 * mrSges(7,3) - t154 / 0.2e1 - t122 - t102 * t179 - t104 * t178 - t19 * t105 - t100 * t177;
t11 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t187 = -m(6) * (-t13 + t142) + m(7) * (t137 * t5 - t140 * t6 + t13) + t11;
t185 = m(6) / 0.2e1;
t184 = m(7) / 0.2e1;
t134 = qJD(1) * qJD(5);
t119 = t78 * t134;
t183 = Ifges(7,4) * t182 + Ifges(7,2) * t181 - Ifges(7,6) * t119 / 0.2e1;
t171 = Ifges(6,4) * t76;
t170 = Ifges(6,4) * t78;
t164 = t22 * t76;
t43 = t136 * t73 + t160;
t38 = t43 - t71;
t163 = t38 * t73;
t162 = t53 * Ifges(7,6);
t161 = t54 * Ifges(7,5);
t159 = t63 * Ifges(7,3);
t158 = t72 * t78;
t156 = t75 * t78;
t153 = t77 * t78;
t152 = t78 * t11;
t150 = t72 * t79 + t68;
t149 = Ifges(6,5) * qJD(5);
t148 = Ifges(6,6) * qJD(5);
t52 = -qJ(4) + t150;
t132 = mrSges(6,3) * t147;
t120 = t73 * t135;
t109 = -qJD(6) * t51 * t76 - qJD(4) + t85;
t106 = -mrSges(6,1) * t76 - mrSges(6,2) * t78;
t103 = Ifges(7,1) * t75 + t167;
t101 = Ifges(7,2) * t77 + t168;
t99 = Ifges(7,5) * t75 + Ifges(7,6) * t77;
t23 = -mrSges(7,1) * t119 - mrSges(7,3) * t33;
t24 = mrSges(7,2) * t119 + mrSges(7,3) * t34;
t97 = -t75 * t23 + t77 * t24;
t57 = -t70 + t120;
t62 = -qJD(4) + t144;
t95 = t38 * t62 + t57 * t52;
t94 = -t42 * t72 + t43 * t73;
t58 = -qJD(5) * mrSges(6,2) + t132;
t92 = t35 * t77 - t36 * t75 + t58;
t91 = (-mrSges(6,1) * t78 + mrSges(6,2) * t76) * qJD(5);
t90 = (t76 * Ifges(6,2) - t170) * qJD(1);
t89 = Ifges(7,5) * t33 + Ifges(7,6) * t34 - Ifges(7,3) * t119;
t39 = t110 + t52;
t84 = qJD(6) * t39 + t191;
t82 = m(6) * (t12 - t143) + m(7) * (qJD(5) * t19 + t108) + t97 + t83 * qJD(6);
t80 = qJD(1) ^ 2;
t66 = Ifges(6,4) * t147;
t56 = t111 * qJD(1);
t55 = t106 * qJD(1);
t49 = qJD(1) * t91;
t48 = -Ifges(6,1) * t146 + t149 + t66;
t47 = t90 + t148;
t16 = t159 - t161 + t162;
t15 = t155 * t51 + t39 * t75;
t14 = -t157 * t51 + t39 * t77;
t10 = t21 * t77 + t56 * t75;
t9 = -t21 * t75 + t56 * t77;
t8 = t33 * Ifges(7,1) + t34 * Ifges(7,4) - Ifges(7,5) * t119;
t4 = t109 * t77 - t75 * t84;
t3 = t109 * t75 + t77 * t84;
t7 = [(-qJD(1) * t62 - t57) * mrSges(5,3) + 0.2e1 * t193 * t121 + (-m(7) * t19 - t151) * (-t51 * t141 + t78 * t145) + 0.2e1 * t192 * t135 + t191 * t58 + t195 * mrSges(6,3) + t62 * t55 + t52 * t49 + t4 * t36 + qJD(5) ^ 2 * (Ifges(6,5) * t76 + Ifges(6,6) * t78) / 0.2e1 + t3 * t35 + t14 * t23 + t15 * t24 + m(7) * (t1 * t15 + t14 * t2 - t51 * t165 + t3 * t6 + t4 * t5) - (Ifges(6,1) * t76 + t170) * t119 - t105 * t165 + m(6) * ((t21 * t78 + t164) * t145 + ((-t21 * t76 + t22 * t78) * qJD(5) + t98) * t51 + t95) - t51 * t152 - t8 * t153 / 0.2e1 + t2 * (-mrSges(7,1) * t76 + mrSges(7,3) * t153) + t1 * (mrSges(7,2) * t76 + mrSges(7,3) * t156) + m(5) * ((qJD(1) * t112 + t37) * t145 + t95) + t5 * (-mrSges(7,1) * t139 - mrSges(7,3) * t88) + t6 * (mrSges(7,2) * t139 - mrSges(7,3) * t87) + 0.2e1 * mrSges(4,2) * t120 + t122 * t141 + m(4) * ((-t114 * t72 + t150 * t73) * qJD(1) + t94) * qJD(2) + t76 * (Ifges(6,2) * t78 + t171) * t134 + t38 * t91 + (t99 * t138 + (-Ifges(7,3) * t78 + t100 * t76) * qJD(5)) * t177 + (t103 * t138 + (-Ifges(7,5) * t78 + t104 * t76) * qJD(5)) * t178 + (t101 * t138 + (-Ifges(7,6) * t78 + t102 * t76) * qJD(5)) * t179 + (-Ifges(7,6) * t76 - t102 * t78) * t181 + (-Ifges(7,5) * t76 - t104 * t78) * t182 + t156 * t183 + t19 * (mrSges(7,1) * t87 + mrSges(7,2) * t88) - t76 * t89 / 0.2e1 + t57 * t106 + (t77 * t17 + t75 * t18) * t138 / 0.2e1 + (t90 + t47) * t139 / 0.2e1 - (t16 + qJD(1) * (-Ifges(7,3) * t76 - t100 * t78)) * t139 / 0.2e1 + (t48 + qJD(1) * (-t78 * Ifges(6,1) + t171) + t154) * t141 / 0.2e1; t45 * t23 - t93 * t24 - t192 * t80 + t197 * t36 + t196 * t35 + (t152 + (-mrSges(4,2) + mrSges(5,3)) * t80 + (-t151 * t76 - t78 * t58) * qJD(5)) * t73 + 0.2e1 * (-t195 * t185 + (-t141 * t19 + t165) * t184) * t73 + (-t193 * t80 + t49 + 0.2e1 * (t185 + m(5) / 0.2e1) * t57) * t72 + (-t1 * t93 + t196 * t6 + t197 * t5 + t2 * t45) * m(7) + ((t151 * t78 - t76 * t58) * t72 - t55 * t73 - m(6) * (t158 * t21 + t164 * t72 + t163) + 0.2e1 * t19 * t158 * t184 - m(4) * t94 - (t163 + (t144 + t37) * t72) * m(5)) * qJD(1); ((-t92 + t132) * qJD(5) + t187) * t76 + ((t124 + t151) * qJD(5) + t82) * t78; -t80 * mrSges(5,3) + (qJD(5) * t92 - t187) * t78 + (qJD(5) * t151 + t82) * t76 + (m(6) * t38 + t55 + (t38 + t145) * m(5) - t83) * qJD(1); t103 * t182 + t101 * t181 + t77 * t183 + t75 * t8 / 0.2e1 - t21 * t58 - t10 * t35 - t9 * t36 - pkin(5) * t11 - t12 * mrSges(6,2) - t151 * t22 + (-mrSges(7,1) * t77 + mrSges(7,2) * t75 - mrSges(6,1)) * t13 + t97 * pkin(8) + t108 * mrSges(7,3) + (pkin(8) * t83 - t188) * qJD(6) + ((-qJD(5) * t99 / 0.2e1 - t22 * mrSges(6,3) + t148 / 0.2e1 - t47 / 0.2e1 + t16 / 0.2e1 - t6 * mrSges(7,2) + t38 * mrSges(6,1) + t159 / 0.2e1 + t162 / 0.2e1 - t161 / 0.2e1 + t5 * mrSges(7,1) + Ifges(6,4) * t146 / 0.2e1) * t78 + (-t66 / 0.2e1 - t48 / 0.2e1 - t38 * mrSges(6,2) + t149 / 0.2e1 + t21 * mrSges(6,3) + (-Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t146 + t188) * t76) * qJD(1) + (-pkin(5) * t13 + pkin(8) * t108 - t10 * t6 - t19 * t22 - t5 * t9) * m(7); -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t19 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t54 * (Ifges(7,1) * t53 + t169) / 0.2e1 + t17 * t178 - t63 * (Ifges(7,5) * t53 + Ifges(7,6) * t54) / 0.2e1 - t5 * t35 + t6 * t36 + (t5 * t53 - t54 * t6) * mrSges(7,3) + t89 - (Ifges(7,2) * t54 + t18 + t50) * t53 / 0.2e1;];
tauc  = t7(:);
