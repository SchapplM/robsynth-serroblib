% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:47
% EndTime: 2019-03-09 02:58:53
% DurationCPUTime: 2.84s
% Computational Cost: add. (2177->368), mult. (4412->463), div. (0->0), fcn. (1937->4), ass. (0->172)
t89 = sin(qJ(6));
t91 = cos(qJ(6));
t119 = mrSges(7,1) * t89 + mrSges(7,2) * t91;
t93 = -pkin(3) - pkin(4);
t86 = -pkin(8) + t93;
t90 = sin(qJ(3));
t92 = cos(qJ(3));
t103 = pkin(5) * t92 + t86 * t90 - qJ(2);
t151 = qJD(1) * t92;
t80 = qJ(4) * t151;
t141 = qJD(5) + t80;
t21 = qJD(1) * t103 + t141;
t140 = qJ(5) * qJD(1);
t94 = -pkin(1) - pkin(7);
t75 = qJD(1) * t94 + qJD(2);
t45 = (t75 + t140) * t92;
t198 = t45 - qJD(4);
t26 = qJD(3) * t86 - t198;
t5 = t21 * t91 - t26 * t89;
t6 = t21 * t89 + t26 * t91;
t123 = t5 * t91 + t6 * t89;
t142 = t89 * qJD(3);
t152 = qJD(1) * t90;
t60 = t152 * t91 - t142;
t175 = Ifges(7,4) * t60;
t147 = qJD(3) * t91;
t59 = -t152 * t89 - t147;
t77 = qJD(6) + t151;
t16 = Ifges(7,2) * t59 + Ifges(7,6) * t77 + t175;
t58 = Ifges(7,4) * t59;
t17 = Ifges(7,1) * t60 + Ifges(7,5) * t77 + t58;
t183 = -t91 / 0.2e1;
t184 = -t89 / 0.2e1;
t65 = t90 * t75;
t56 = qJD(3) * qJ(4) + t65;
t79 = t90 * t140;
t40 = -t79 - t56;
t32 = qJD(3) * pkin(5) - t40;
t207 = -t123 * mrSges(7,3) + t32 * t119 + t16 * t184 - t17 * t183;
t170 = Ifges(7,6) * t89;
t171 = Ifges(7,5) * t91;
t114 = -t170 + t171;
t173 = Ifges(7,4) * t91;
t116 = -Ifges(7,2) * t89 + t173;
t174 = Ifges(7,4) * t89;
t118 = Ifges(7,1) * t91 - t174;
t185 = t77 / 0.2e1;
t187 = t60 / 0.2e1;
t189 = t59 / 0.2e1;
t201 = -qJD(3) / 0.2e1;
t202 = -qJD(1) / 0.2e1;
t203 = Ifges(5,6) / 0.2e1;
t121 = t90 * t93 - qJ(2);
t39 = qJD(1) * t121 + t141;
t133 = pkin(3) * t90 + qJ(2);
t49 = qJD(1) * t133 - t80;
t81 = Ifges(5,5) * t151;
t206 = t39 * mrSges(6,2) + t49 * mrSges(5,1) + qJD(3) * t203 + Ifges(5,3) * t152 / 0.2e1 + t81 / 0.2e1 + (Ifges(4,4) * t92 - t90 * Ifges(4,2)) * t202 + (t90 * Ifges(6,1) - Ifges(6,4) * t92) * qJD(1) / 0.2e1 - t40 * mrSges(6,3) - t56 * mrSges(5,2) + t116 * t189 + t118 * t187 + t114 * t185 + (Ifges(4,6) + Ifges(6,5)) * t201 + t207;
t205 = (qJ(2) * (m(3) + m(4)) + mrSges(3,3)) * qJD(1);
t204 = (-mrSges(4,1) - mrSges(5,1)) * qJD(3);
t62 = (mrSges(6,1) * t92 + mrSges(6,2) * t90) * qJD(1);
t63 = (mrSges(5,1) * t90 - mrSges(5,3) * t92) * qJD(1);
t200 = -t62 + t63;
t199 = -t5 * t89 + t6 * t91;
t84 = t92 * qJD(4);
t78 = qJD(1) * t84;
t88 = qJ(4) + pkin(5);
t102 = t86 * t92 - t88 * t90;
t98 = qJD(3) * t102 - qJD(2);
t12 = qJD(1) * t98 + t78;
t148 = qJD(3) * t90;
t135 = t75 * t148;
t145 = qJD(5) * t92;
t27 = t135 + (qJ(5) * t148 - t145) * qJD(1);
t1 = qJD(6) * t5 + t12 * t89 + t27 * t91;
t2 = -qJD(6) * t6 + t12 * t91 - t27 * t89;
t125 = -t1 * t91 + t2 * t89;
t137 = qJD(3) * qJD(6);
t144 = qJD(6) * t90;
t33 = -t91 * t137 + (-t144 * t89 + t147 * t92) * qJD(1);
t34 = t89 * t137 + (-t142 * t92 - t144 * t91) * qJD(1);
t197 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t33 + Ifges(7,6) * t34;
t150 = qJD(3) * mrSges(6,1);
t159 = mrSges(7,1) * t59 - mrSges(7,2) * t60 - mrSges(6,3) * t152 - t150;
t74 = -mrSges(5,2) * t152 + qJD(3) * mrSges(5,3);
t136 = t74 - t159;
t196 = m(5) * t65 - t136;
t36 = -mrSges(7,2) * t77 + mrSges(7,3) * t59;
t37 = mrSges(7,1) * t77 - mrSges(7,3) * t60;
t110 = t89 * t36 + t91 * t37;
t138 = qJD(1) * qJD(3);
t134 = t90 * t138;
t19 = -mrSges(7,1) * t134 - mrSges(7,3) * t33;
t20 = mrSges(7,2) * t134 + mrSges(7,3) * t34;
t112 = t89 * t19 - t91 * t20;
t193 = t110 * qJD(6) + t112;
t192 = t33 / 0.2e1;
t191 = t34 / 0.2e1;
t190 = -t59 / 0.2e1;
t188 = -t60 / 0.2e1;
t186 = -t77 / 0.2e1;
t172 = Ifges(7,5) * t89;
t169 = Ifges(7,6) * t91;
t146 = qJD(5) * t90;
t23 = qJD(1) * t146 + (qJD(4) + t45) * qJD(3);
t153 = qJ(5) + t94;
t66 = t153 * t90;
t168 = t23 * t66;
t167 = t23 * t90;
t164 = t75 * t92;
t162 = t89 * t92;
t161 = t91 * t92;
t160 = -mrSges(5,2) + mrSges(6,3);
t149 = qJD(3) * mrSges(4,2);
t70 = -mrSges(4,3) * t152 - t149;
t158 = t70 + t74;
t157 = (mrSges(5,2) + mrSges(4,3)) * t151 + t204;
t156 = qJ(2) * mrSges(4,1);
t155 = qJ(2) * mrSges(4,2);
t154 = qJ(4) * t90;
t139 = qJD(1) * qJD(2);
t48 = -qJD(3) * pkin(3) + qJD(4) - t164;
t132 = -t48 + t164;
t131 = qJD(3) * t153;
t129 = -0.3e1 / 0.2e1 * Ifges(6,4) - 0.3e1 / 0.2e1 * Ifges(4,4) + 0.3e1 / 0.2e1 * Ifges(5,5);
t128 = -Ifges(6,5) / 0.2e1 + t203 - Ifges(4,6) / 0.2e1;
t127 = -Ifges(6,6) / 0.2e1 - Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t124 = -t1 * t89 - t2 * t91;
t120 = mrSges(7,1) * t91 - mrSges(7,2) * t89;
t117 = -Ifges(7,1) * t89 - t173;
t115 = -Ifges(7,2) * t91 - t174;
t113 = pkin(3) * t92 + t154;
t111 = t91 * t36 - t89 * t37;
t85 = t92 * qJ(4);
t38 = t85 + t103;
t67 = t153 * t92;
t13 = t38 * t91 + t67 * t89;
t14 = t38 * t89 - t67 * t91;
t109 = t48 * t90 + t56 * t92;
t71 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t151;
t108 = t111 + t71;
t105 = t92 * t93 - t154;
t101 = qJD(3) * t113 + qJD(2);
t100 = qJD(3) * t105 - qJD(2);
t99 = -m(7) * t123 - t110;
t30 = qJD(3) * t93 - t198;
t82 = Ifges(6,4) * t152;
t97 = t30 * mrSges(6,3) + t49 * mrSges(5,3) + t6 * mrSges(7,2) - t77 * Ifges(7,3) - t60 * Ifges(7,5) - t59 * Ifges(7,6) - Ifges(6,2) * t151 / 0.2e1 + t82 / 0.2e1 - t39 * mrSges(6,1) - t5 * mrSges(7,1) + ((Ifges(4,1) + Ifges(5,1)) * t92 + (-Ifges(4,4) + Ifges(5,5)) * t90) * t202 + (Ifges(6,6) + Ifges(5,4) + Ifges(4,5)) * t201;
t76 = t92 * mrSges(6,2) * t138;
t68 = t133 - t85;
t64 = (mrSges(4,1) * t90 + mrSges(4,2) * t92) * qJD(1);
t61 = t113 * qJD(1);
t57 = t85 + t121;
t47 = (qJD(4) + t164) * qJD(3);
t46 = t105 * qJD(1);
t44 = t65 + t79;
t43 = -t84 + t101;
t42 = t131 * t92 + t146;
t41 = t131 * t90 - t145;
t35 = t84 + t100;
t31 = qJD(1) * t101 - t78;
t25 = t102 * qJD(1);
t22 = qJD(1) * t100 + t78;
t18 = t84 + t98;
t11 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t10 = t25 * t89 + t44 * t91;
t9 = t25 * t91 - t44 * t89;
t8 = t33 * Ifges(7,1) + t34 * Ifges(7,4) - Ifges(7,5) * t134;
t7 = t33 * Ifges(7,4) + t34 * Ifges(7,2) - Ifges(7,6) * t134;
t4 = -qJD(6) * t14 + t18 * t91 - t41 * t89;
t3 = qJD(6) * t13 + t18 * t89 + t41 * t91;
t15 = [t66 * t11 + t13 * t19 + t14 * t20 + t3 * t36 + t35 * t62 + t4 * t37 + t41 * t71 + t43 * t63 + t57 * t76 - t159 * t42 + m(7) * (t1 * t14 + t13 * t2 + t3 * t6 + t32 * t42 + t4 * t5 + t168) + m(6) * (t22 * t57 - t27 * t67 + t30 * t41 + t35 * t39 - t40 * t42 + t168) + m(5) * (t31 * t68 + t43 * t49) + (t64 + 0.2e1 * t205) * qJD(2) + (-t31 * mrSges(5,3) + t22 * mrSges(6,1) - t27 * mrSges(6,3) + mrSges(4,2) * t139 + ((t68 * mrSges(5,1) + t66 * mrSges(6,3) + t129 * t92 + 0.2e1 * t156) * qJD(1) + t128 * qJD(3) + (m(5) * t56 + t158) * t94 + t206) * qJD(3) + t197) * t92 + (t31 * mrSges(5,1) + t22 * mrSges(6,2) + t91 * t8 / 0.2e1 + t7 * t184 + t116 * t191 + t118 * t192 + mrSges(4,1) * t139 + (m(5) * t94 - mrSges(5,2)) * t47 + (mrSges(6,3) + t119) * t23 + t124 * mrSges(7,3) + ((-t169 - t172) * t185 + t115 * t189 + t117 * t187 + t32 * t120 + t16 * t183 + t17 * t184 - t199 * mrSges(7,3)) * qJD(6) + (((-t171 / 0.2e1 + t170 / 0.2e1 - t129) * t90 - t67 * mrSges(6,3) + t68 * mrSges(5,3) - 0.2e1 * t155 - t57 * mrSges(6,1) + (-Ifges(7,3) - 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t92) * qJD(1) + t132 * mrSges(5,2) + t127 * qJD(3) + (-m(5) * t132 + t157) * t94 + t97) * qJD(3)) * t90; t193 * t92 + m(7) * (t123 * qJD(6) * t92 - t1 * t161 + t2 * t162 + t167) + m(6) * (-t27 * t92 + t167) + (m(5) * t109 + (-m(6) * t40 + m(7) * t32 - t196 + t70) * t92) * qJD(3) + (t11 + m(5) * t47 + (m(6) * t30 + m(7) * t199 + t108 + t157) * qJD(3)) * t90 + (-m(5) * t49 + m(6) * t39 - t200 - t205 - t64 - t99) * qJD(1); (t114 * t186 + t116 * t190 + t118 * t188 + t86 * t99 - t207) * qJD(6) + (((t172 / 0.2e1 + t169 / 0.2e1 + t93 * mrSges(6,3) + pkin(3) * mrSges(5,2) + t127) * qJD(3) + (t155 + (-Ifges(4,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t90) * qJD(1) - t82 / 0.2e1 + t48 * mrSges(5,2) - t97) * t90 + ((qJ(4) * t160 + t128) * qJD(3) + (-t156 + (Ifges(6,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t92) * qJD(1) - t81 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t152 - t206) * t92) * qJD(1) + ((-t149 - t158) * t92 + (t204 - t157) * t90) * t75 + t159 * t45 + t136 * qJD(4) + t88 * t11 - t44 * t71 - t46 * t62 - t61 * t63 + (mrSges(6,1) + t120) * t23 + t125 * mrSges(7,3) - t112 * t86 + t47 * mrSges(5,3) - t10 * t36 - t9 * t37 + t27 * mrSges(6,2) + t7 * t183 + t8 * t184 + t115 * t191 + t117 * t192 + (-t10 * t6 - t125 * t86 - t198 * t32 + t23 * t88 - t5 * t9) * m(7) + (qJ(4) * t23 + t198 * t40 + t27 * t93 - t30 * t44 - t39 * t46) * m(6) + (-pkin(3) * t135 + qJ(4) * t47 + qJD(4) * t56 - t109 * t75 - t49 * t61) * m(5); t196 * qJD(3) + (t160 * t148 + (-t110 + t200) * t92) * qJD(1) - m(5) * (qJD(3) * t56 - t151 * t49) + (-qJD(3) * t32 - t77 * t123 - t125) * m(7) + (qJD(3) * t40 - t151 * t39 + t27) * m(6) - t193; t91 * t19 + t89 * t20 + t76 + t111 * qJD(6) + m(6) * t22 + m(7) * (qJD(6) * t199 - t124) + (t108 * t92 + (-t150 + t159) * t90 - m(6) * (-t30 * t92 - t40 * t90) - m(7) * (-t161 * t6 + t162 * t5 + t32 * t90)) * qJD(1); -Ifges(7,3) * t134 - t32 * (mrSges(7,1) * t60 + mrSges(7,2) * t59) + (Ifges(7,1) * t59 - t175) * t188 + t16 * t187 + (Ifges(7,5) * t59 - Ifges(7,6) * t60) * t186 - t5 * t36 + t6 * t37 + (t5 * t59 + t6 * t60) * mrSges(7,3) + (-Ifges(7,2) * t60 + t17 + t58) * t190 + t197;];
tauc  = t15(:);
