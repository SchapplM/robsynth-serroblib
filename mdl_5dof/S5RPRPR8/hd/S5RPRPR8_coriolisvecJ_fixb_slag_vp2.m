% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR8
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:18
% EndTime: 2019-12-31 18:21:26
% DurationCPUTime: 3.07s
% Computational Cost: add. (2469->358), mult. (6194->529), div. (0->0), fcn. (3907->8), ass. (0->168)
t216 = qJD(3) / 0.2e1;
t133 = sin(qJ(5));
t135 = cos(qJ(5));
t129 = sin(pkin(9));
t131 = cos(pkin(9));
t145 = t129 * t133 - t131 * t135;
t134 = sin(qJ(3));
t136 = cos(qJ(3));
t172 = t131 * t136;
t144 = pkin(4) * t134 - pkin(7) * t172;
t149 = pkin(3) * t134 - qJ(4) * t136;
t114 = t149 * qJD(1);
t125 = sin(pkin(8)) * pkin(1) + pkin(6);
t119 = t125 * qJD(1);
t112 = t134 * t119;
t93 = qJD(2) * t136 - t112;
t45 = t131 * t114 - t129 * t93;
t31 = qJD(1) * t144 + t45;
t174 = t129 * t136;
t165 = pkin(7) * t174;
t46 = t129 * t114 + t131 * t93;
t35 = -qJD(1) * t165 + t46;
t191 = pkin(7) + qJ(4);
t115 = t191 * t129;
t116 = t191 * t131;
t65 = -t115 * t135 - t116 * t133;
t215 = -qJD(4) * t145 + qJD(5) * t65 - t133 * t31 - t135 * t35;
t111 = t129 * t135 + t131 * t133;
t66 = -t115 * t133 + t116 * t135;
t214 = -qJD(4) * t111 - qJD(5) * t66 + t133 * t35 - t135 * t31;
t171 = qJD(1) * t134;
t213 = t171 / 0.2e1;
t159 = -Ifges(4,6) * qJD(3) / 0.2e1;
t212 = Ifges(4,5) * t216;
t108 = t131 * qJD(3) - t129 * t171;
t109 = t129 * qJD(3) + t131 * t171;
t155 = t135 * t108 - t109 * t133;
t100 = t145 * qJD(5);
t141 = t144 * qJD(3);
t168 = qJD(3) * t136;
t127 = qJD(2) * t168;
t73 = t127 + (qJD(4) - t112) * qJD(3);
t98 = t149 * qJD(3) - t134 * qJD(4);
t82 = t98 * qJD(1);
t29 = -t129 * t73 + t131 * t82;
t19 = qJD(1) * t141 + t29;
t166 = qJD(1) * qJD(3);
t157 = t136 * t166;
t154 = t129 * t157;
t30 = t129 * t82 + t131 * t73;
t23 = -pkin(7) * t154 + t30;
t170 = qJD(1) * t136;
t167 = t134 * qJD(2);
t94 = t119 * t136 + t167;
t76 = qJD(3) * qJ(4) + t94;
t161 = -cos(pkin(8)) * pkin(1) - pkin(2);
t106 = -pkin(3) * t136 - t134 * qJ(4) + t161;
t81 = t106 * qJD(1);
t32 = -t129 * t76 + t131 * t81;
t18 = -pkin(4) * t170 - t109 * pkin(7) + t32;
t33 = t129 * t81 + t131 * t76;
t22 = pkin(7) * t108 + t33;
t5 = -t133 * t22 + t135 * t18;
t1 = qJD(5) * t5 + t133 * t19 + t135 * t23;
t6 = t133 * t18 + t135 * t22;
t2 = -qJD(5) * t6 - t133 * t23 + t135 * t19;
t142 = t145 * t136;
t139 = qJD(3) * t142;
t27 = -qJD(1) * t139 + t155 * qJD(5);
t143 = t111 * t136;
t140 = qJD(3) * t143;
t57 = t108 * t133 + t109 * t135;
t28 = -qJD(1) * t140 - qJD(5) * t57;
t211 = -t2 * mrSges(6,1) + t1 * mrSges(6,2) - Ifges(6,5) * t27 - Ifges(6,6) * t28;
t121 = t161 * qJD(1);
t128 = Ifges(4,4) * t170;
t185 = Ifges(5,2) * t129;
t188 = Ifges(5,4) * t131;
t150 = -t185 + t188;
t189 = Ifges(5,4) * t129;
t151 = Ifges(5,1) * t131 - t189;
t190 = mrSges(5,2) * t131;
t152 = mrSges(5,1) * t129 + t190;
t198 = t131 / 0.2e1;
t199 = -t129 / 0.2e1;
t74 = -qJD(3) * pkin(3) + qJD(4) - t93;
t210 = -(t33 * t129 + t32 * t131) * mrSges(5,3) + t121 * mrSges(4,2) + t74 * t152 + Ifges(4,1) * t213 + t128 / 0.2e1 + t212 + t108 * t150 / 0.2e1 + t109 * t151 / 0.2e1 + (Ifges(5,4) * t109 + Ifges(5,2) * t108 - Ifges(5,6) * t170) * t199 + (Ifges(5,1) * t109 + Ifges(5,4) * t108 - Ifges(5,5) * t170) * t198 - t93 * mrSges(4,3);
t209 = t27 / 0.2e1;
t208 = t28 / 0.2e1;
t207 = -t155 / 0.2e1;
t206 = t155 / 0.2e1;
t205 = -t57 / 0.2e1;
t204 = t57 / 0.2e1;
t89 = t111 * t134;
t203 = -t89 / 0.2e1;
t90 = t145 * t134;
t202 = -t90 / 0.2e1;
t124 = qJD(5) - t170;
t201 = -t124 / 0.2e1;
t200 = t124 / 0.2e1;
t85 = mrSges(5,1) * t154 + t157 * t190;
t9 = -t28 * mrSges(6,1) + t27 * mrSges(6,2);
t197 = -t85 - t9;
t196 = Ifges(6,4) * t57;
t195 = pkin(4) * t129;
t186 = Ifges(5,5) * t131;
t183 = Ifges(5,6) * t129;
t84 = t94 * qJD(3);
t181 = t136 * t84;
t80 = qJD(1) * t142;
t180 = -t100 + t80;
t101 = t111 * qJD(5);
t79 = qJD(1) * t143;
t179 = -t101 + t79;
t169 = qJD(3) * t134;
t160 = t125 * t169;
t52 = t129 * t160 + t131 * t98;
t59 = t129 * t106 + t125 * t172;
t162 = mrSges(4,3) * t171;
t178 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t108 + mrSges(5,2) * t109 + t162;
t173 = t131 * t134;
t164 = t84 * t134 * t125;
t163 = mrSges(4,3) * t170;
t158 = t134 * t166;
t156 = t125 + t195;
t148 = -t129 * t29 + t131 * t30;
t146 = -t129 * t32 + t131 * t33;
t92 = t131 * t106;
t41 = -pkin(7) * t173 + t92 + (-t125 * t129 - pkin(4)) * t136;
t44 = -pkin(7) * t129 * t134 + t59;
t12 = -t133 * t44 + t135 * t41;
t13 = t133 * t41 + t135 * t44;
t67 = t167 + (qJD(1) * t195 + t119) * t136;
t138 = t121 * mrSges(4,1) + t32 * mrSges(5,1) + t5 * mrSges(6,1) + t159 - (Ifges(4,4) * t134 + Ifges(4,2) * t136) * qJD(1) / 0.2e1 + t124 * Ifges(6,3) + t57 * Ifges(6,5) + t155 * Ifges(6,6) - Ifges(5,3) * t170 / 0.2e1 + Ifges(5,6) * t108 + Ifges(5,5) * t109 - t33 * mrSges(5,2) - t6 * mrSges(6,2) - t94 * mrSges(4,3);
t126 = -pkin(4) * t131 - pkin(3);
t123 = Ifges(6,3) * t158;
t122 = -qJD(3) * mrSges(4,2) + t163;
t97 = t156 * t134;
t96 = (mrSges(5,1) * t134 - mrSges(5,3) * t172) * t166;
t95 = (-mrSges(5,2) * t134 - mrSges(5,3) * t174) * t166;
t88 = t156 * t168;
t86 = t129 * t98;
t83 = -t119 * t169 + t127;
t78 = -mrSges(5,1) * t170 - t109 * mrSges(5,3);
t77 = mrSges(5,2) * t170 + t108 * mrSges(5,3);
t64 = (Ifges(5,5) * t134 + t136 * t151) * t166;
t63 = (Ifges(5,6) * t134 + t136 * t150) * t166;
t60 = t67 * qJD(3);
t58 = -t125 * t174 + t92;
t53 = -t131 * t160 + t86;
t51 = Ifges(6,4) * t155;
t47 = -t108 * pkin(4) + t74;
t43 = t134 * t100 - t140;
t42 = -t101 * t134 - t139;
t40 = t86 + (-t125 * t173 - t165) * qJD(3);
t39 = mrSges(6,1) * t124 - mrSges(6,3) * t57;
t38 = -mrSges(6,2) * t124 + mrSges(6,3) * t155;
t34 = t141 + t52;
t21 = -mrSges(6,2) * t158 + mrSges(6,3) * t28;
t20 = mrSges(6,1) * t158 - mrSges(6,3) * t27;
t17 = -mrSges(6,1) * t155 + mrSges(6,2) * t57;
t16 = Ifges(6,1) * t57 + Ifges(6,5) * t124 + t51;
t15 = Ifges(6,2) * t155 + Ifges(6,6) * t124 + t196;
t8 = t27 * Ifges(6,1) + t28 * Ifges(6,4) + Ifges(6,5) * t158;
t7 = t27 * Ifges(6,4) + t28 * Ifges(6,2) + Ifges(6,6) * t158;
t4 = -qJD(5) * t13 - t133 * t40 + t135 * t34;
t3 = qJD(5) * t12 + t133 * t34 + t135 * t40;
t10 = [(Ifges(6,1) * t42 + Ifges(6,4) * t43) * t204 + (Ifges(6,4) * t42 + Ifges(6,2) * t43) * t206 + (Ifges(6,5) * t42 + Ifges(6,6) * t43) * t200 + t8 * t202 + t7 * t203 + (-Ifges(6,4) * t90 - Ifges(6,2) * t89) * t208 + (-Ifges(6,1) * t90 - Ifges(6,4) * t89) * t209 + (-t1 * t89 + t2 * t90 - t42 * t5 + t43 * t6) * mrSges(6,3) + t60 * (mrSges(6,1) * t89 - mrSges(6,2) * t90) + ((t138 + (-m(4) * t94 - t122) * t125 + t159) * t134 + (t212 + (-m(4) * t93 + m(5) * t74 + t178) * t125 + t210) * t136) * qJD(3) + t59 * t95 + t58 * t96 + t97 * t9 + t88 * t17 + t53 * t77 + t52 * t78 + t47 * (-mrSges(6,1) * t43 + mrSges(6,2) * t42) + t3 * t38 + t4 * t39 + t42 * t16 / 0.2e1 + t43 * t15 / 0.2e1 + t12 * t20 + t13 * t21 + m(4) * (t125 * t136 * t83 + t164) + m(5) * (t29 * t58 + t30 * t59 + t32 * t52 + t33 * t53 + t164) + m(6) * (t1 * t13 + t12 * t2 + t3 * t6 + t4 * t5 + t47 * t88 + t60 * t97) + (t64 * t198 + t63 * t199 + t125 * t85 + (mrSges(4,3) + t152) * t84 + (-t129 * t30 - t131 * t29) * mrSges(5,3)) * t134 + ((t161 * mrSges(4,1) + Ifges(6,5) * t202 + Ifges(6,6) * t203 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t186 / 0.2e1 - t183 / 0.2e1) * t134) * t134 + (t161 * mrSges(4,2) + (-0.3e1 / 0.2e1 * t186 + 0.3e1 / 0.2e1 * t183 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t136 + (-0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(6,3) / 0.2e1 + Ifges(5,1) * t131 ^ 2 / 0.2e1 + (-t188 + t185 / 0.2e1) * t129) * t134) * t136) * t166 + (t83 * mrSges(4,3) + t30 * mrSges(5,2) - t29 * mrSges(5,1) - t123 / 0.2e1 + t211) * t136; -t89 * t20 - t90 * t21 + t42 * t38 + t43 * t39 + t197 * t136 + (-t129 * t96 + t131 * t95) * t134 + ((-t129 * t78 + t131 * t77 + t122 - t163) * t136 + (t17 - t162 + t178) * t134) * qJD(3) + m(4) * (t83 * t134 - t181 + (-t134 * t93 + t136 * t94) * qJD(3)) + m(6) * (-t1 * t90 - t136 * t60 + t47 * t169 - t2 * t89 + t6 * t42 + t5 * t43) + m(5) * (-t181 + t148 * t134 + (t134 * t74 + t136 * t146) * qJD(3)); (-Ifges(6,1) * t100 - Ifges(6,4) * t101) * t204 + (-Ifges(6,4) * t100 - Ifges(6,2) * t101) * t206 + (-Ifges(6,5) * t100 - Ifges(6,6) * t101) * t200 + t215 * t38 + (t1 * t66 + t126 * t60 + t2 * t65 + t214 * t5 + t215 * t6 - t47 * t67) * m(6) + ((-t128 / 0.2e1 + (-t183 + t186) * t170 / 0.2e1 + (Ifges(4,5) / 0.2e1 + (Ifges(5,1) * t129 + t188) * t198 + (Ifges(5,2) * t131 + t189) * t199) * qJD(3) - t210) * t136 + (-t138 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t170 + Ifges(4,4) * t213 + t159 + (Ifges(5,5) * t129 + Ifges(6,5) * t111 + Ifges(5,6) * t131 - Ifges(6,6) * t145) * t216) * t134) * qJD(1) + t214 * t39 + t126 * t9 - t93 * t122 + t111 * t8 / 0.2e1 - t83 * mrSges(4,2) - t84 * mrSges(4,1) - pkin(3) * t85 - t178 * t94 + (-t179 * mrSges(6,1) + t180 * mrSges(6,2)) * t47 + t65 * t20 + t66 * t21 - t67 * t17 - t46 * t77 - t45 * t78 + (-t84 * mrSges(5,1) + t63 / 0.2e1 + qJ(4) * t95 + qJD(4) * t77 + t30 * mrSges(5,3)) * t131 + (t84 * mrSges(5,2) + t64 / 0.2e1 - qJ(4) * t96 - qJD(4) * t78 - t29 * mrSges(5,3)) * t129 + (-pkin(3) * t84 + qJ(4) * t148 + qJD(4) * t146 - t32 * t45 - t33 * t46 - t74 * t94) * m(5) + (Ifges(6,4) * t111 - Ifges(6,2) * t145) * t208 + (Ifges(6,1) * t111 - Ifges(6,4) * t145) * t209 + t60 * (mrSges(6,1) * t145 + mrSges(6,2) * t111) + (-t1 * t145 - t111 * t2 + t179 * t6 - t180 * t5) * mrSges(6,3) - t145 * t7 / 0.2e1 + (-t101 / 0.2e1 + t79 / 0.2e1) * t15 + (-t100 / 0.2e1 + t80 / 0.2e1) * t16 + (-Ifges(6,1) * t80 - Ifges(6,4) * t79) * t205 + (-Ifges(6,4) * t80 - Ifges(6,2) * t79) * t207 + (-Ifges(6,5) * t80 - Ifges(6,6) * t79) * t201; -t108 * t77 + t109 * t78 - t155 * t38 + t57 * t39 - t197 + (-t155 * t6 + t5 * t57 + t60) * m(6) + (-t108 * t33 + t109 * t32 + t84) * m(5); t123 - t47 * (mrSges(6,1) * t57 + mrSges(6,2) * t155) + (Ifges(6,1) * t155 - t196) * t205 + t15 * t204 + (Ifges(6,5) * t155 - Ifges(6,6) * t57) * t201 - t5 * t38 + t6 * t39 + (t155 * t5 + t57 * t6) * mrSges(6,3) + (-Ifges(6,2) * t57 + t16 + t51) * t207 - t211;];
tauc = t10(:);
