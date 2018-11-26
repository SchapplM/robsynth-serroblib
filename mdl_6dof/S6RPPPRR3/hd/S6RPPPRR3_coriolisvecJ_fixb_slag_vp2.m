% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:05
% EndTime: 2018-11-23 15:38:08
% DurationCPUTime: 3.27s
% Computational Cost: add. (3812->319), mult. (7965->460), div. (0->0), fcn. (5127->8), ass. (0->154)
t108 = sin(pkin(10));
t110 = cos(pkin(10));
t114 = sin(qJ(5));
t116 = cos(qJ(5));
t220 = -t108 * t114 + t116 * t110;
t207 = t220 * qJD(1);
t219 = qJD(6) + t207;
t192 = t219 / 0.2e1;
t113 = sin(qJ(6));
t115 = cos(qJ(6));
t87 = t108 * t116 + t110 * t114;
t83 = t87 * qJD(1);
t60 = qJD(5) * t113 - t115 * t83;
t194 = t60 / 0.2e1;
t59 = qJD(5) * t115 + t113 * t83;
t195 = t59 / 0.2e1;
t223 = Ifges(7,5) * t194 + Ifges(7,6) * t195 + Ifges(7,3) * t192 - Ifges(6,6) * qJD(5) / 0.2e1;
t222 = -t83 / 0.2e1;
t111 = cos(pkin(9));
t64 = t111 * t207;
t109 = sin(pkin(9));
t85 = t87 * qJD(5);
t65 = t109 * t85;
t221 = -t64 - t65;
t155 = t108 ^ 2 + t110 ^ 2;
t151 = qJD(2) * t111;
t95 = -qJD(4) + t151;
t91 = t95 * qJD(1);
t137 = t155 * t91;
t160 = Ifges(6,5) * qJD(5);
t218 = t160 / 0.2e1 + Ifges(6,1) * t222;
t217 = Ifges(6,4) * t83 / 0.2e1 + t223;
t216 = t207 * Ifges(6,2);
t74 = t220 * t109;
t123 = t111 * t113 - t115 * t74;
t154 = qJD(1) * t109;
t213 = qJD(6) * t123 - t113 * t221 - t115 * t154;
t53 = -t111 * t115 - t113 * t74;
t212 = qJD(6) * t53 - t113 * t154 + t115 * t221;
t63 = t111 * t83;
t206 = t220 * qJD(5);
t66 = t206 * t109;
t211 = t63 - t66;
t170 = qJD(5) * mrSges(6,1) + mrSges(7,1) * t59 - mrSges(7,2) * t60 + t83 * mrSges(6,3);
t210 = mrSges(5,3) * t155;
t209 = -m(3) * qJ(2) - mrSges(3,3);
t75 = qJD(5) * t207;
t35 = qJD(6) * t59 - t115 * t75;
t36 = -qJD(6) * t60 + t113 * t75;
t15 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t208 = -t75 * mrSges(6,3) + t15;
t191 = t115 / 0.2e1;
t186 = Ifges(7,4) * t60;
t23 = t59 * Ifges(7,2) + Ifges(7,6) * t219 + t186;
t56 = Ifges(7,4) * t59;
t24 = t60 * Ifges(7,1) + Ifges(7,5) * t219 + t56;
t205 = t24 * t191 - t113 * t23 / 0.2e1;
t102 = t110 * qJD(3);
t104 = t111 * qJ(2);
t117 = -pkin(1) - pkin(2);
t93 = qJD(1) * t117 + qJD(2);
t78 = qJD(1) * t104 + t109 * t93;
t62 = -qJD(1) * qJ(4) + t78;
t52 = t108 * qJD(3) + t110 * t62;
t204 = t110 * t52 - t108 * (-t108 * t62 + t102);
t45 = t102 + (pkin(7) * qJD(1) - t62) * t108;
t153 = qJD(1) * t110;
t46 = -pkin(7) * t153 + t52;
t18 = -t114 * t46 + t116 * t45;
t11 = qJD(5) * t18 + t220 * t91;
t148 = qJD(1) * qJD(2);
t139 = t109 * t148;
t76 = qJD(1) * t85;
t37 = -pkin(5) * t76 + pkin(8) * t75 + t139;
t19 = t114 * t45 + t116 * t46;
t17 = qJD(5) * pkin(8) + t19;
t77 = -qJ(2) * t154 + t111 * t93;
t61 = qJD(1) * pkin(3) + qJD(4) - t77;
t55 = pkin(4) * t153 + t61;
t27 = pkin(5) * t207 + pkin(8) * t83 + t55;
t5 = -t113 * t17 + t115 * t27;
t1 = qJD(6) * t5 + t11 * t115 + t113 * t37;
t6 = t113 * t27 + t115 * t17;
t2 = -qJD(6) * t6 - t11 * t113 + t115 * t37;
t134 = t1 * t115 - t113 * t2;
t128 = Ifges(7,5) * t115 - Ifges(7,6) * t113;
t168 = Ifges(7,4) * t115;
t129 = -Ifges(7,2) * t113 + t168;
t169 = Ifges(7,4) * t113;
t130 = Ifges(7,1) * t115 - t169;
t131 = mrSges(7,1) * t113 + mrSges(7,2) * t115;
t133 = t113 * t6 + t115 * t5;
t16 = -qJD(5) * pkin(5) - t18;
t202 = -t133 * mrSges(7,3) + t128 * t192 + t129 * t195 + t130 * t194 + t131 * t16 + t205;
t201 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t79 = Ifges(6,4) * t207;
t200 = t55 * mrSges(6,2) - t18 * mrSges(6,3) - t79 / 0.2e1 + t218;
t199 = t5 * mrSges(7,1) + t55 * mrSges(6,1) - t19 * mrSges(6,3) + t216 / 0.2e1 - t6 * mrSges(7,2) + t217;
t198 = t35 / 0.2e1;
t197 = t36 / 0.2e1;
t193 = -t76 / 0.2e1;
t161 = t109 * t117 + t104;
t88 = -qJ(4) + t161;
t190 = pkin(7) - t88;
t188 = m(5) * t109;
t12 = qJD(5) * t19 + t87 * t91;
t67 = t190 * t108;
t68 = t190 * t110;
t33 = -t114 * t68 - t116 * t67;
t183 = t12 * t33;
t73 = t87 * t109;
t182 = t12 * t73;
t181 = t12 * t220;
t174 = t76 * mrSges(6,3);
t165 = t113 * t206;
t164 = t113 * t87;
t163 = t115 * t206;
t162 = t115 * t87;
t152 = qJD(2) * t109;
t150 = qJD(6) * t113;
t149 = qJD(6) * t115;
t147 = Ifges(7,5) * t35 + Ifges(7,6) * t36 - Ifges(7,3) * t76;
t41 = -t76 * mrSges(6,1) - t75 * mrSges(6,2);
t138 = -t151 / 0.2e1;
t136 = -t109 * qJ(2) + t111 * t117;
t135 = pkin(3) - t136;
t126 = -t77 * t109 + t78 * t111;
t25 = -mrSges(7,1) * t76 - mrSges(7,3) * t35;
t26 = mrSges(7,2) * t76 + mrSges(7,3) * t36;
t125 = -t113 * t25 + t115 * t26;
t34 = t114 * t67 - t116 * t68;
t81 = t110 * pkin(4) + t135;
t38 = pkin(5) * t220 + pkin(8) * t87 + t81;
t14 = t113 * t38 + t115 * t34;
t13 = -t113 * t34 + t115 * t38;
t39 = -mrSges(7,2) * t219 + mrSges(7,3) * t59;
t40 = mrSges(7,1) * t219 - mrSges(7,3) * t60;
t124 = -t113 * t39 - t115 * t40;
t122 = t149 * t87 + t165;
t121 = t150 * t87 - t163;
t89 = qJD(1) * (mrSges(5,1) * t110 - mrSges(5,2) * t108);
t118 = qJD(1) ^ 2;
t69 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t207;
t48 = -pkin(5) * t83 + pkin(8) * t207;
t47 = mrSges(6,1) * t207 - mrSges(6,2) * t83;
t42 = -pkin(5) * t85 + pkin(8) * t206 + t152;
t20 = -qJD(5) * t33 + t220 * t95;
t10 = t35 * Ifges(7,1) + t36 * Ifges(7,4) - t76 * Ifges(7,5);
t9 = t35 * Ifges(7,4) + t36 * Ifges(7,2) - t76 * Ifges(7,6);
t8 = t113 * t48 + t115 * t18;
t7 = -t113 * t18 + t115 * t48;
t4 = -qJD(6) * t14 - t113 * t20 + t115 * t42;
t3 = qJD(6) * t13 + t113 * t42 + t115 * t20;
t21 = [(t1 * t164 - t5 * t121 + t6 * t122 + t2 * t162) * mrSges(7,3) + (-t131 - mrSges(6,3)) * t12 * t87 + (t113 * t24 + t115 * t23) * qJD(6) * t87 / 0.2e1 + (-t76 * t87 + t75 * t220 + t207 * t206 / 0.2e1 + t85 * t222) * Ifges(6,4) + (-mrSges(6,2) * t139 + Ifges(6,1) * t75 - t128 * t193 - t129 * t197 - t130 * t198) * t87 - (-mrSges(6,1) * t139 + Ifges(6,2) * t76 + t11 * mrSges(6,3) - t147 / 0.2e1 - Ifges(7,3) * t193 - Ifges(7,6) * t197 - Ifges(7,5) * t198 + t201) * t220 + (-t216 / 0.2e1 - t199 - t223) * t85 - 0.2e1 * mrSges(5,3) * t137 + m(6) * (t11 * t34 + t19 * t20 + t183) - (t200 + t205 + t218) * t206 + m(5) * (t88 * t137 + t204 * t95) + (Ifges(7,1) * t121 + Ifges(7,4) * t122) * t194 + m(4) * ((-t109 * t136 + t111 * t161) * qJD(1) + t126) * qJD(2) + (-m(6) * t18 + m(7) * t16 - t170) * (qJD(5) * t34 + t87 * t95) + t208 * t33 + 0.2e1 * (mrSges(4,2) * t111 - t209) * t148 + (Ifges(7,5) * t121 + Ifges(7,6) * t122) * t192 + m(7) * (t1 * t14 + t13 * t2 + t3 * t6 + t4 * t5 + t183) + (m(6) * (qJD(1) * t81 + t55) + m(5) * (qJD(1) * t135 + t61) + 0.2e1 * t89 + t47) * t152 + (Ifges(7,4) * t121 + Ifges(7,2) * t122) * t195 + t81 * t41 + t20 * t69 + t3 * t39 + t4 * t40 + t13 * t25 + t14 * t26 + t9 * t164 / 0.2e1 - t10 * t162 / 0.2e1 + 0.2e1 * mrSges(4,1) * t139 + t16 * (-mrSges(7,1) * t122 + mrSges(7,2) * t121) + t34 * t174; -t111 * t41 + t73 * t15 + t53 * t25 - t123 * t26 + t221 * t69 + t213 * t40 + t212 * t39 + (-t47 - t89) * t154 + (-t73 * t75 + t74 * t76) * mrSges(6,3) + (-mrSges(4,1) * t109 + (-mrSges(4,2) + t210) * t111 + t209) * t118 + m(6) * (t11 * t74 - t18 * t66 - t19 * t65 + t182) - m(6) * (-t18 * t63 + t19 * t64) + t137 * t188 + 0.2e1 * (-m(5) * (t61 * t109 + t111 * t204) / 0.2e1 + t138 * t188 - m(4) * t126 / 0.2e1 + (t138 - t55 / 0.2e1) * t109 * m(6)) * qJD(1) + t170 * t211 + (-t1 * t123 - t16 * t211 + t2 * t53 + t212 * t6 + t213 * t5 + t182) * m(7); -t208 * t220 - t170 * t85 - (t113 * t40 - t115 * t39 - t69) * t206 + m(6) * (-t18 * t85 + t19 * t206 - t181) + m(7) * (t16 * t85 + t6 * t163 - t5 * t165 - t181) + (t174 + t124 * qJD(6) + m(6) * t11 + m(7) * (-t149 * t5 - t150 * t6 + t134) + t125) * t87; t207 * t69 - t170 * t83 + (m(5) + m(6)) * t139 - t118 * t210 + (t219 * t39 + t25) * t115 + (-t219 * t40 + t26) * t113 - m(6) * (t18 * t83 - t19 * t207) + m(5) * t204 * qJD(1) + t41 + (t1 * t113 + t115 * t2 + t16 * t83 + t219 * (-t113 * t5 + t115 * t6)) * m(7); (Ifges(7,5) * t113 + Ifges(7,6) * t115) * t193 + (Ifges(7,2) * t115 + t169) * t197 + (Ifges(7,1) * t113 + t168) * t198 + t9 * t191 + t113 * t10 / 0.2e1 - t18 * t69 - Ifges(6,5) * t75 + Ifges(6,6) * t76 - t8 * t39 - t7 * t40 - t11 * mrSges(6,2) - pkin(5) * t15 + t170 * t19 + (-mrSges(7,1) * t115 + mrSges(7,2) * t113 - mrSges(6,1)) * t12 + t134 * mrSges(7,3) + (t199 + t217) * t83 - (t79 / 0.2e1 - t160 / 0.2e1 + (-Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t83 - t200 - t202) * t207 + t202 * qJD(6) + (-pkin(5) * t12 - t16 * t19 - t5 * t7 - t6 * t8) * m(7) + (t125 + m(7) * t134 + (-m(7) * t133 + t124) * qJD(6)) * pkin(8); -t16 * (mrSges(7,1) * t60 + mrSges(7,2) * t59) - t60 * (Ifges(7,1) * t59 - t186) / 0.2e1 + t23 * t194 - t219 * (Ifges(7,5) * t59 - Ifges(7,6) * t60) / 0.2e1 - t5 * t39 + t6 * t40 + (t5 * t59 + t6 * t60) * mrSges(7,3) + t147 - (-Ifges(7,2) * t60 + t24 + t56) * t59 / 0.2e1 - t201;];
tauc  = t21(:);
