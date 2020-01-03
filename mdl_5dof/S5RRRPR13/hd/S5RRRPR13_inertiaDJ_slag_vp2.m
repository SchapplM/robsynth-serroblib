% Calculate time derivative of joint inertia matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:20
% EndTime: 2019-12-31 21:43:29
% DurationCPUTime: 3.45s
% Computational Cost: add. (2899->440), mult. (7783->636), div. (0->0), fcn. (6776->8), ass. (0->184)
t226 = Ifges(5,1) + Ifges(4,3);
t225 = Ifges(4,5) - Ifges(5,4);
t224 = -Ifges(4,6) + Ifges(5,5);
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t150 = cos(qJ(3));
t183 = qJD(5) * t150;
t147 = sin(qJ(3));
t185 = qJD(3) * t147;
t156 = t146 * t183 + t149 * t185;
t155 = t146 * t185 - t149 * t183;
t210 = -t146 / 0.2e1;
t223 = t149 / 0.2e1;
t144 = sin(pkin(5));
t148 = sin(qJ(2));
t195 = t144 * t148;
t135 = pkin(7) * t195;
t145 = cos(pkin(5));
t151 = cos(qJ(2));
t207 = pkin(1) * t151;
t97 = t145 * t207 - t135;
t184 = qJD(3) * t150;
t194 = t144 * t151;
t98 = t145 * t148 * pkin(1) + pkin(7) * t194;
t86 = pkin(8) * t145 + t98;
t87 = (-pkin(2) * t151 - pkin(8) * t148 - pkin(1)) * t144;
t187 = qJD(2) * t144;
t91 = (pkin(2) * t148 - pkin(8) * t151) * t187;
t92 = t97 * qJD(2);
t158 = -t147 * t91 - t150 * t92 - t87 * t184 + t86 * t185;
t186 = qJD(2) * t148;
t12 = -(qJ(4) * t186 - qJD(4) * t151) * t144 + t158;
t203 = mrSges(6,3) * t150;
t111 = mrSges(6,1) * t147 + t146 * t203;
t112 = -mrSges(6,2) * t147 - t149 * t203;
t196 = qJ(4) * t147;
t215 = pkin(3) + pkin(9);
t100 = -t215 * t150 - pkin(2) - t196;
t214 = pkin(4) + pkin(8);
t126 = t214 * t147;
t58 = -t100 * t146 + t126 * t149;
t59 = t100 * t149 + t126 * t146;
t222 = m(6) * (-t146 * t58 + t149 * t59) + t149 * t112 - t146 * t111;
t221 = 0.2e1 * m(6);
t220 = -2 * mrSges(3,3);
t219 = 2 * mrSges(3,3);
t178 = t144 * t186;
t177 = t151 * t187;
t96 = t145 * t147 + t150 * t195;
t65 = t96 * qJD(3) + t147 * t177;
t179 = t147 * t195;
t95 = -t145 * t150 + t179;
t67 = t146 * t194 + t149 * t95;
t28 = t67 * qJD(5) + t146 * t65 + t149 * t178;
t218 = t28 / 0.2e1;
t217 = t67 / 0.2e1;
t159 = -t146 * t95 + t149 * t194;
t216 = -t159 / 0.2e1;
t166 = -Ifges(6,5) * t146 - Ifges(6,6) * t149;
t213 = t166 * qJD(5) / 0.2e1;
t212 = Ifges(6,5) * t223 + Ifges(6,6) * t210;
t199 = Ifges(6,4) * t149;
t122 = -Ifges(6,2) * t146 + t199;
t211 = t122 / 0.2e1;
t209 = -t149 / 0.2e1;
t93 = t98 * qJD(2);
t208 = mrSges(4,1) * t93;
t206 = t92 * mrSges(3,2);
t205 = t93 * mrSges(3,1);
t204 = t93 * mrSges(4,2);
t44 = t147 * t87 + t150 * t86;
t202 = Ifges(4,4) * t147;
t201 = Ifges(4,4) * t150;
t200 = Ifges(6,4) * t146;
t198 = Ifges(5,6) * t147;
t197 = Ifges(5,6) * t150;
t124 = Ifges(6,1) * t149 - t200;
t192 = t146 * t124;
t191 = t146 * t215;
t189 = t149 * t122;
t188 = t149 * t215;
t182 = 0.2e1 * t144;
t27 = t159 * qJD(5) - t146 * t178 + t149 * t65;
t66 = -qJD(3) * t179 + (qJD(3) * t145 + t177) * t150;
t3 = Ifges(6,5) * t28 + Ifges(6,6) * t27 + Ifges(6,3) * t66;
t180 = m(5) * pkin(8) + mrSges(5,1);
t127 = t214 * t150;
t172 = t194 / 0.2e1;
t50 = t66 * mrSges(5,1) + mrSges(5,2) * t178;
t43 = -t147 * t86 + t150 * t87;
t171 = pkin(3) * t185 - qJD(4) * t147;
t40 = pkin(3) * t194 - t43;
t22 = pkin(4) * t96 + pkin(9) * t194 + t40;
t85 = t135 + (-pkin(2) - t207) * t145;
t154 = -qJ(4) * t96 + t85;
t29 = t215 * t95 + t154;
t6 = -t146 * t29 + t149 * t22;
t7 = t146 * t22 + t149 * t29;
t170 = t146 * t6 - t149 * t7;
t169 = mrSges(6,1) * t149 - mrSges(6,2) * t146;
t118 = mrSges(6,1) * t146 + mrSges(6,2) * t149;
t117 = t150 * mrSges(5,2) - t147 * mrSges(5,3);
t168 = Ifges(6,1) * t146 + t199;
t167 = Ifges(6,2) * t149 + t200;
t165 = -pkin(3) * t150 - t196;
t114 = qJD(3) * t127;
t81 = (pkin(9) * t147 - qJ(4) * t150) * qJD(3) + t171;
t25 = t58 * qJD(5) + t114 * t146 + t149 * t81;
t26 = -t59 * qJD(5) + t114 * t149 - t146 * t81;
t164 = -t146 * t25 - t149 * t26;
t41 = -mrSges(6,2) * t96 + mrSges(6,3) * t67;
t42 = mrSges(6,1) * t96 + mrSges(6,3) * t159;
t163 = -t146 * t42 + t149 * t41;
t18 = -t147 * t92 + t150 * t91 - t86 * t184 - t87 * t185;
t39 = qJ(4) * t194 - t44;
t161 = t226 * t178 + t224 * t65 + t225 * t66;
t20 = -Ifges(6,4) * t159 + Ifges(6,2) * t67 + Ifges(6,6) * t96;
t21 = -Ifges(6,1) * t159 + Ifges(6,4) * t67 + Ifges(6,5) * t96;
t160 = t20 * t209 + t21 * t210;
t53 = t155 * Ifges(6,5) + t156 * Ifges(6,6) + Ifges(6,3) * t184;
t153 = -qJ(4) * t66 - qJD(4) * t96 + t93;
t142 = Ifges(4,5) * t184;
t141 = Ifges(5,5) * t185;
t130 = Ifges(3,5) * t177;
t125 = Ifges(4,1) * t147 + t201;
t123 = Ifges(4,2) * t150 + t202;
t120 = -Ifges(5,2) * t147 - t197;
t119 = -Ifges(5,3) * t150 - t198;
t115 = -pkin(2) + t165;
t113 = t214 * t185;
t110 = (Ifges(4,1) * t150 - t202) * qJD(3);
t109 = t168 * qJD(5);
t108 = (-Ifges(4,2) * t147 + t201) * qJD(3);
t107 = t167 * qJD(5);
t105 = (-Ifges(5,2) * t150 + t198) * qJD(3);
t104 = (Ifges(5,3) * t147 - t197) * qJD(3);
t103 = (mrSges(4,1) * t147 + mrSges(4,2) * t150) * qJD(3);
t102 = (-mrSges(5,2) * t147 - mrSges(5,3) * t150) * qJD(3);
t101 = t169 * qJD(5);
t99 = t169 * t150;
t94 = -qJ(4) * t184 + t171;
t90 = Ifges(6,5) * t147 - t168 * t150;
t89 = Ifges(6,6) * t147 - t167 * t150;
t88 = Ifges(6,3) * t147 + t166 * t150;
t77 = mrSges(6,1) * t184 - t155 * mrSges(6,3);
t76 = -mrSges(6,2) * t184 + t156 * mrSges(6,3);
t72 = -mrSges(4,1) * t194 - mrSges(4,3) * t96;
t71 = mrSges(4,2) * t194 - mrSges(4,3) * t95;
t70 = mrSges(5,1) * t96 - mrSges(5,2) * t194;
t69 = mrSges(5,1) * t95 + mrSges(5,3) * t194;
t57 = -t156 * mrSges(6,1) + t155 * mrSges(6,2);
t56 = -mrSges(5,2) * t95 - mrSges(5,3) * t96;
t55 = -t124 * t183 + (Ifges(6,5) * t150 + t168 * t147) * qJD(3);
t54 = -t122 * t183 + (Ifges(6,6) * t150 + t167 * t147) * qJD(3);
t52 = mrSges(4,1) * t178 - mrSges(4,3) * t66;
t51 = -mrSges(4,2) * t178 - mrSges(4,3) * t65;
t49 = mrSges(5,1) * t65 - mrSges(5,3) * t178;
t48 = Ifges(4,1) * t96 - Ifges(4,4) * t95 - Ifges(4,5) * t194;
t47 = Ifges(4,4) * t96 - Ifges(4,2) * t95 - Ifges(4,6) * t194;
t46 = -Ifges(5,4) * t194 - Ifges(5,2) * t96 + Ifges(5,6) * t95;
t45 = -Ifges(5,5) * t194 - Ifges(5,6) * t96 + Ifges(5,3) * t95;
t38 = pkin(3) * t95 + t154;
t37 = -mrSges(6,1) * t67 - mrSges(6,2) * t159;
t36 = -mrSges(5,2) * t65 - mrSges(5,3) * t66;
t35 = mrSges(4,1) * t65 + mrSges(4,2) * t66;
t34 = Ifges(4,1) * t66 - Ifges(4,4) * t65 + Ifges(4,5) * t178;
t33 = Ifges(4,4) * t66 - Ifges(4,2) * t65 + Ifges(4,6) * t178;
t32 = Ifges(5,4) * t178 - Ifges(5,2) * t66 + Ifges(5,6) * t65;
t31 = Ifges(5,5) * t178 - Ifges(5,6) * t66 + Ifges(5,3) * t65;
t30 = -pkin(4) * t95 - t39;
t19 = -Ifges(6,5) * t159 + Ifges(6,6) * t67 + Ifges(6,3) * t96;
t16 = pkin(3) * t65 + t153;
t15 = -pkin(3) * t178 - t18;
t14 = mrSges(6,1) * t66 - mrSges(6,3) * t28;
t13 = -mrSges(6,2) * t66 + mrSges(6,3) * t27;
t11 = t215 * t65 + t153;
t10 = -pkin(4) * t65 - t12;
t9 = pkin(4) * t66 - t215 * t178 - t18;
t8 = -mrSges(6,1) * t27 + mrSges(6,2) * t28;
t5 = Ifges(6,1) * t28 + Ifges(6,4) * t27 + Ifges(6,5) * t66;
t4 = Ifges(6,4) * t28 + Ifges(6,2) * t27 + Ifges(6,6) * t66;
t2 = -t7 * qJD(5) - t11 * t146 + t149 * t9;
t1 = t6 * qJD(5) + t11 * t149 + t146 * t9;
t17 = [(t31 - t33 + 0.2e1 * t208) * t95 + (t3 - t32 + t34 + 0.2e1 * t204) * t96 + (t130 - 0.2e1 * t205 - 0.2e1 * t206) * t145 + 0.2e1 * m(4) * (-t158 * t44 + t18 * t43 + t85 * t93) - 0.2e1 * t158 * t71 - t159 * t5 + 0.2e1 * m(5) * (t12 * t39 + t15 * t40 + t16 * t38) + 0.2e1 * m(3) * (t92 * t98 - t93 * t97) + 0.2e1 * t7 * t13 + 0.2e1 * t6 * t14 + (t19 + t48 - t46) * t66 + (t45 - t47) * t65 + (t1 * t7 + t10 * t30 + t2 * t6) * t221 + (t93 * t148 * t219 + (t92 * t219 - t161) * t151 + ((t97 * t220 + Ifges(3,5) * t145 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t151) * t182) * t151 + (t98 * t220 - 0.2e1 * Ifges(3,6) * t145 + t225 * t96 + t224 * t95 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t148) * t182 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t226) * t194) * t148) * qJD(2)) * t144 + t27 * t20 + t28 * t21 + 0.2e1 * t30 * t8 + 0.2e1 * t10 * t37 + 0.2e1 * t38 * t36 + 0.2e1 * t1 * t41 + 0.2e1 * t2 * t42 + 0.2e1 * t39 * t49 + 0.2e1 * t40 * t50 + 0.2e1 * t44 * t51 + 0.2e1 * t43 * t52 + 0.2e1 * t16 * t56 + t67 * t4 + 0.2e1 * t12 * t69 + 0.2e1 * t15 * t70 + 0.2e1 * t18 * t72 + 0.2e1 * t85 * t35; -t205 - t206 + (t15 * mrSges(5,1) - t18 * mrSges(4,3) + t3 / 0.2e1 - t32 / 0.2e1 + t34 / 0.2e1 + t204) * t147 + ((-t142 / 0.2e1 - t141 / 0.2e1) * t151 + (-Ifges(3,6) + (-Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1) * t150 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t147) * t186) * t144 + ((-t46 / 0.2e1 + t48 / 0.2e1 + t19 / 0.2e1 + Ifges(5,4) * t172 + t40 * mrSges(5,1) - t43 * mrSges(4,3)) * t150 + (t39 * mrSges(5,1) - t44 * mrSges(4,3) + t45 / 0.2e1 - t47 / 0.2e1 + Ifges(4,6) * t172 - t160) * t147) * qJD(3) + ((-t49 + t51) * t150 + (t50 - t52) * t147 + ((t70 - t72) * t150 + (t69 - t71) * t147) * qJD(3) + m(4) * (-t147 * t18 - t150 * t158 - t43 * t184 - t44 * t185) + m(5) * (-t12 * t150 + t147 * t15 + t40 * t184 + t39 * t185)) * pkin(8) + (-t12 * mrSges(5,1) - t158 * mrSges(4,3) + t5 * t210 + t4 * t209 - t31 / 0.2e1 + t33 / 0.2e1 - t208 + (t146 * t20 / 0.2e1 + t21 * t209) * qJD(5)) * t150 + m(6) * (t1 * t59 + t10 * t127 - t113 * t30 + t2 * t58 + t25 * t7 + t26 * t6) + m(5) * (t115 * t16 + t38 * t94) + (-m(4) * t93 - t35) * pkin(2) + (t53 / 0.2e1 - t105 / 0.2e1 + t110 / 0.2e1) * t96 + (t104 / 0.2e1 - t108 / 0.2e1) * t95 + (t88 / 0.2e1 - t120 / 0.2e1 + t125 / 0.2e1) * t66 + (t119 / 0.2e1 - t123 / 0.2e1) * t65 + t130 + t55 * t216 + t54 * t217 + t90 * t218 + t25 * t41 + t26 * t42 + t30 * t57 + t58 * t14 + t59 * t13 + t7 * t76 + t6 * t77 + t27 * t89 / 0.2e1 + t94 * t56 + t10 * t99 + t38 * t102 + t85 * t103 + t2 * t111 + t1 * t112 - t113 * t37 + t115 * t36 + t16 * t117 + t127 * t8; (-t113 * t127 + t25 * t59 + t26 * t58) * t221 + 0.2e1 * t59 * t76 + 0.2e1 * t58 * t77 - 0.2e1 * pkin(2) * t103 + 0.2e1 * t26 * t111 + 0.2e1 * t25 * t112 - 0.2e1 * t113 * t99 + 0.2e1 * t94 * t117 + 0.2e1 * t127 * t57 + 0.2e1 * (m(5) * t94 + t102) * t115 + (-t105 + t110 + t53 + (t146 * t90 + t149 * t89 + t119 - t123) * qJD(3)) * t147 + (-t146 * t55 - t149 * t54 - t104 + t108 + (t146 * t89 - t149 * t90) * qJD(5) + (-t120 + t125 + t88) * qJD(3)) * t150; (t170 * mrSges(6,3) - (-m(6) * t170 + t163) * t215 + t160) * qJD(5) + t161 - t12 * mrSges(5,3) + t15 * mrSges(5,2) + t158 * mrSges(4,2) + t18 * mrSges(4,1) + (-t2 * mrSges(6,3) - t215 * t14 + t5 / 0.2e1) * t149 + (-t1 * mrSges(6,3) - t215 * t13 - t4 / 0.2e1) * t146 + (t37 - t69) * qJD(4) + (t8 - t49) * qJ(4) + m(6) * (qJ(4) * t10 + qJD(4) * t30 - t1 * t191 - t2 * t188) + m(5) * (-pkin(3) * t15 - qJ(4) * t12 - qJD(4) * t39) - pkin(3) * t50 + t30 * t101 + t96 * t213 - t107 * t217 - t109 * t216 + t10 * t118 + t66 * t212 + t27 * t211 + t124 * t218; t141 + t142 - t76 * t191 - t77 * t188 + m(6) * (-qJ(4) * t113 + qJD(4) * t127 - t26 * t188 - t25 * t191) + qJ(4) * t57 + qJD(4) * t99 - t113 * t118 + t127 * t101 + t54 * t210 + t147 * t213 + t55 * t223 + t164 * mrSges(6,3) + (t180 * qJD(4) - t107 * t209 - t109 * t210) * t150 + ((-t89 / 0.2e1 - t59 * mrSges(6,3) - t150 * t124 / 0.2e1) * t149 + (t58 * mrSges(6,3) - t90 / 0.2e1 + t150 * t211) * t146 - t222 * t215) * qJD(5) + ((-pkin(3) * mrSges(5,1) - Ifges(5,4) + t212) * t150 + (t192 / 0.2e1 - qJ(4) * mrSges(5,1) - Ifges(4,6) + t189 / 0.2e1) * t147 + (m(5) * t165 - t150 * mrSges(4,1) + t147 * mrSges(4,2) + t117) * pkin(8)) * qJD(3); 0.2e1 * qJ(4) * t101 + t107 * t146 - t109 * t149 + (-t189 - t192) * qJD(5) + 0.2e1 * (mrSges(5,3) + t118 + (m(5) + m(6)) * qJ(4)) * qJD(4); t146 * t13 + t149 * t14 + t163 * qJD(5) + m(6) * (-t170 * qJD(5) + t1 * t146 + t149 * t2) + m(5) * t15 + t50; -m(6) * t164 + t222 * qJD(5) + t146 * t76 + t149 * t77 + t180 * t184; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t26 - mrSges(6,2) * t25 + t53; ((mrSges(6,2) * t215 - Ifges(6,6)) * t149 + (mrSges(6,1) * t215 - Ifges(6,5)) * t146) * qJD(5); -t118 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
