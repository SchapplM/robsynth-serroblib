% Calculate time derivative of joint inertia matrix for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:22
% EndTime: 2019-03-09 10:45:29
% DurationCPUTime: 3.07s
% Computational Cost: add. (4842->340), mult. (10262->495), div. (0->0), fcn. (9486->8), ass. (0->150)
t95 = sin(qJ(6));
t98 = cos(qJ(6));
t147 = t95 ^ 2 + t98 ^ 2;
t204 = m(4) * pkin(7) + mrSges(4,2);
t100 = cos(qJ(2));
t178 = pkin(7) - pkin(8);
t131 = t178 * t100;
t97 = sin(qJ(2));
t198 = t178 * t97;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t46 = t99 * t131 + t198 * t96;
t140 = qJD(6) * t98;
t111 = t100 * t96 - t97 * t99;
t45 = -t96 * t131 + t198 * t99;
t104 = qJ(5) * t111 + t45;
t144 = cos(pkin(10));
t112 = t100 * t99 + t96 * t97;
t34 = -qJ(5) * t112 + t46;
t94 = sin(pkin(10));
t18 = -t104 * t144 + t34 * t94;
t194 = qJD(2) - qJD(4);
t42 = t194 * t111;
t43 = t194 * t112;
t26 = t144 * t42 + t43 * t94;
t27 = t144 * t43 - t94 * t42;
t130 = qJD(4) * t178;
t119 = t100 * t130;
t120 = t97 * t130;
t28 = -qJD(2) * t45 - t119 * t96 + t99 * t120;
t106 = t46 * qJD(2);
t29 = -qJD(4) * t46 + t106;
t36 = -t111 * t94 + t112 * t144;
t116 = mrSges(7,1) * t95 + mrSges(7,2) * t98;
t67 = t116 * qJD(6);
t110 = -t43 * qJ(5) + qJD(5) * t111;
t16 = -t42 * qJ(5) - qJD(5) * t112 + t28;
t7 = t144 * t16 + (-t119 * t99 - t120 * t96 + t106 + t110) * t94;
t76 = Ifges(7,5) * t95 + Ifges(7,6) * t98;
t141 = qJD(6) * t95;
t83 = Ifges(7,6) * t141;
t203 = -(t76 / 0.2e1 - Ifges(6,6)) * t26 + t28 * mrSges(5,2) + t7 * mrSges(6,2) - Ifges(5,5) * t43 - Ifges(6,5) * t27 + Ifges(5,6) * t42 - t18 * t67 - t29 * mrSges(5,1) - t36 * (Ifges(7,5) * t140 - t83) / 0.2e1;
t37 = -t111 * t144 - t112 * t94;
t199 = -t100 * pkin(2) - t97 * qJ(3);
t74 = -pkin(1) + t199;
t59 = t100 * pkin(3) - t74;
t44 = pkin(4) * t112 + t59;
t17 = pkin(5) * t36 - pkin(9) * t37 + t44;
t19 = t104 * t94 + t144 * t34;
t10 = t17 * t95 + t19 * t98;
t139 = t10 * qJD(6);
t101 = -pkin(2) - pkin(3);
t143 = qJD(2) * t97;
t138 = qJD(2) * t100;
t148 = qJ(3) * t138 + t97 * qJD(3);
t47 = t101 * t143 + t148;
t30 = pkin(4) * t42 + t47;
t8 = pkin(5) * t26 - pkin(9) * t27 + t30;
t2 = -t7 * t95 + t8 * t98 - t139;
t200 = t2 + t139;
t133 = t37 * t141;
t153 = t98 * t27;
t108 = t133 - t153;
t109 = t37 * t140 + t95 * t27;
t4 = -Ifges(7,1) * t108 - Ifges(7,4) * t109 + Ifges(7,5) * t26;
t202 = t4 / 0.2e1 - t200 * mrSges(7,3);
t107 = t144 * t99 - t94 * t96;
t201 = (mrSges(5,1) * t96 + mrSges(5,2) * t99) * qJD(4) + t107 * t67;
t75 = -t98 * mrSges(7,1) + mrSges(7,2) * t95;
t151 = mrSges(6,1) - t75;
t181 = m(6) * pkin(4);
t197 = t181 * t94;
t72 = -qJ(3) * t96 + t99 * t101;
t11 = mrSges(7,1) * t26 + mrSges(7,3) * t108;
t12 = -mrSges(7,2) * t26 - mrSges(7,3) * t109;
t196 = -t95 * t11 + t98 * t12;
t125 = t144 * pkin(4);
t81 = -t125 - pkin(5);
t192 = m(7) * t81 - t144 * t181 - t151;
t177 = pkin(4) * t94;
t80 = pkin(9) + t177;
t191 = -mrSges(6,2) + t197 + (m(7) * t80 + mrSges(7,3)) * t147;
t190 = 2 * m(5);
t189 = 0.2e1 * m(6);
t188 = 0.2e1 * m(7);
t187 = -0.2e1 * pkin(1);
t6 = t16 * t94 - t144 * (t110 + t29);
t186 = 0.2e1 * t6;
t185 = 0.2e1 * t18;
t184 = 0.2e1 * t30;
t183 = 0.2e1 * t74;
t176 = t18 * t6;
t172 = Ifges(7,4) * t95;
t171 = Ifges(7,4) * t98;
t170 = Ifges(7,5) * t98;
t49 = t99 * qJD(3) + qJD(4) * t72;
t73 = t99 * qJ(3) + t101 * t96;
t50 = -t96 * qJD(3) - qJD(4) * t73;
t31 = -t144 * t50 + t49 * t94;
t169 = t18 * t31;
t168 = t26 * mrSges(6,3);
t167 = t31 * t107;
t166 = t37 * t95;
t71 = -pkin(4) + t72;
t40 = t144 * t71 - t94 * t73;
t38 = pkin(5) - t40;
t165 = t38 * t67;
t164 = t49 * mrSges(5,2);
t163 = t50 * mrSges(5,1);
t62 = t144 * t96 + t94 * t99;
t53 = t62 * qJD(4);
t162 = t53 * t107;
t160 = t81 * t67;
t155 = t98 * mrSges(7,3);
t152 = qJD(6) / 0.2e1;
t150 = Ifges(4,5) - Ifges(3,4);
t149 = Ifges(7,5) * t153 + Ifges(7,3) * t26;
t41 = t144 * t73 + t94 * t71;
t9 = t17 * t98 - t19 * t95;
t145 = t9 * qJD(6);
t39 = -pkin(9) + t41;
t142 = qJD(6) * t39;
t136 = 0.2e1 * t97;
t129 = t147 * mrSges(7,3);
t128 = t147 * t39;
t127 = t147 * t62;
t123 = t26 * mrSges(6,1) + t27 * mrSges(6,2);
t122 = -Ifges(7,6) * t95 - (2 * Ifges(6,4));
t118 = t10 * t98 - t9 * t95;
t117 = -t107 * t6 + t18 * t53;
t21 = -mrSges(7,2) * t36 - mrSges(7,3) * t166;
t22 = mrSges(7,1) * t36 - t155 * t37;
t115 = t98 * t21 - t95 * t22;
t69 = Ifges(7,4) * t140 - Ifges(7,2) * t141;
t70 = Ifges(7,1) * t140 - Ifges(7,4) * t141;
t114 = t98 * t69 + t95 * t70;
t113 = -t100 * mrSges(4,1) - t97 * mrSges(4,3);
t1 = t7 * t98 + t8 * t95 + t145;
t105 = t1 * t98 - t2 * t95 + (-t10 * t95 - t9 * t98) * qJD(6);
t77 = Ifges(7,2) * t98 + t172;
t78 = Ifges(7,1) * t95 + t171;
t103 = (-t77 * t95 + t78 * t98) * qJD(6) + t114;
t54 = t107 * qJD(4);
t32 = t144 * t49 + t94 * t50;
t20 = t116 * t37;
t14 = Ifges(7,5) * t36 + (Ifges(7,1) * t98 - t172) * t37;
t13 = Ifges(7,6) * t36 + (-Ifges(7,2) * t95 + t171) * t37;
t5 = mrSges(7,1) * t109 - mrSges(7,2) * t108;
t3 = -Ifges(7,4) * t108 - Ifges(7,2) * t109 + Ifges(7,6) * t26;
t15 = [0.2e1 * t47 * (mrSges(5,1) * t112 - mrSges(5,2) * t111) + 0.2e1 * t42 * Ifges(5,2) * t112 - 0.2e1 * t43 * t111 * Ifges(5,1) + 0.2e1 * t59 * (mrSges(5,1) * t42 + mrSges(5,2) * t43) + 0.2e1 * t44 * t123 + t20 * t186 + 0.2e1 * t1 * t21 + 0.2e1 * t2 * t22 + t5 * t185 + 0.2e1 * t9 * t11 + 0.2e1 * t10 * t12 - 0.2e1 * t19 * t168 + (mrSges(6,3) * t185 - t95 * t13 + t98 * t14) * t27 + (t1 * t10 + t2 * t9 + t176) * t188 + (t19 * t7 + t30 * t44 + t176) * t189 + (t28 * t46 + t29 * t45 + t47 * t59) * t190 + (mrSges(6,1) * t184 - 0.2e1 * t7 * mrSges(6,3) + t122 * t27 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t26 + t149) * t36 + (mrSges(6,2) * t184 + mrSges(6,3) * t186 + 0.2e1 * Ifges(6,1) * t27 - t95 * t3 + t98 * t4 + (t122 + t170) * t26 + (-t98 * t13 - t95 * t14 - t36 * t76) * qJD(6)) * t37 + 0.2e1 * (t111 * t42 - t112 * t43) * Ifges(5,4) + 0.2e1 * (t111 * t29 - t112 * t28 - t42 * t46 - t43 * t45) * mrSges(5,3) + (m(4) * t183 + 0.2e1 * t113) * (pkin(2) * t143 - t148) + ((mrSges(3,1) * t187 + mrSges(4,1) * t183 + t136 * t150) * t97 + (mrSges(3,2) * t187 - 0.2e1 * t74 * mrSges(4,3) + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t136 - 0.2e1 * t150 * t100) * t100) * qJD(2); (-t41 * t26 - t40 * t27 + t31 * t37 - t32 * t36) * mrSges(6,3) + (t111 * t50 - t112 * t49 - t73 * t42 - t72 * t43) * mrSges(5,3) + m(6) * (t19 * t32 - t40 * t6 + t41 * t7 + t169) + m(5) * (t28 * t73 + t29 * t72 + t45 * t50 + t46 * t49) + ((-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t97 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t100 + (m(4) * t199 - t100 * mrSges(3,1) + t97 * mrSges(3,2) + t113) * pkin(7)) * qJD(2) + (-t22 * t142 - t3 / 0.2e1 - qJD(6) * t14 / 0.2e1 - t27 * t78 / 0.2e1 + m(7) * (t1 * t39 + t10 * t32 - t142 * t9) + t32 * t21 + t39 * t12 + (t77 * t152 - t70 / 0.2e1) * t37 + (-t1 + t145) * mrSges(7,3)) * t98 + (-t21 * t142 + t13 * t152 + t27 * t77 / 0.2e1 + (t78 * t152 + t69 / 0.2e1) * t37 + (-m(7) * t9 - t22) * t32 + (-m(7) * t200 - t11) * t39 - t202) * t95 + t38 * t5 + t31 * t20 + t151 * t6 + m(7) * (t38 * t6 + t169) + t204 * qJD(3) * t100 + t203; -0.2e1 * t163 + 0.2e1 * t164 - 0.2e1 * t165 + (t49 * t73 + t50 * t72) * t190 + t103 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) + (t188 * t38 - t189 * t40 + 0.2e1 * t151) * t31 + (t128 * t188 + t189 * t41 + 0.2e1 * mrSges(6,2) - 0.2e1 * t129) * t32; t53 * t20 - t107 * t5 + t115 * t54 + t204 * t138 + ((-t21 * t95 - t22 * t98) * qJD(6) + t196) * t62 + m(7) * (t105 * t62 + t118 * t54 + t117) + m(6) * (t19 * t54 + t62 * t7 + t117) + m(5) * (t28 * t96 + t29 * t99 + (-t45 * t96 + t46 * t99) * qJD(4)) + (-t107 * t27 - t26 * t62 - t36 * t54 + t37 * t53) * mrSges(6,3) + (-t96 * t42 - t99 * t43 + (-t111 * t96 - t112 * t99) * qJD(4)) * mrSges(5,3); t151 * t53 + (mrSges(6,2) - t129) * t54 + m(7) * (t127 * t32 + t128 * t54 + t38 * t53 - t167) + m(6) * (t32 * t62 - t40 * t53 + t41 * t54 - t167) + m(5) * (t49 * t96 + t50 * t99 + (-t72 * t96 + t73 * t99) * qJD(4)) + t201; 0.2e1 * m(6) * (t54 * t62 - t162) + 0.2e1 * m(7) * (t127 * t54 - t162); -t27 * mrSges(6,3) * t125 - t69 * t166 / 0.2e1 + t78 * t153 / 0.2e1 + t81 * t5 - t168 * t177 + t1 * t155 + t7 * t197 - t109 * t77 / 0.2e1 + (t37 * t70 + t3) * t98 / 0.2e1 + (t14 / 0.2e1 - t9 * mrSges(7,3)) * t140 - (t37 * t78 + t13) * t141 / 0.2e1 + t192 * t6 + (m(7) * t105 - t22 * t140 - t21 * t141 + t196) * t80 + t202 * t95 - t203; -t78 * t140 + t77 * t141 + t191 * t32 + t192 * t31 - t114 - t160 + t163 - t164 + t165; t191 * t54 + t192 * t53 - t201; t103 + 0.2e1 * t160; t98 * t11 + t95 * t12 + t115 * qJD(6) + m(7) * (qJD(6) * t118 + t1 * t95 + t2 * t98) + m(6) * t30 + t123; 0; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,5) * t133 - Ifges(7,6) * t109 + t149; t83 - t116 * t32 + (t39 * t75 - t170) * qJD(6); (t141 * t62 - t54 * t98) * mrSges(7,2) + (-t140 * t62 - t54 * t95) * mrSges(7,1); -t83 + (t75 * t80 + t170) * qJD(6); -t67; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
