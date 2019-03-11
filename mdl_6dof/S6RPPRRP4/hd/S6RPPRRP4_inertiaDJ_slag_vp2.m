% Calculate time derivative of joint inertia matrix for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:12
% EndTime: 2019-03-09 02:05:18
% DurationCPUTime: 2.80s
% Computational Cost: add. (1591->330), mult. (3312->474), div. (0->0), fcn. (2420->6), ass. (0->148)
t96 = cos(qJ(5));
t144 = qJD(5) * t96;
t95 = sin(qJ(4));
t126 = t95 * t144;
t97 = cos(qJ(4));
t147 = qJD(4) * t97;
t94 = sin(qJ(5));
t104 = t94 * t147 + t126;
t145 = qJD(5) * t95;
t127 = t94 * t145;
t128 = t96 * t147;
t197 = -t128 + t127;
t188 = m(7) / 0.2e1 + m(6) / 0.2e1;
t198 = 0.2e1 * t188;
t152 = t94 ^ 2 + t96 ^ 2;
t89 = t95 ^ 2;
t91 = t97 ^ 2;
t196 = qJD(4) * (t89 - t91);
t148 = qJD(4) * t95;
t93 = cos(pkin(9));
t150 = qJD(2) * t93;
t132 = t97 * t150;
t143 = qJD(5) * t97;
t92 = sin(pkin(9));
t98 = -pkin(1) - pkin(2);
t122 = -t92 * qJ(2) + t93 * t98;
t44 = pkin(3) - t122;
t30 = t97 * pkin(4) + t95 * pkin(8) + t44;
t151 = qJD(2) * t92;
t39 = t151 + (-pkin(4) * t95 + pkin(8) * t97) * qJD(4);
t153 = t93 * qJ(2) + t92 * t98;
t45 = -pkin(7) + t153;
t4 = -(qJD(5) * t30 - t148 * t45 + t132) * t94 - (t143 * t45 - t39) * t96;
t2 = pkin(5) * t148 - t4;
t195 = m(7) * t2;
t194 = m(7) + m(6);
t193 = mrSges(7,2) + mrSges(6,3);
t192 = Ifges(5,1) - Ifges(5,2);
t191 = -Ifges(7,2) - Ifges(6,3);
t187 = m(7) * qJ(6) + mrSges(7,3);
t165 = t95 * t96;
t54 = mrSges(6,1) * t97 + mrSges(6,3) * t165;
t55 = -mrSges(7,1) * t97 - mrSges(7,2) * t165;
t157 = -t54 + t55;
t27 = mrSges(6,2) * t148 + mrSges(6,3) * t104;
t28 = mrSges(7,2) * t104 - mrSges(7,3) * t148;
t160 = t27 + t28;
t186 = t157 * qJD(5) + t160;
t164 = t96 * t97;
t10 = t45 * t164 + t94 * t30;
t167 = t94 * t95;
t53 = -mrSges(6,2) * t97 + mrSges(6,3) * t167;
t56 = mrSges(7,2) * t167 + mrSges(7,3) * t97;
t158 = t53 + t56;
t25 = -mrSges(6,1) * t148 - mrSges(6,3) * t197;
t26 = mrSges(7,1) * t148 + mrSges(7,2) * t197;
t161 = -t25 + t26;
t185 = t161 + m(6) * (-t10 * qJD(5) - t4) - t158 * qJD(5);
t184 = 0.2e1 * t45;
t180 = pkin(5) * t94;
t179 = mrSges(5,2) * t95;
t176 = Ifges(6,4) * t94;
t175 = Ifges(6,4) * t96;
t174 = Ifges(7,5) * t94;
t173 = Ifges(7,5) * t96;
t172 = Ifges(7,6) * t94;
t130 = t96 * t148;
t166 = t94 * t97;
t40 = t166 * t92 + t93 * t96;
t20 = -qJD(5) * t40 - t130 * t92;
t171 = t20 * t96;
t140 = t92 * t164;
t146 = qJD(5) * t94;
t21 = -t148 * t92 * t94 + qJD(5) * t140 - t146 * t93;
t170 = t21 * t94;
t169 = t30 * t96;
t168 = t92 * t95;
t163 = t97 * Ifges(6,6);
t16 = -mrSges(7,1) * t104 - mrSges(7,3) * t197;
t17 = -mrSges(6,1) * t104 + mrSges(6,2) * t197;
t162 = t16 + t17;
t114 = mrSges(7,1) * t94 - mrSges(7,3) * t96;
t46 = t114 * qJD(5);
t115 = mrSges(6,1) * t94 + mrSges(6,2) * t96;
t47 = t115 * qJD(5);
t159 = t46 + t47;
t63 = -t96 * mrSges(6,1) + t94 * mrSges(6,2);
t156 = t63 - mrSges(5,1);
t142 = qJD(6) * t96;
t155 = qJ(6) * t128 + t95 * t142;
t154 = t152 * pkin(8) * t147;
t62 = -t96 * mrSges(7,1) - t94 * mrSges(7,3);
t137 = t62 + t156;
t136 = t96 * t132 + t30 * t144 + t94 * t39;
t135 = Ifges(6,5) * t127 + Ifges(6,6) * t104;
t133 = t92 * t150;
t131 = t92 * t147;
t129 = t95 * t147;
t64 = -Ifges(7,3) * t96 + t174;
t65 = Ifges(6,2) * t96 + t176;
t124 = -t64 / 0.2e1 + t65 / 0.2e1;
t66 = Ifges(7,1) * t94 - t173;
t67 = Ifges(6,1) * t94 + t175;
t123 = t66 / 0.2e1 + t67 / 0.2e1;
t42 = t114 * t95;
t43 = t115 * t95;
t121 = (-t42 - t43) * qJD(4);
t120 = t147 * t184;
t118 = -m(6) * pkin(4) + t156;
t113 = Ifges(6,1) * t96 - t176;
t112 = Ifges(7,1) * t96 + t174;
t111 = -Ifges(6,2) * t94 + t175;
t110 = -Ifges(7,4) * t96 - t172;
t109 = Ifges(7,3) * t94 + t173;
t108 = -pkin(5) * t96 - qJ(6) * t94;
t107 = -qJ(6) * t96 + t180;
t106 = Ifges(7,4) * t127 - Ifges(7,6) * t126;
t105 = t108 * qJD(5);
t85 = Ifges(7,4) * t144;
t84 = Ifges(6,5) * t144;
t83 = Ifges(7,6) * t146;
t58 = -pkin(4) + t108;
t57 = t89 * t133;
t52 = t113 * qJD(5);
t51 = t112 * qJD(5);
t50 = t111 * qJD(5);
t49 = t109 * qJD(5);
t48 = (-t95 * mrSges(5,1) - t97 * mrSges(5,2)) * qJD(4);
t41 = -t93 * t94 + t140;
t38 = qJD(5) * t107 - qJD(6) * t94;
t36 = Ifges(6,5) * t97 - t113 * t95;
t35 = Ifges(7,4) * t97 - t112 * t95;
t34 = -t111 * t95 + t163;
t33 = Ifges(7,6) * t97 - t109 * t95;
t31 = t89 * t45 * t150;
t19 = (-t107 + t45) * t95;
t18 = pkin(8) * t171;
t14 = t67 * t145 + (-Ifges(6,5) * t95 - t113 * t97) * qJD(4);
t13 = t66 * t145 + (-Ifges(7,4) * t95 - t112 * t97) * qJD(4);
t12 = t65 * t145 + (-Ifges(6,6) * t95 - t111 * t97) * qJD(4);
t11 = t64 * t145 + (-Ifges(7,6) * t95 - t109 * t97) * qJD(4);
t9 = -t166 * t45 + t169;
t7 = -t169 + (t45 * t94 - pkin(5)) * t97;
t6 = qJ(6) * t97 + t10;
t5 = (t45 - t180) * t147 + (t105 + t150) * t95 + t155;
t3 = (-t143 * t94 - t130) * t45 + t136;
t1 = (-t146 * t45 + qJD(6)) * t97 + (-t45 * t96 - qJ(6)) * t148 + t136;
t8 = [-t197 * (-t36 - t35) + t104 * (-t33 + t34) + 0.2e1 * m(7) * (t1 * t6 + t19 * t5 + t2 * t7) + (-t11 + t12) * t167 + (-t13 - t14) * t165 + t97 * ((-Ifges(6,5) * t164 - Ifges(6,3) * t95) * qJD(4) + t135) - 0.2e1 * t5 * t42 + 0.2e1 * t44 * t48 + 0.2e1 * t3 * t53 + 0.2e1 * t4 * t54 + 0.2e1 * t2 * t55 + 0.2e1 * t1 * t56 + 0.2e1 * (mrSges(5,1) * t97 + mrSges(4,1) - t179) * t151 + 0.2e1 * (m(4) * (-t122 * t92 + t153 * t93) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2) + 0.2e1 * (-t43 * t95 + mrSges(4,2)) * t150 + 0.2e1 * m(5) * (t31 + (t45 * t91 * t93 + t44 * t92) * qJD(2)) + 0.2e1 * m(6) * (t129 * t45 ^ 2 + t10 * t3 + t4 * t9 + t31) + 0.2e1 * t19 * t16 + 0.2e1 * t9 * t25 + 0.2e1 * t7 * t26 + 0.2e1 * t10 * t27 + 0.2e1 * t6 * t28 + ((t191 + t192) * t97 + (Ifges(6,5) * t96 - Ifges(6,6) * t94 - (2 * Ifges(5,4)) - t110) * t95) * t148 + t97 * ((-Ifges(7,2) * t95 + t110 * t97) * qJD(4) + t106) + t95 * t17 * t184 + (0.2e1 * Ifges(5,4) * t97 + t192 * t95) * t147 - t43 * t120 - 0.2e1 * (t91 + t89) * mrSges(5,3) * t150; -t93 * t48 + t160 * t41 + t161 * t40 + t157 * t21 + t158 * t20 + (t97 * t121 + t162 * t95) * t92 + m(7) * (t1 * t41 + t2 * t40 + t20 * t6 + t21 * t7 + (t147 * t19 + t5 * t95) * t92) + m(6) * (t10 * t20 + t120 * t168 - t21 * t9 + t3 * t41 - t4 * t40 + t57) + m(5) * (t57 + (t91 - 0.1e1) * t133); 0.2e1 * t194 * (t92 ^ 2 * t129 + t41 * t20 + t21 * t40); m(6) * t45 * t196 + (-m(7) * t5 + (t158 * t96 + t157 * t94 + m(7) * (t6 * t96 + t7 * t94) + m(6) * (t10 * t96 - t9 * t94)) * qJD(4) - t162) * t97 + (t121 + m(7) * (qJD(4) * t19 + t144 * t7 - t146 * t6) + m(6) * (-t144 * t9 - t132) + (m(6) * t3 + m(7) * t1 + t186) * t96 + (t185 + t195) * t94) * t95; t94 * (-t145 * t41 + t147 * t40 + t21 * t95) * t198 + t194 * (t126 * t40 + t41 * t128 + t20 * t165 + t196 * t92); 0.4e1 * t188 * (-0.1e1 + t152) * t129; t58 * t16 + t5 * t62 - t38 * t42 + t19 * t46 - pkin(4) * t17 + m(7) * (t19 * t38 + t5 * t58) + (-mrSges(5,2) * t150 + t84 / 0.2e1 + t85 / 0.2e1 + t83 / 0.2e1 + (t118 * t45 - Ifges(5,5)) * qJD(4)) * t97 + (t13 / 0.2e1 + t14 / 0.2e1 + t2 * mrSges(7,2) - t4 * mrSges(6,3) + t124 * t147 + (-t10 * mrSges(6,3) - t6 * mrSges(7,2) - t163 / 0.2e1 + t33 / 0.2e1 - t34 / 0.2e1) * qJD(5) + (m(7) * (-t6 * qJD(5) + t2) + t185) * pkin(8)) * t94 + (-t11 / 0.2e1 + t12 / 0.2e1 + t1 * mrSges(7,2) + t3 * mrSges(6,3) - t123 * t147 + (t7 * mrSges(7,2) - t9 * mrSges(6,3) + t35 / 0.2e1 + t36 / 0.2e1) * qJD(5) + (m(7) * (t7 * qJD(5) + t1) + m(6) * (-t9 * qJD(5) + t3) + t186) * pkin(8)) * t96 + (t45 * t47 + (-t51 / 0.2e1 - t52 / 0.2e1) * t96 + (-t49 / 0.2e1 + t50 / 0.2e1) * t94 + (t45 * mrSges(5,2) + Ifges(5,6) + (-Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t96 + (-Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1) * t94) * qJD(4) + t118 * t150 + (t123 * t94 + t124 * t96) * qJD(5)) * t95; (t159 * t95 + (t137 * t97 + t179) * qJD(4)) * t92 + m(7) * (t131 * t58 + t168 * t38 + t18) + m(6) * (-pkin(4) * t131 + t18) + t193 * (t171 + t170 + (t40 * t96 - t41 * t94) * qJD(5)) + pkin(8) * (t144 * t40 - t146 * t41 + t170) * t198; t137 * t148 + m(6) * (-pkin(4) * t148 + t154) + m(7) * (t148 * t58 + t154) + (-m(7) * t38 + (t193 * t152 - mrSges(5,2)) * qJD(4) - t159) * t97; -0.2e1 * pkin(4) * t47 + 0.2e1 * t46 * t58 + (-t49 + t50) * t96 + (t51 + t52) * t94 + 0.2e1 * (m(7) * t58 + t62) * t38 + ((t66 + t67) * t96 + (t64 - t65) * t94) * qJD(5); t1 * mrSges(7,3) + qJD(6) * t56 + qJ(6) * t28 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) - t3 * mrSges(6,2) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - pkin(5) * t26 + (t191 * t95 + (-t172 + (-Ifges(7,4) - Ifges(6,5)) * t96) * t97) * qJD(4) + t106 + t135; m(7) * qJD(6) * t41 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t21 + (-mrSges(6,2) + t187) * t20; m(7) * (-pkin(5) * t104 - qJ(6) * t127 + t155) + t162; -Ifges(6,6) * t146 + t83 + t84 + t85 + (t105 + t142) * mrSges(7,2) + (m(7) * t142 + (m(7) * t108 + t62 + t63) * qJD(5)) * pkin(8); 0.2e1 * t187 * qJD(6); t26 + t195; m(7) * t21; t104 * m(7); (m(7) * pkin(8) + mrSges(7,2)) * t144; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
