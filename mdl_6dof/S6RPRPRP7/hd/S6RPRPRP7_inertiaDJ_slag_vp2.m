% Calculate time derivative of joint inertia matrix for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:13
% EndTime: 2019-03-09 03:21:21
% DurationCPUTime: 3.84s
% Computational Cost: add. (2390->300), mult. (4942->433), div. (0->0), fcn. (4291->6), ass. (0->141)
t216 = Ifges(6,4) + Ifges(7,4);
t217 = Ifges(6,1) + Ifges(7,1);
t215 = Ifges(6,2) + Ifges(7,2);
t207 = Ifges(6,5) + Ifges(7,5);
t144 = cos(pkin(9));
t178 = cos(qJ(3));
t92 = sin(pkin(9));
t94 = sin(qJ(3));
t63 = t144 * t94 + t178 * t92;
t220 = t207 * t63;
t95 = cos(qJ(5));
t219 = t216 * t95;
t93 = sin(qJ(5));
t218 = t216 * t93;
t154 = Ifges(6,6) + Ifges(7,6);
t214 = -t93 ^ 2 - t95 ^ 2;
t213 = -t215 * t93 + t219;
t212 = t217 * t95 - t218;
t141 = qJD(5) * t95;
t109 = t144 * t178;
t62 = t92 * t94 - t109;
t129 = t62 * t141;
t58 = t63 * qJD(3);
t159 = t93 * t58;
t100 = t129 + t159;
t142 = qJD(5) * t93;
t130 = t62 * t142;
t156 = t95 * t58;
t99 = t130 - t156;
t211 = t154 * t95 + t207 * t93;
t143 = qJD(3) * t94;
t57 = -qJD(3) * t109 + t143 * t92;
t210 = t100 * t215 - t154 * t57 + t216 * t99;
t209 = t216 * t100 - t207 * t57 + t217 * t99;
t164 = t58 * t62;
t208 = -Ifges(4,1) + Ifges(4,2);
t206 = t154 * t63 - t213 * t62;
t205 = -t212 * t62 + t220;
t83 = t94 * pkin(3) + qJ(2);
t37 = pkin(4) * t63 + pkin(8) * t62 + t83;
t96 = -pkin(1) - pkin(7);
t146 = qJ(4) - t96;
t70 = t146 * t94;
t98 = t146 * t178;
t44 = -t144 * t70 - t92 * t98;
t42 = t95 * t44;
t16 = t93 * t37 + t42;
t204 = t213 * qJD(5);
t203 = t212 * qJD(5);
t202 = t215 * t95 + t218;
t201 = t217 * t93 + t219;
t200 = t207 * t141;
t117 = qJD(3) * t178;
t199 = -mrSges(4,1) * t143 - mrSges(4,2) * t117;
t15 = t95 * t37 - t44 * t93;
t198 = -t15 * t141 - t16 * t142;
t147 = qJ(6) * t62;
t11 = t147 * t93 + t16;
t5 = pkin(5) * t63 + t147 * t95 + t15;
t197 = -t11 * t142 - t5 * t141;
t196 = t144 * t58 + t57 * t92;
t82 = -pkin(3) * t144 - pkin(4);
t71 = -t95 * pkin(5) + t82;
t72 = -mrSges(7,1) * t95 + mrSges(7,2) * t93;
t194 = m(7) * t71 + t72;
t192 = m(6) * t82 - mrSges(6,1) * t95 + mrSges(6,2) * t93 - mrSges(5,1);
t157 = t95 * mrSges(7,3);
t40 = mrSges(7,1) * t63 + t157 * t62;
t158 = t95 * mrSges(6,3);
t41 = mrSges(6,1) * t63 + t158 * t62;
t150 = t40 + t41;
t160 = t93 * mrSges(7,3);
t38 = -mrSges(7,2) * t63 + t160 * t62;
t161 = t93 * mrSges(6,3);
t39 = -mrSges(6,2) * t63 + t161 * t62;
t151 = t38 + t39;
t191 = t150 * t93 - t151 * t95;
t190 = 2 * m(5);
t189 = 0.2e1 * m(7);
t188 = -2 * mrSges(5,3);
t187 = -0.2e1 * mrSges(7,3);
t186 = 2 * qJD(2);
t185 = m(5) * pkin(3);
t184 = m(7) * pkin(5);
t177 = mrSges(6,2) * t95;
t53 = -qJD(3) * t98 - t94 * qJD(4);
t97 = -qJD(4) * t178 + t143 * t146;
t29 = -t144 * t97 + t53 * t92;
t43 = t144 * t98 - t70 * t92;
t171 = t29 * t43;
t170 = t43 * t58;
t166 = t57 * t93;
t165 = t57 * t95;
t163 = t62 * t93;
t162 = t62 * t95;
t155 = mrSges(6,2) + mrSges(7,2);
t17 = -mrSges(7,1) * t57 - mrSges(7,3) * t99;
t18 = -mrSges(6,1) * t57 - mrSges(6,3) * t99;
t153 = t17 + t18;
t19 = mrSges(7,2) * t57 + mrSges(7,3) * t100;
t20 = mrSges(6,2) * t57 + mrSges(6,3) * t100;
t152 = t19 + t20;
t149 = -Ifges(7,5) * t156 - Ifges(7,3) * t57;
t148 = -Ifges(6,5) * t156 - Ifges(6,3) * t57;
t64 = mrSges(7,1) * t142 + mrSges(7,2) * t141;
t81 = pkin(3) * t92 + pkin(8);
t145 = qJ(6) + t81;
t78 = pkin(3) * t117 + qJD(2);
t140 = 0.2e1 * t57 * mrSges(5,3);
t139 = t62 * t188;
t30 = t144 * t53 + t92 * t97;
t31 = -pkin(4) * t57 + pkin(8) * t58 + t78;
t138 = t37 * t141 + t95 * t30 + t93 * t31;
t137 = pkin(5) * t142;
t124 = mrSges(7,1) + t184;
t122 = t57 * t214;
t119 = -t57 * mrSges(5,1) - t58 * mrSges(5,2);
t118 = -t30 * t93 + t95 * t31;
t115 = qJD(5) * t145;
t114 = mrSges(6,1) + t124;
t3 = -t142 * t44 + t138;
t4 = -t16 * qJD(5) + t118;
t110 = t3 * t95 - t4 * t93;
t108 = mrSges(6,1) * t93 + t177;
t103 = t29 * t62 + t170;
t59 = t145 * t93;
t60 = t145 * t95;
t102 = -t59 * t93 - t60 * t95;
t101 = qJ(6) * t58 + qJD(6) * t62;
t13 = -t100 * mrSges(7,1) + t99 * mrSges(7,2);
t65 = t108 * qJD(5);
t46 = -qJD(6) * t93 - t115 * t95;
t45 = qJD(6) * t95 - t115 * t93;
t36 = t108 * t62;
t35 = (-mrSges(7,1) * t93 - mrSges(7,2) * t95) * t62;
t26 = -pkin(5) * t163 + t43;
t14 = -mrSges(6,1) * t100 + mrSges(6,2) * t99;
t12 = -pkin(5) * t100 + t29;
t2 = qJ(6) * t129 + (-qJD(5) * t44 + t101) * t93 + t138;
t1 = -pkin(5) * t57 + t101 * t95 + (-t42 + (-t37 - t147) * t93) * qJD(5) + t118;
t6 = [0.2e1 * t78 * (mrSges(5,1) * t63 - mrSges(5,2) * t62) + 0.2e1 * t12 * t35 + 0.2e1 * t2 * t38 + 0.2e1 * t3 * t39 + 0.2e1 * t1 * t40 + 0.2e1 * t4 * t41 + 0.2e1 * t43 * t14 + (((m(4) + m(3)) * t186) + 0.2e1 * (mrSges(4,1) * t178 - mrSges(4,2) * t94) * qJD(3)) * qJ(2) + 0.2e1 * t26 * t13 + 0.2e1 * t5 * t17 + 0.2e1 * t15 * t18 + 0.2e1 * t11 * t19 + 0.2e1 * t16 * t20 + (-0.2e1 * t36 + t139) * t29 + 0.2e1 * (-t57 * t62 + t58 * t63) * Ifges(5,4) + 0.2e1 * m(6) * (t15 * t4 + t16 * t3 + t171) + (t205 + t220) * t130 - t209 * t162 + t210 * t163 + 0.2e1 * Ifges(5,1) * t164 + ((-t154 * t93 + t207 * t95) * t62 + (-(2 * Ifges(5,2)) - Ifges(6,3) - Ifges(7,3)) * t63) * t57 + (-0.2e1 * Ifges(4,4) * t178 + t208 * t94) * t117 + (0.2e1 * Ifges(4,4) * t94 + t178 * t208) * t143 - t205 * t156 + (t30 * t63 + t170) * t188 + (t94 * mrSges(4,1) + mrSges(4,2) * t178 + mrSges(3,3)) * t186 + 0.2e1 * t83 * t119 + t63 * (Ifges(7,6) * t100 + t149) + t63 * (Ifges(6,6) * t100 + t148) + t44 * t140 + (t1 * t5 + t11 * t2 + t12 * t26) * t189 + (t30 * t44 + t78 * t83 + t171) * t190 + t100 * t206; (t13 + t14) * t62 + (t35 - t36 + t139) * t58 + t191 * t57 + m(7) * (-t11 * t165 + t12 * t62 + t166 * t5 + t26 * t58) + m(6) * (t15 * t166 - t16 * t165 + t103) + m(5) * (-t44 * t57 + t103) + (t140 + t152 * t95 - t153 * t93 + (-t150 * t95 - t151 * t93) * qJD(5) + m(7) * (-t1 * t93 + t2 * t95 + t197) + m(6) * (t110 + t198) + m(5) * t30) * t63; (-t57 * t63 + t164) * t190 + 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t122 * t63 + t164); t71 * t13 + t12 * t72 + t82 * t14 + t26 * t64 + t43 * t65 - t59 * t17 + t60 * t19 - Ifges(5,5) * t58 + t45 * t38 + t46 * t40 + m(7) * (-t1 * t59 + t11 * t45 + t12 * t71 + t137 * t26 + t2 * t60 + t46 * t5) + t196 * mrSges(5,3) * pkin(3) - t1 * t160 - t4 * t161 - Ifges(4,5) * t143 + (-t93 * t18 + m(6) * ((-t15 * t95 - t16 * t93) * qJD(5) + t110) - t39 * t142 - t41 * t141 + t95 * t20) * t81 + (t202 * t62 + t205) * t141 / 0.2e1 + (t201 * qJD(5) + t204) * t163 / 0.2e1 + t209 * t93 / 0.2e1 + t210 * t95 / 0.2e1 - t206 * t142 / 0.2e1 + (-t142 * t154 + t200) * t63 / 0.2e1 - t201 * t156 / 0.2e1 + t202 * t159 / 0.2e1 - t203 * t162 / 0.2e1 + t197 * mrSges(7,3) + t198 * mrSges(6,3) + t199 * t96 + (-t144 * t185 + t192) * t29 + (t185 * t92 - mrSges(5,2)) * t30 + t35 * t137 - Ifges(4,6) * t117 + (Ifges(5,6) - t211 / 0.2e1) * t57 + t2 * t157 + t3 * t158; m(7) * (pkin(5) * t130 + (t45 * t95 - t46 * t93 + (t59 * t95 - t60 * t93) * qJD(5)) * t63) + m(6) * t122 * t81 - t196 * t185 + (t65 + t64) * t62 + (t192 + t194) * t58 + t199 + (m(7) * t102 + mrSges(5,2) + (mrSges(7,3) + mrSges(6,3)) * t214) * t57; 0.2e1 * t82 * t65 + 0.2e1 * t71 * t64 + (t45 * t60 - t46 * t59) * t189 + (t46 * t187 + (0.2e1 * pkin(5) * t194 + t60 * t187 - t202) * qJD(5) + t203) * t93 + (0.2e1 * t45 * mrSges(7,3) + (-t187 * t59 + t201) * qJD(5) + t204) * t95; t153 * t95 + t152 * t93 - t191 * qJD(5) + m(7) * (t1 * t95 + t2 * t93 + (t11 * t95 - t5 * t93) * qJD(5)) + m(6) * (t3 * t93 + t4 * t95 + (-t15 * t93 + t16 * t95) * qJD(5)) + m(5) * t78 + t119; 0; m(7) * (-qJD(5) * t102 + t45 * t93 + t46 * t95); 0; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 + t154 * t159 + (m(7) * t1 + t17) * pkin(5) + t211 * t62 * qJD(5) + t148 + t149; (t114 * t93 + t155 * t95) * t57 + (-t114 * t95 + t155 * t93) * t63 * qJD(5); -mrSges(7,2) * t45 + t124 * t46 + ((-mrSges(6,1) * t81 - mrSges(7,3) * pkin(5)) * t95 + (mrSges(6,2) * t81 - t154) * t93) * qJD(5) + t200; (-t177 + (-mrSges(6,1) - t184) * t93) * qJD(5) - t64; 0; m(7) * t12 + t13; m(7) * t58; m(7) * t137 + t64; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
