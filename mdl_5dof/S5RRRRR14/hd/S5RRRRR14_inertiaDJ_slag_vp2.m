% Calculate time derivative of joint inertia matrix for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:24
% EndTime: 2024-09-27 18:42:26
% DurationCPUTime: 1.72s
% Computational Cost: add. (5269->299), mult. (14235->437), div. (0->0), fcn. (12802->10), ass. (0->168)
t139 = sin(pkin(5));
t147 = cos(qJ(3));
t179 = qJD(3) * t147;
t169 = t139 * t179;
t223 = Ifges(4,5) * t169;
t142 = sin(qJ(4));
t143 = sin(qJ(3));
t146 = cos(qJ(4));
t104 = (-t142 * t143 + t146 * t147) * t139;
t220 = qJD(3) + qJD(4);
t68 = t220 * t104;
t105 = (t142 * t147 + t143 * t146) * t139;
t69 = t220 * t105;
t193 = Ifges(5,5) * t68 - Ifges(5,6) * t69;
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t63 = t104 * t145 - t105 * t141;
t28 = qJD(5) * t63 - t141 * t69 + t145 * t68;
t64 = t104 * t141 + t105 * t145;
t29 = -qJD(5) * t64 - t141 * t68 - t145 * t69;
t194 = Ifges(6,5) * t28 + Ifges(6,6) * t29;
t177 = qJD(4) * t146;
t178 = qJD(4) * t142;
t140 = cos(pkin(5));
t144 = sin(qJ(2));
t189 = pkin(1) * qJD(2);
t172 = t144 * t189;
t162 = t140 * t172;
t124 = pkin(1) * t144 + pkin(8) * t139;
t166 = -pkin(9) * t139 - t124;
t148 = cos(qJ(2));
t135 = pkin(1) * t148 + pkin(2);
t168 = t140 * t179;
t181 = t147 * t148 * t189 + t135 * t168;
t53 = (qJD(3) * t166 - t162) * t143 + t181;
t185 = t140 * t143;
t122 = t135 * t185;
t184 = t140 * t147;
t151 = (-t143 * t148 - t144 * t184) * t189;
t54 = t151 + (t147 * t166 - t122) * qJD(3);
t123 = t135 * t184;
t137 = t140 * pkin(3);
t76 = t143 * t166 + t123 + t137;
t186 = t139 * t147;
t130 = pkin(9) * t186;
t96 = t147 * t124 + t122;
t80 = t130 + t96;
t18 = t142 * t54 + t146 * t53 + t76 * t177 - t178 * t80;
t67 = t69 * pkin(10);
t11 = t18 - t67;
t46 = t142 * t76 + t146 * t80;
t19 = -qJD(4) * t46 - t142 * t53 + t146 * t54;
t204 = pkin(10) * t68;
t12 = t19 - t204;
t165 = t140 * pkin(4) - pkin(10) * t105;
t45 = -t142 * t80 + t146 * t76;
t35 = t165 + t45;
t103 = t104 * pkin(10);
t36 = t103 + t46;
t15 = -t141 * t36 + t145 * t35;
t2 = qJD(5) * t15 + t11 * t145 + t12 * t141;
t16 = t141 * t35 + t145 * t36;
t3 = -qJD(5) * t16 - t11 * t141 + t12 * t145;
t222 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t132 = pkin(2) * t185;
t115 = pkin(8) * t186 + t132;
t100 = t130 + t115;
t133 = pkin(2) * t184;
t171 = t139 * (-pkin(8) - pkin(9));
t160 = t143 * t171;
t90 = t133 + t137 + t160;
t55 = -t100 * t142 + t146 * t90;
t44 = t165 + t55;
t56 = t146 * t100 + t142 * t90;
t47 = t103 + t56;
t20 = -t141 * t47 + t145 * t44;
t128 = pkin(2) * t168;
t91 = qJD(3) * t160 + t128;
t92 = (t147 * t171 - t132) * qJD(3);
t31 = -t100 * t178 + t142 * t92 + t146 * t91 + t90 * t177;
t22 = t31 - t67;
t32 = -qJD(4) * t56 - t142 * t91 + t146 * t92;
t23 = t32 - t204;
t5 = qJD(5) * t20 + t141 * t23 + t145 * t22;
t21 = t141 * t44 + t145 * t47;
t6 = -qJD(5) * t21 - t141 * t22 + t145 * t23;
t221 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t127 = t139 * t172;
t219 = t19 * mrSges(5,1) - t18 * mrSges(5,2) + t222;
t218 = t32 * mrSges(5,1) - t31 * mrSges(5,2) + t221;
t217 = 2 * m(4);
t216 = 2 * m(5);
t215 = 2 * m(6);
t10 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t214 = 0.2e1 * t10;
t37 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t213 = 0.2e1 * t37;
t41 = mrSges(5,1) * t69 + mrSges(5,2) * t68;
t212 = 0.2e1 * t41;
t58 = -mrSges(6,2) * t140 + mrSges(6,3) * t63;
t211 = 0.2e1 * t58;
t59 = mrSges(6,1) * t140 - mrSges(6,3) * t64;
t210 = 0.2e1 * t59;
t93 = -mrSges(5,2) * t140 + mrSges(5,3) * t104;
t209 = 0.2e1 * t93;
t94 = mrSges(5,1) * t140 - mrSges(5,3) * t105;
t208 = 0.2e1 * t94;
t138 = t139 ^ 2;
t106 = (mrSges(4,1) * t143 + mrSges(4,2) * t147) * t139 * qJD(3);
t207 = -0.2e1 * t106;
t187 = t139 * t143;
t118 = mrSges(4,1) * t140 - mrSges(4,3) * t187;
t206 = 0.2e1 * t118;
t119 = -mrSges(4,2) * t140 + mrSges(4,3) * t186;
t205 = 0.2e1 * t119;
t201 = pkin(3) * t147;
t200 = pkin(4) * t104;
t198 = t28 * mrSges(6,3);
t197 = t29 * mrSges(6,3);
t134 = pkin(3) * t146 + pkin(4);
t175 = qJD(5) * t145;
t176 = qJD(5) * t141;
t183 = t141 * t142;
t78 = t134 * t175 + (-t142 * t176 + (t145 * t146 - t183) * qJD(4)) * pkin(3);
t195 = t78 * mrSges(6,2);
t192 = mrSges(4,3) * t147;
t191 = Ifges(4,4) * t143;
t190 = Ifges(4,4) * t147;
t180 = qJD(3) * t143;
t170 = t139 * t180;
t126 = pkin(3) * t170;
t107 = t127 + t126;
t71 = -mrSges(5,1) * t104 + mrSges(5,2) * t105;
t188 = t107 * t71;
t182 = t142 * t145;
t174 = 0.2e1 * mrSges(5,3);
t173 = 0.2e1 * mrSges(6,3);
t79 = -t134 * t176 + (-t142 * t175 + (-t141 * t146 - t182) * qJD(4)) * pkin(3);
t74 = t79 * mrSges(6,1);
t167 = t74 - t195;
t57 = pkin(4) * t69 + t126;
t164 = t193 + t194;
t163 = t138 * t172;
t161 = (-mrSges(4,1) * t147 + mrSges(4,2) * t143) * t127;
t121 = (-pkin(2) - t201) * t139;
t110 = (-t135 - t201) * t139;
t159 = t164 + (t141 * t197 + t175 * t58) * pkin(4);
t156 = -Ifges(4,6) * t170 + t223;
t155 = -t146 * t68 * mrSges(5,3) - t178 * t94;
t154 = -t145 * t198 - t176 * t59;
t153 = (-mrSges(3,1) * t144 - mrSges(3,2) * t148) * t189;
t152 = (-mrSges(5,1) * t142 - mrSges(5,2) * t146) * qJD(4) * pkin(3);
t150 = 0.2e1 * t28 * t64 * Ifges(6,1) + 0.2e1 * t63 * Ifges(6,2) * t29 - 0.2e1 * t104 * Ifges(5,2) * t69 + 0.2e1 * t105 * t68 * Ifges(5,1) + (Ifges(4,1) * t143 + t190) * t139 * t169 + ((Ifges(4,1) * t147 - t191) * t180 + (-Ifges(4,2) * t143 + t190) * t179) * t138 + (t156 + 0.2e1 * t193 + 0.2e1 * t194 + t223) * t140 + 0.2e1 * (t28 * t63 + t29 * t64) * Ifges(6,4) + 0.2e1 * (t104 * t68 - t105 * t69) * Ifges(5,4);
t111 = -pkin(3) * t183 + t134 * t145;
t112 = pkin(3) * t182 + t134 * t141;
t149 = -t111 * t198 + t112 * t197 + t78 * t58 + t79 * t59 + t156 + t164 + (-mrSges(5,3) * t142 * t69 + t177 * t93) * pkin(3);
t116 = (-mrSges(6,1) * t141 - mrSges(6,2) * t145) * qJD(5) * pkin(4);
t114 = -pkin(8) * t187 + t133;
t109 = t115 * qJD(3);
t108 = -pkin(8) * t170 + t128;
t102 = Ifges(4,6) * t140 + (Ifges(4,2) * t147 + t191) * t139;
t95 = -t124 * t143 + t123;
t81 = t121 - t200;
t75 = t110 - t200;
t61 = -qJD(3) * t96 + t151;
t60 = (-qJD(3) * t124 - t162) * t143 + t181;
t48 = t127 + t57;
t1 = [(t15 * t3 + t16 * t2 + t48 * t75) * t215 + (t107 * t110 + t18 * t46 + t19 * t45) * t216 + (-t135 * t163 + t60 * t96 + t61 * t95) * t217 + (0.2e1 * t161 + t135 * t207 + (-t102 * t143 + 0.2e1 * (-t143 * t96 - t147 * t95) * mrSges(4,3)) * qJD(3)) * t139 + 0.2e1 * t153 + (-t15 * t28 + t16 * t29) * t173 + (-t45 * t68 - t46 * t69) * t174 + t150 + t60 * t205 + 0.2e1 * t188 + t110 * t212 + t61 * t206 + t18 * t209 + t19 * t208 + t75 * t214 + t2 * t211 + t3 * t210 + t48 * t213; (t161 + (-pkin(2) - t135) * t106 + ((-t114 - t95) * t192 + (-t102 + (m(5) * t110 + t71) * pkin(3) + (-t96 - t115) * mrSges(4,3)) * t143) * qJD(3)) * t139 + (-(t46 + t56) * t69 + (-t45 - t55) * t68) * mrSges(5,3) + m(6) * (t15 * t6 + t16 * t5 + t2 * t21 + t20 * t3 + t48 * t81 + t57 * t75) + m(4) * (-pkin(2) * t163 + t108 * t96 - t109 * t95 + t114 * t61 + t115 * t60) + (t32 + t19) * t94 + (t31 + t18) * t93 + (t6 + t3) * t59 + (t5 + t2) * t58 + t153 + t150 + (t121 + t110) * t41 + (t57 + t48) * t37 + ((t16 + t21) * t29 + (-t15 - t20) * t28) * mrSges(6,3) + m(5) * (t107 * t121 + t18 * t56 + t19 * t55 + t31 * t46 + t32 * t45) + t188 + (t108 + t60) * t119 + (-t109 + t61) * t118 + (t75 + t81) * t10; (-t20 * t28 + t21 * t29) * t173 + (-t55 * t68 - t56 * t69) * t174 + (pkin(2) * t207 + (-0.2e1 * t114 * t192 + (-0.2e1 * mrSges(4,3) * t115 + 0.2e1 * pkin(3) * t71 - t102) * t143) * qJD(3)) * t139 + (t20 * t6 + t21 * t5 + t57 * t81) * t215 + (t108 * t115 - t109 * t114) * t217 + (t121 * t126 + t31 * t56 + t32 * t55) * t216 + t150 + t108 * t205 + t121 * t212 - t109 * t206 + t31 * t209 + t32 * t208 + t81 * t214 + t57 * t213 + t5 * t211 + t6 * t210; m(6) * (t111 * t3 + t112 * t2 + t15 * t79 + t16 * t78) + (m(5) * (t142 * t18 + t146 * t19 + t177 * t46 - t178 * t45) + t155) * pkin(3) + t149 - t60 * mrSges(4,2) + t61 * mrSges(4,1) + t219; m(6) * (t111 * t6 + t112 * t5 + t20 * t79 + t21 * t78) + (m(5) * (t142 * t31 + t146 * t32 + t177 * t56 - t178 * t55) + t155) * pkin(3) + t149 - t108 * mrSges(4,2) - t109 * mrSges(4,1) + t218; 0.2e1 * t74 - 0.2e1 * t195 + (t111 * t79 + t112 * t78) * t215 + 0.2e1 * t152; (m(6) * (t141 * t2 + t145 * t3 - t15 * t176 + t16 * t175) + t154) * pkin(4) + t159 + t219; (m(6) * (t141 * t5 + t145 * t6 + t175 * t21 - t176 * t20) + t154) * pkin(4) + t159 + t218; t152 + (-mrSges(6,1) * t176 + m(6) * (-t111 * t176 + t112 * t175 + t141 * t78 + t145 * t79) - mrSges(6,2) * t175) * pkin(4) + t167; 0.2e1 * t116; t194 + t222; t194 + t221; t167; t116; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
