% Calculate time derivative of joint inertia matrix for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:03
% EndTime: 2019-03-09 06:57:12
% DurationCPUTime: 4.44s
% Computational Cost: add. (7200->429), mult. (15366->625), div. (0->0), fcn. (14291->10), ass. (0->184)
t172 = sin(qJ(4));
t173 = sin(qJ(3));
t176 = cos(qJ(4));
t177 = cos(qJ(3));
t144 = t172 * t173 - t176 * t177;
t267 = qJD(3) + qJD(4);
t115 = t267 * t144;
t171 = sin(qJ(5));
t175 = cos(qJ(5));
t219 = t171 ^ 2 + t175 ^ 2;
t197 = t219 * t115;
t269 = t219 * t176;
t146 = t172 * t177 + t176 * t173;
t218 = qJD(5) * t171;
t205 = t146 * t218;
t228 = t115 * t175;
t183 = t205 + t228;
t217 = qJD(5) * t175;
t229 = t115 * t171;
t184 = t146 * t217 - t229;
t38 = t184 * mrSges(6,1) - t183 * mrSges(6,2);
t160 = sin(pkin(11)) * pkin(1) + pkin(7);
t247 = pkin(8) + t160;
t200 = qJD(3) * t247;
t134 = t173 * t200;
t196 = t177 * t200;
t137 = t247 * t173;
t138 = t247 * t177;
t96 = -t137 * t172 + t138 * t176;
t55 = qJD(4) * t96 - t134 * t172 + t176 * t196;
t279 = m(6) * t55 + t38;
t270 = -t176 * t137 - t138 * t172;
t54 = qJD(4) * t270 - t176 * t134 - t172 * t196;
t116 = t267 * t146;
t251 = pkin(4) * t116;
t63 = pkin(3) * qJD(3) * t173 + pkin(9) * t115 + t251;
t207 = -cos(pkin(11)) * pkin(1) - pkin(2);
t153 = -pkin(3) * t177 + t207;
t92 = t144 * pkin(4) - t146 * pkin(9) + t153;
t14 = t171 * t63 + t175 * t54 + t92 * t217 - t218 * t96;
t201 = -t171 * t54 + t175 * t63;
t91 = t175 * t96;
t52 = t171 * t92 + t91;
t15 = -qJD(5) * t52 + t201;
t278 = t14 * t175 - t15 * t171;
t170 = sin(qJ(6));
t174 = cos(qJ(6));
t187 = t170 * t171 - t174 * t175;
t266 = qJD(5) + qJD(6);
t113 = t266 * t187;
t145 = t170 * t175 + t171 * t174;
t114 = t266 * t145;
t220 = -Ifges(7,5) * t113 - Ifges(7,6) * t114;
t277 = Ifges(6,5) * t217 + t220;
t51 = -t171 * t96 + t175 * t92;
t276 = (-t171 * t52 - t175 * t51) * qJD(5) + t278;
t27 = -t114 * t146 + t115 * t187;
t94 = t187 * t146;
t28 = t145 * t115 + t266 * t94;
t10 = -t28 * mrSges(7,1) + t27 * mrSges(7,2);
t274 = t10 + t38;
t162 = pkin(3) * t172 + pkin(9);
t246 = -pkin(10) - t162;
t199 = qJD(5) * t246;
t242 = pkin(3) * qJD(4);
t210 = t176 * t242;
t128 = t171 * t199 + t175 * t210;
t129 = -t171 * t210 + t175 * t199;
t141 = t246 * t171;
t166 = t175 * pkin(10);
t142 = t162 * t175 + t166;
t97 = t141 * t174 - t142 * t170;
t43 = qJD(6) * t97 + t128 * t174 + t129 * t170;
t98 = t141 * t170 + t142 * t174;
t44 = -qJD(6) * t98 - t128 * t170 + t129 * t174;
t273 = t44 * mrSges(7,1) - t43 * mrSges(7,2);
t257 = -pkin(10) - pkin(9);
t158 = t257 * t171;
t159 = pkin(9) * t175 + t166;
t125 = t158 * t174 - t159 * t170;
t206 = qJD(5) * t257;
t150 = t171 * t206;
t151 = t175 * t206;
t73 = qJD(6) * t125 + t150 * t174 + t151 * t170;
t126 = t158 * t170 + t159 * t174;
t74 = -qJD(6) * t126 - t150 * t170 + t151 * t174;
t272 = t74 * mrSges(7,1) - t73 * mrSges(7,2);
t271 = -Ifges(6,5) * t228 + Ifges(6,3) * t116;
t268 = -t171 * t51 + t175 * t52;
t265 = 0.2e1 * m(6);
t264 = 2 * m(7);
t263 = 0.2e1 * pkin(3);
t262 = 0.2e1 * t55;
t65 = mrSges(7,1) * t114 - mrSges(7,2) * t113;
t261 = 0.2e1 * t65;
t119 = mrSges(7,1) * t187 + mrSges(7,2) * t145;
t260 = 0.2e1 * t119;
t259 = 0.2e1 * t153;
t258 = m(5) / 0.2e1;
t252 = pkin(3) * t176;
t249 = t55 * t270;
t245 = Ifges(6,4) * t171;
t244 = Ifges(6,4) * t175;
t243 = Ifges(6,6) * t171;
t241 = pkin(5) * qJD(6);
t240 = t114 * mrSges(7,3);
t238 = t187 * mrSges(7,3);
t235 = t172 * mrSges(5,1);
t232 = t176 * mrSges(5,2);
t86 = t144 * t116;
t226 = t144 * t172;
t225 = t146 * t171;
t224 = t146 * t175;
t155 = -mrSges(6,1) * t175 + mrSges(6,2) * t171;
t221 = t172 * t155;
t216 = qJD(6) * t170;
t215 = qJD(6) * t174;
t214 = 0.2e1 * mrSges(7,3);
t213 = 0.2e1 * t177;
t212 = mrSges(7,3) * t241;
t211 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t116;
t209 = pkin(5) * t218;
t208 = t174 * t113 * mrSges(7,3);
t164 = -pkin(5) * t175 - pkin(4);
t203 = -t218 / 0.2e1;
t202 = -(2 * Ifges(5,4)) - t243;
t198 = t116 * mrSges(5,1) - t115 * mrSges(5,2);
t195 = mrSges(6,3) * t269;
t194 = mrSges(6,1) * t171 + mrSges(6,2) * t175;
t193 = Ifges(6,1) * t175 - t245;
t192 = -Ifges(6,2) * t171 + t244;
t191 = Ifges(6,5) * t171 + Ifges(6,6) * t175;
t190 = -t116 * t270 + t144 * t55;
t37 = pkin(5) * t144 - pkin(10) * t224 + t51;
t41 = -pkin(10) * t225 + t52;
t16 = -t170 * t41 + t174 * t37;
t17 = t170 * t37 + t174 * t41;
t48 = mrSges(6,1) * t116 + mrSges(6,3) * t183;
t49 = -mrSges(6,2) * t116 - mrSges(6,3) * t184;
t189 = -t171 * t48 + t175 * t49;
t104 = -mrSges(6,2) * t144 - mrSges(6,3) * t225;
t105 = mrSges(6,1) * t144 - mrSges(6,3) * t224;
t188 = t175 * t104 - t171 * t105;
t11 = -pkin(10) * t184 + t14;
t9 = pkin(10) * t228 + pkin(5) * t116 + (-t91 + (pkin(10) * t146 - t92) * t171) * qJD(5) + t201;
t3 = qJD(6) * t16 + t11 * t174 + t170 * t9;
t4 = -qJD(6) * t17 - t11 * t170 + t174 * t9;
t186 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t211;
t185 = -t174 * t187 * t212 + (-pkin(5) * t240 + t145 * t212) * t170 + t277;
t120 = Ifges(7,4) * t145 - Ifges(7,2) * t187;
t121 = Ifges(7,1) * t145 - Ifges(7,4) * t187;
t148 = t192 * qJD(5);
t149 = t193 * qJD(5);
t156 = Ifges(6,2) * t175 + t245;
t157 = Ifges(6,1) * t171 + t244;
t66 = -Ifges(7,4) * t113 - Ifges(7,2) * t114;
t67 = -Ifges(7,1) * t113 - Ifges(7,4) * t114;
t181 = -t113 * t121 - t114 * t120 + t145 * t67 + t175 * t148 + t171 * t149 - t156 * t218 + t157 * t217 - t187 * t66;
t147 = t194 * qJD(5);
t93 = t145 * t146;
t180 = -t27 * t238 + t94 * t240 + (-t113 * t93 - t145 * t28) * mrSges(7,3) - t198 + (t147 + t65) * t144 + (t119 + t155) * t116 - t197 * mrSges(6,3);
t179 = m(6) * (-t217 * t51 - t218 * t52 + t278) - t104 * t218 - t105 * t217 + t189;
t31 = -Ifges(6,4) * t183 - Ifges(6,2) * t184 + Ifges(6,6) * t116;
t32 = -Ifges(6,1) * t183 - Ifges(6,4) * t184 + Ifges(6,5) * t116;
t33 = pkin(5) * t184 + t55;
t45 = -Ifges(7,4) * t94 - Ifges(7,2) * t93 + Ifges(7,6) * t144;
t46 = -Ifges(7,1) * t94 - Ifges(7,4) * t93 + Ifges(7,5) * t144;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t116;
t79 = pkin(5) * t225 - t270;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t116;
t80 = Ifges(6,6) * t144 + t146 * t192;
t81 = Ifges(6,5) * t144 + t146 * t193;
t178 = (t146 * t203 - t228 / 0.2e1) * t157 + (t16 * t113 - t4 * t145) * mrSges(7,3) + (Ifges(7,5) * t145 - Ifges(7,6) * t187 + t191) * t116 / 0.2e1 - t187 * t7 / 0.2e1 - t184 * t156 / 0.2e1 + t276 * mrSges(6,3) + (-Ifges(6,6) * t218 + t277) * t144 / 0.2e1 - t3 * t238 - t17 * t240 + t80 * t203 - t270 * t147 + (-mrSges(5,1) + t155) * t55 + t81 * t217 / 0.2e1 + t149 * t224 / 0.2e1 - t148 * t225 / 0.2e1 - t54 * mrSges(5,2) + t79 * t65 - t93 * t66 / 0.2e1 - t94 * t67 / 0.2e1 - t113 * t46 / 0.2e1 - t114 * t45 / 0.2e1 - Ifges(5,5) * t115 - Ifges(5,6) * t116 + t33 * t119 + t28 * t120 / 0.2e1 + t27 * t121 / 0.2e1 + t145 * t8 / 0.2e1 + t171 * t32 / 0.2e1 + t175 * t31 / 0.2e1;
t163 = -pkin(4) - t252;
t154 = t164 - t252;
t152 = t172 * t242 + t209;
t140 = (-mrSges(7,1) * t170 - mrSges(7,2) * t174) * t241;
t100 = t194 * t146;
t78 = mrSges(7,1) * t144 + mrSges(7,3) * t94;
t77 = -mrSges(7,2) * t144 - mrSges(7,3) * t93;
t56 = mrSges(7,1) * t93 - mrSges(7,2) * t94;
t19 = -mrSges(7,2) * t116 + mrSges(7,3) * t28;
t18 = mrSges(7,1) * t116 - mrSges(7,3) * t27;
t1 = [(t14 * t52 + t15 * t51 - t249) * t265 + 0.2e1 * m(5) * (t54 * t96 - t249) + (-0.2e1 * t54 * mrSges(5,3) - t202 * t115 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t116 + t211 + t271) * t144 + (mrSges(5,3) * t262 - 0.2e1 * Ifges(5,1) * t115 - t171 * t31 + t175 * t32 + (Ifges(6,5) * t175 + t202) * t116 + (-t144 * t191 - t171 * t81 - t175 * t80) * qJD(5)) * t146 + t116 * (-Ifges(7,5) * t94 - Ifges(7,6) * t93) + t198 * t259 + t100 * t262 + (t16 * t4 + t17 * t3 + t33 * t79) * t264 - t81 * t228 + 0.2e1 * (t115 * t270 - t116 * t96) * mrSges(5,3) - 0.2e1 * t270 * t38 + t80 * t229 + ((t207 * mrSges(4,2) + Ifges(4,4) * t177) * t213 + (m(5) * pkin(3) * t259 + (mrSges(5,1) * t144 + mrSges(5,2) * t146) * t263 + 0.2e1 * t207 * mrSges(4,1) - 0.2e1 * Ifges(4,4) * t173 + (Ifges(4,1) - Ifges(4,2)) * t213) * t173) * qJD(3) + 0.2e1 * t16 * t18 + 0.2e1 * t17 * t19 + t28 * t45 + t27 * t46 + 0.2e1 * t51 * t48 + 0.2e1 * t52 * t49 + 0.2e1 * t33 * t56 + 0.2e1 * t3 * t77 + 0.2e1 * t4 * t78 + 0.2e1 * t79 * t10 - t93 * t7 - t94 * t8 + 0.2e1 * t14 * t104 + 0.2e1 * t15 * t105; -t93 * t18 - t94 * t19 + t27 * t77 + t28 * t78 + t274 * t144 + (t100 + t56) * t116 - t188 * t115 + ((-t171 * t104 - t175 * t105) * qJD(5) + t189) * t146 + m(7) * (t116 * t79 + t144 * t33 + t16 * t28 + t17 * t27 - t3 * t94 - t4 * t93) + m(6) * (-t268 * t115 + t146 * t276 + t190) + m(5) * (-t115 * t96 + t146 * t54 + t190); 0.2e1 * m(7) * (-t27 * t94 - t28 * t93 + t86) + 0.2e1 * m(6) * (-t146 * t197 + t86) + 0.2e1 * m(5) * (-t115 * t146 + t86); t178 + t179 * t162 + (m(5) * (t172 * t54 - t176 * t55) + (t176 * t115 - t172 * t116) * mrSges(5,3) + ((m(5) * t96 + m(6) * t268 - t144 * mrSges(5,3) + t188) * t176 + (t146 * mrSges(5,3) + t100 - (m(6) + m(5)) * t270) * t172) * qJD(4)) * pkin(3) + m(7) * (t152 * t79 + t154 * t33 + t16 * t44 + t17 * t43 + t3 * t98 + t4 * t97) + (Ifges(4,5) * t177 - Ifges(4,6) * t173 + (-mrSges(4,1) * t177 + mrSges(4,2) * t173) * t160) * qJD(3) + t43 * t77 + t44 * t78 + t97 * t18 + t98 * t19 + t152 * t56 + t154 * t10 + t279 * t163; t180 + ((-t115 * t172 - t116 * t176) * t258 + (m(6) * (t269 * t146 + t226) / 0.2e1 + (t146 * t176 + t226) * t258) * qJD(4)) * t263 + m(6) * (t116 * t163 - t162 * t197) + m(7) * (t116 * t154 + t144 * t152 + t27 * t98 + t28 * t97 - t43 * t94 - t44 * t93) + (-mrSges(4,1) * t173 - mrSges(4,2) * t177) * qJD(3); t152 * t260 + t154 * t261 + (t152 * t154 + t43 * t98 + t44 * t97) * t264 + 0.2e1 * t163 * t147 + (t113 * t97 - t114 * t98 - t145 * t44 - t187 * t43) * t214 + (-0.2e1 * t232 - 0.2e1 * t235 + 0.2e1 * t221 + (t162 * t269 + t163 * t172) * t265 + 0.2e1 * t195) * t242 + t181; t178 + t179 * pkin(9) + m(7) * (t125 * t4 + t126 * t3 + t16 * t74 + t164 * t33 + t17 * t73 + t209 * t79) + t56 * t209 + t73 * t77 + t74 * t78 + t125 * t18 + t126 * t19 + t164 * t10 - t279 * pkin(4); t180 + m(6) * (-pkin(9) * t197 - t251) + m(7) * (t116 * t164 + t125 * t28 + t126 * t27 + t144 * t209 - t73 * t94 - t74 * t93); m(7) * (t125 * t44 + t126 * t43 + t152 * t164 + t154 * t209 + t73 * t98 + t74 * t97) + (t154 + t164) * t65 + (-pkin(4) + t163) * t147 + (t152 + t209) * t119 + (m(6) * (-pkin(4) * t172 + pkin(9) * t269) - t235 + t221 - t232 + t195) * t242 + ((-t44 - t74) * t145 - (t43 + t73) * t187 - (t126 + t98) * t114 - (-t125 - t97) * t113) * mrSges(7,3) + t181; t209 * t260 + t164 * t261 - 0.2e1 * pkin(4) * t147 + (t125 * t74 + t126 * t73 + t164 * t209) * t264 + (t113 * t125 - t114 * t126 - t145 * t74 - t187 * t73) * t214 + t181; -Ifges(6,5) * t205 + t15 * mrSges(6,1) - t14 * mrSges(6,2) - t184 * Ifges(6,6) + (-t78 * t216 + t174 * t18 + m(7) * (-t16 * t216 + t17 * t215 + t170 * t3 + t174 * t4) + t77 * t215 + t170 * t19) * pkin(5) + t186 + t271; m(7) * (t170 * t27 + t174 * t28 + (t170 * t93 - t174 * t94) * qJD(6)) * pkin(5) - t274; -t194 * t210 + (t155 * t162 - t243) * qJD(5) + (m(7) * (t170 * t43 + t174 * t44 + t215 * t98 - t216 * t97) + t208) * pkin(5) + t185 + t273; (pkin(9) * t155 - t243) * qJD(5) + (m(7) * (-t125 * t216 + t126 * t215 + t170 * t73 + t174 * t74) + t208) * pkin(5) + t185 + t272; 0.2e1 * t140; t186; -t10; t220 + t273; t220 + t272; t140; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
