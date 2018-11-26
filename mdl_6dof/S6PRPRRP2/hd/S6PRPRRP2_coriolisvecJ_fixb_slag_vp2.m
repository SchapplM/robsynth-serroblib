% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2018-11-23 15:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:00:08
% EndTime: 2018-11-23 15:00:13
% DurationCPUTime: 4.98s
% Computational Cost: add. (3464->435), mult. (8894->572), div. (0->0), fcn. (6165->10), ass. (0->211)
t134 = cos(qJ(4));
t190 = qJD(2) * t134;
t124 = Ifges(5,4) * t190;
t268 = -t124 / 0.2e1;
t265 = Ifges(6,1) + Ifges(7,1);
t254 = Ifges(7,4) + Ifges(6,5);
t131 = sin(qJ(4));
t129 = cos(pkin(6));
t117 = qJD(1) * t129 + qJD(3);
t198 = t117 * t131;
t135 = cos(qJ(2));
t127 = sin(pkin(6));
t192 = qJD(1) * t127;
t174 = t135 * t192;
t112 = qJD(2) * pkin(2) + t174;
t126 = sin(pkin(11));
t128 = cos(pkin(11));
t132 = sin(qJ(2));
t175 = t132 * t192;
t70 = t126 * t112 + t128 * t175;
t68 = qJD(2) * pkin(8) + t70;
t43 = t134 * t68 + t198;
t193 = t43 * qJD(4);
t146 = t126 * t132 - t128 * t135;
t249 = qJD(2) * t127;
t91 = t146 * t249;
t79 = qJD(1) * t91;
t15 = -t131 * t79 + t193;
t267 = m(5) * (-t15 + t193);
t266 = qJD(4) / 0.2e1;
t264 = Ifges(6,6) - Ifges(7,6);
t173 = Ifges(5,5) * t266;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t184 = t133 * qJD(4);
t191 = qJD(2) * t131;
t106 = t130 * t191 - t184;
t189 = qJD(4) * t130;
t107 = t133 * t191 + t189;
t42 = t117 * t134 - t131 * t68;
t39 = -qJD(4) * pkin(4) - t42;
t13 = pkin(5) * t106 - qJ(6) * t107 + t39;
t40 = qJD(4) * pkin(9) + t43;
t145 = -pkin(4) * t134 - pkin(9) * t131 - pkin(3);
t113 = t126 * t175;
t69 = t112 * t128 - t113;
t45 = qJD(2) * t145 - t69;
t10 = -t130 * t40 + t133 * t45;
t11 = t130 * t45 + t133 * t40;
t149 = t10 * t133 + t11 * t130;
t206 = Ifges(7,5) * t133;
t153 = Ifges(7,3) * t130 + t206;
t209 = Ifges(6,4) * t133;
t157 = -Ifges(6,2) * t130 + t209;
t162 = mrSges(7,1) * t130 - mrSges(7,3) * t133;
t164 = mrSges(6,1) * t130 + mrSges(6,2) * t133;
t120 = qJD(5) - t190;
t251 = qJD(6) - t10;
t8 = -pkin(5) * t120 + t251;
t9 = qJ(6) * t120 + t11;
t166 = t9 * t130 - t8 * t133;
t229 = t133 / 0.2e1;
t231 = t130 / 0.2e1;
t232 = -t130 / 0.2e1;
t235 = t107 / 0.2e1;
t237 = t106 / 0.2e1;
t238 = -t106 / 0.2e1;
t105 = Ifges(6,4) * t106;
t208 = Ifges(7,5) * t106;
t252 = t107 * t265 + t120 * t254 - t105 + t208;
t255 = t120 / 0.2e1;
t207 = Ifges(7,5) * t130;
t210 = Ifges(6,4) * t130;
t258 = t133 * t265 + t207 - t210;
t104 = Ifges(7,5) * t107;
t47 = Ifges(7,6) * t120 + Ifges(7,3) * t106 + t104;
t211 = Ifges(6,4) * t107;
t50 = -Ifges(6,2) * t106 + Ifges(6,6) * t120 + t211;
t246 = t166 * mrSges(7,2) + t149 * mrSges(6,3) - t13 * t162 - t153 * t237 - t157 * t238 - t39 * t164 - t47 * t231 - t50 * t232 - t258 * t235 - (-t130 * t264 + t133 * t254) * t255 - t252 * t229;
t262 = t191 / 0.2e1;
t67 = -qJD(2) * pkin(3) - t69;
t263 = -t67 * mrSges(5,2) + t42 * mrSges(5,3) - Ifges(5,1) * t262 - t173 + t246 + t268;
t183 = qJD(2) * qJD(4);
t171 = t131 * t183;
t182 = qJD(4) * qJD(5);
t186 = qJD(5) * t130;
t80 = t133 * t182 + (-t131 * t186 + t134 * t184) * qJD(2);
t185 = qJD(5) * t133;
t187 = qJD(4) * t134;
t81 = t130 * t182 + (t130 * t187 + t131 * t185) * qJD(2);
t261 = (-Ifges(6,4) + Ifges(7,5)) * t81 + t265 * t80 + t254 * t171;
t259 = t42 * qJD(4);
t257 = t130 * t265 - t206 + t209;
t256 = -m(5) * t42 + m(6) * t39;
t172 = -Ifges(5,6) * qJD(4) / 0.2e1;
t250 = t130 * t254 + t133 * t264;
t14 = -t134 * t79 + t259;
t169 = pkin(4) * t131 - pkin(9) * t134;
t111 = t169 * qJD(4);
t93 = (t126 * t135 + t128 * t132) * t127;
t88 = qJD(1) * t93;
t46 = (t111 + t88) * qJD(2);
t3 = t130 * t46 + t133 * t14 + t45 * t185 - t186 * t40;
t4 = -qJD(5) * t11 - t130 * t14 + t133 * t46;
t167 = -t130 * t4 + t133 * t3;
t1 = qJ(6) * t171 + qJD(6) * t120 + t3;
t2 = -pkin(5) * t171 - t4;
t168 = t1 * t133 + t130 * t2;
t228 = pkin(2) * t128;
t103 = t145 - t228;
t122 = pkin(2) * t126 + pkin(8);
t194 = t133 * t134;
t204 = t130 * t103 + t122 * t194;
t247 = qJD(5) * t204 - t111 * t133;
t178 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t179 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t180 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t245 = t179 * t106 - t180 * t107 - t178 * t120 - t10 * mrSges(6,1) - t67 * mrSges(5,1) - t9 * mrSges(7,3) - Ifges(6,6) * t238 - Ifges(7,6) * t237 - t172 + (Ifges(5,4) * t131 + t134 * Ifges(5,2)) * qJD(2) / 0.2e1 + t11 * mrSges(6,2) + t43 * mrSges(5,3) + t8 * mrSges(7,1) - (Ifges(6,3) + Ifges(7,2)) * t255 - t254 * t235;
t244 = 0.2e1 * m(5);
t243 = t80 / 0.2e1;
t242 = -t81 / 0.2e1;
t241 = t81 / 0.2e1;
t236 = -t107 / 0.2e1;
t234 = -t120 / 0.2e1;
t230 = -t133 / 0.2e1;
t147 = t129 * t134 - t131 * t93;
t223 = t15 * t147;
t89 = qJD(2) * t93;
t78 = qJD(1) * t89;
t92 = t146 * t127;
t220 = t78 * t92;
t35 = mrSges(7,1) * t81 - mrSges(7,3) * t80;
t36 = mrSges(6,1) * t81 + mrSges(6,2) * t80;
t218 = t35 + t36;
t110 = t169 * qJD(2);
t21 = t130 * t110 + t133 * t42;
t55 = -mrSges(7,2) * t81 + mrSges(7,3) * t171;
t58 = -mrSges(6,2) * t171 - mrSges(6,3) * t81;
t217 = t55 + t58;
t56 = mrSges(6,1) * t171 - mrSges(6,3) * t80;
t57 = -mrSges(7,1) * t171 + t80 * mrSges(7,2);
t216 = -t56 + t57;
t213 = mrSges(6,3) * t106;
t83 = -mrSges(6,2) * t120 - t213;
t86 = -mrSges(7,2) * t106 + mrSges(7,3) * t120;
t215 = t83 + t86;
t212 = mrSges(6,3) * t107;
t84 = mrSges(6,1) * t120 - t212;
t85 = -mrSges(7,1) * t120 + mrSges(7,2) * t107;
t214 = t84 - t85;
t205 = t103 * t185 + t130 * t111;
t177 = mrSges(5,3) * t191;
t203 = qJD(4) * mrSges(5,1) - mrSges(6,1) * t106 - mrSges(6,2) * t107 - t177;
t200 = t103 * t133;
t196 = t122 * t130;
t195 = t130 * t134;
t188 = qJD(4) * t131;
t63 = mrSges(7,1) * t106 - mrSges(7,3) * t107;
t181 = -t63 + t203;
t176 = mrSges(5,3) * t190;
t170 = pkin(5) + t196;
t165 = mrSges(6,1) * t133 - mrSges(6,2) * t130;
t163 = mrSges(7,1) * t133 + mrSges(7,3) * t130;
t156 = Ifges(6,2) * t133 + t210;
t152 = -Ifges(7,3) * t133 + t207;
t151 = pkin(5) * t133 + qJ(6) * t130;
t150 = pkin(5) * t130 - qJ(6) * t133;
t66 = t129 * t131 + t134 * t93;
t32 = t130 * t92 + t133 * t66;
t31 = t130 * t66 - t133 * t92;
t20 = t110 * t133 - t130 * t42;
t144 = t122 + t150;
t142 = t130 * t216 + t133 * t217;
t141 = -t130 * t215 - t133 * t214;
t140 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t123 = -pkin(3) - t228;
t119 = Ifges(7,2) * t171;
t118 = Ifges(6,3) * t171;
t116 = -qJD(4) * mrSges(5,2) + t176;
t114 = -pkin(4) - t151;
t108 = (-mrSges(5,1) * t134 + mrSges(5,2) * t131) * qJD(2);
t101 = (mrSges(5,1) * t131 + mrSges(5,2) * t134) * t183;
t96 = qJD(5) * t150 - qJD(6) * t130;
t90 = t128 * t174 - t113;
t82 = t144 * t131;
t75 = Ifges(7,4) * t80;
t74 = Ifges(6,5) * t80;
t73 = Ifges(6,6) * t81;
t72 = Ifges(7,6) * t81;
t62 = pkin(5) * t107 + qJ(6) * t106;
t60 = -t122 * t195 + t200;
t54 = t134 * t170 - t200;
t53 = -qJ(6) * t134 + t204;
t37 = (qJD(5) * t151 - qJD(6) * t133) * t131 + t144 * t187;
t34 = t130 * t88 + t194 * t90;
t33 = -t133 * t88 + t195 * t90;
t30 = qJD(4) * t66 - t131 * t91;
t29 = qJD(4) * t147 - t134 * t91;
t26 = Ifges(6,4) * t80 - Ifges(6,2) * t81 + Ifges(6,6) * t171;
t25 = Ifges(7,5) * t80 + Ifges(7,6) * t171 + Ifges(7,3) * t81;
t24 = t198 + (qJD(2) * t150 + t68) * t134;
t23 = t188 * t196 - t247;
t22 = (-t131 * t184 - t134 * t186) * t122 + t205;
t19 = -t170 * t188 + t247;
t18 = -pkin(5) * t191 - t20;
t17 = qJ(6) * t191 + t21;
t16 = (-t122 * t186 - qJD(6)) * t134 + (-t122 * t133 + qJ(6)) * t188 + t205;
t7 = pkin(5) * t81 - qJ(6) * t80 - qJD(6) * t107 + t15;
t6 = qJD(5) * t32 + t130 * t29 - t133 * t89;
t5 = -qJD(5) * t31 + t130 * t89 + t133 * t29;
t12 = [t92 * t101 + t89 * t108 + t29 * t116 - t218 * t147 - t214 * t6 + t215 * t5 + t217 * t32 + t216 * t31 - t181 * t30 + m(5) * (t14 * t66 + t29 * t43 - t30 * t42 + t67 * t89 + t220 - t223) + m(7) * (t1 * t32 + t13 * t30 - t147 * t7 + t2 * t31 + t5 * t9 + t6 * t8) + m(6) * (-t10 * t6 + t11 * t5 + t3 * t32 + t30 * t39 - t31 * t4 - t223) + m(4) * (-t69 * t89 - t70 * t91 - t79 * t93 + t220) + (-t89 * mrSges(4,1) + t91 * mrSges(4,2) + (-t131 * t66 - t134 * t147) * qJD(4) * mrSges(5,3) + (-mrSges(3,1) * t132 - mrSges(3,2) * t135) * t249) * qJD(2); t204 * t58 + m(6) * (t10 * t23 + t11 * t22 + t204 * t3 + t4 * t60) + m(7) * (t1 * t53 + t13 * t37 + t16 * t9 + t19 * t8 + t2 * t54 + t7 * t82) - m(7) * (t33 * t8 + t34 * t9) - m(6) * (-t10 * t33 + t11 * t34) + t123 * t101 - t88 * t108 + (t25 * t231 + t26 * t232 + t153 * t241 + t157 * t242 + t78 * mrSges(5,2) + t7 * t162 + (mrSges(5,3) + t164) * t15 + (-t3 * t130 - t4 * t133) * mrSges(6,3) + (-t1 * t130 + t2 * t133) * mrSges(7,2) + (((-0.3e1 / 0.2e1 * Ifges(5,4) + t180 * t133 - t179 * t130) * t131 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - t178) * t134) * qJD(2) + t172 - t245) * qJD(4) + (t156 * t237 + t39 * t165 + t13 * t163 + t152 * t238 + t50 * t230 + (t10 * t130 - t11 * t133) * mrSges(6,3) + (-t8 * t130 - t9 * t133) * mrSges(7,2) + t257 * t236 + t250 * t234 + t252 * t232) * qJD(5) + t258 * t243 + (t47 * qJD(5) + t261) * t229 + (m(6) * t15 - qJD(4) * t116 - t267 + t36) * t122 + (-m(7) * t13 + t181 - t256) * t90) * t131 + t23 * t84 + t19 * t85 + t16 * t86 + t82 * t35 + t22 * t83 + t37 * t63 + (-t118 / 0.2e1 - t119 / 0.2e1 - t74 / 0.2e1 - t75 / 0.2e1 - t72 / 0.2e1 + t73 / 0.2e1 - t78 * mrSges(5,1) + t14 * mrSges(5,3) - t90 * t116 + t179 * t81 - t180 * t80 + (-t43 * t90 / 0.2e1 + t122 * t14 / 0.2e1) * t244 + ((-t203 + t256) * t122 + t173 + 0.3e1 / 0.2e1 * t124 - t263) * qJD(4) + t140) * t134 + t53 * t55 + t54 * t57 + t60 * t56 + (qJD(2) * t88 - t78) * mrSges(4,1) + t214 * t33 - t215 * t34 + (-t67 * t88 / 0.2e1 + t123 * t78 / 0.2e1) * t244 + (qJD(2) * t90 + t79) * mrSges(4,2) + (t69 * t88 - t70 * t90 + (-t126 * t79 - t128 * t78) * pkin(2)) * m(4); ((-t130 * t214 + t133 * t215 + t116 - t176) * qJD(4) + m(6) * (-t10 * t189 + t11 * t184 - t15) + m(7) * (t184 * t9 + t189 * t8 - t7) + t267 - t218) * t134 + (t141 * qJD(5) + (-t177 - t181) * qJD(4) + m(6) * (qJD(4) * t39 - t10 * t185 - t11 * t186 + t167) + m(7) * (qJD(4) * t13 + t185 * t8 - t186 * t9 + t168) + m(5) * (t14 - t259) + t142) * t131; m(7) * (pkin(9) * t168 + t114 * t7 + t13 * t96) + m(6) * (-pkin(4) * t15 + pkin(9) * t167) + ((Ifges(5,4) * t262 + t250 * t266 + t172 + t245) * t131 + (t268 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t191 + t173 + t263) * t134) * qJD(2) + t142 * pkin(9) + (-t246 + (-m(6) * t149 - m(7) * t166 + t141) * pkin(9)) * qJD(5) + (t96 - t24) * t63 + t114 * t35 - t42 * t116 + t261 * t231 - t18 * t85 - t17 * t86 - t21 * t83 - t20 * t84 + t257 * t243 + t26 * t229 - m(6) * (t10 * t20 + t11 * t21 + t39 * t43) - m(7) * (t13 * t24 + t17 * t9 + t18 * t8) - pkin(4) * t36 - t14 * mrSges(5,2) + t167 * mrSges(6,3) + t168 * mrSges(7,2) - t7 * t163 + (-mrSges(5,1) - t165) * t15 + t25 * t230 + t203 * t43 + t152 * t241 + t156 * t242; -t140 + t118 + t119 + t74 + t75 + t72 - t73 - t39 * (mrSges(6,1) * t107 - mrSges(6,2) * t106) - t13 * (mrSges(7,1) * t107 + mrSges(7,3) * t106) + qJD(6) * t86 - t62 * t63 + qJ(6) * t55 - pkin(5) * t57 + (t106 * t8 + t107 * t9) * mrSges(7,2) + (t212 + t214) * t11 + (-t213 - t215) * t10 + t50 * t235 + (Ifges(7,3) * t107 - t208) * t238 + (-t106 * t254 - t107 * t264) * t234 + (-pkin(5) * t2 + qJ(6) * t1 - t11 * t8 - t13 * t62 + t251 * t9) * m(7) + (-Ifges(6,2) * t107 - t105 + t252) * t237 + (-t106 * t265 + t104 - t211 + t47) * t236; t107 * t63 - t120 * t86 + 0.2e1 * (t2 / 0.2e1 + t13 * t235 + t9 * t234) * m(7) + t57;];
tauc  = t12(:);
