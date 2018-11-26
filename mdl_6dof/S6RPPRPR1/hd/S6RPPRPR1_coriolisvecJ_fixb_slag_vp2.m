% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2018-11-23 15:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:39:37
% EndTime: 2018-11-23 15:39:41
% DurationCPUTime: 4.04s
% Computational Cost: add. (5760->413), mult. (14437->576), div. (0->0), fcn. (10662->10), ass. (0->186)
t259 = Ifges(5,1) / 0.2e1;
t163 = cos(pkin(10));
t228 = cos(qJ(4));
t189 = t228 * t163;
t150 = qJD(1) * t189;
t166 = sin(qJ(4));
t160 = sin(pkin(10));
t193 = qJD(1) * t160;
t133 = t166 * t193 - t150;
t258 = -t133 / 0.2e1;
t257 = Ifges(5,2) + Ifges(6,3);
t185 = qJD(1) * (t160 ^ 2 + t163 ^ 2);
t256 = mrSges(4,3) * t185;
t145 = t160 * t228 + t166 * t163;
t134 = t145 * qJD(1);
t159 = sin(pkin(11));
t162 = cos(pkin(11));
t115 = qJD(4) * t159 + t134 * t162;
t165 = sin(qJ(6));
t167 = cos(qJ(6));
t186 = t162 * qJD(4) - t134 * t159;
t255 = -t115 * t165 + t167 * t186;
t73 = t115 * t167 + t165 * t186;
t172 = -t166 * t160 + t189;
t137 = t172 * qJD(4);
t127 = qJD(1) * t137;
t175 = t159 * t165 - t162 * t167;
t38 = qJD(6) * t255 - t127 * t175;
t241 = t38 / 0.2e1;
t144 = t159 * t167 + t162 * t165;
t39 = -qJD(6) * t73 - t127 * t144;
t240 = t39 / 0.2e1;
t138 = t145 * qJD(4);
t128 = qJD(1) * t138;
t233 = t128 / 0.2e1;
t254 = Ifges(5,4) * t258;
t253 = t133 / 0.2e1;
t252 = t134 * t259;
t222 = pkin(8) + qJ(5);
t146 = t222 * t159;
t147 = t222 * t162;
t112 = -t146 * t165 + t147 * t167;
t226 = pkin(8) * t162;
t108 = pkin(4) * t134 + qJ(5) * t133;
t151 = sin(pkin(9)) * pkin(1) + qJ(3);
t148 = t151 * qJD(1);
t126 = t160 * qJD(2) + t163 * t148;
t216 = pkin(7) * qJD(1);
t118 = t163 * t216 + t126;
t110 = t166 * t118;
t156 = t163 * qJD(2);
t117 = t156 + (-t148 - t216) * t160;
t190 = t228 * t117;
t75 = -t110 + t190;
t42 = t162 * t108 - t159 * t75;
t21 = pkin(5) * t134 + t133 * t226 + t42;
t199 = t133 * t159;
t43 = t159 * t108 + t162 * t75;
t28 = pkin(8) * t199 + t43;
t251 = -qJD(5) * t144 - qJD(6) * t112 + t165 * t28 - t167 * t21;
t111 = -t146 * t167 - t147 * t165;
t250 = -qJD(5) * t175 + qJD(6) * t111 - t165 * t21 - t167 * t28;
t13 = -t39 * mrSges(7,1) + t38 * mrSges(7,2);
t200 = t127 * t162;
t201 = t127 * t159;
t87 = mrSges(6,1) * t201 + mrSges(6,2) * t200;
t249 = t13 + t87;
t136 = t144 * qJD(6);
t92 = t144 * t133;
t205 = -t92 - t136;
t135 = t175 * qJD(6);
t93 = t175 * t133;
t204 = -t93 - t135;
t223 = pkin(7) + t151;
t139 = t223 * t160;
t140 = t223 * t163;
t246 = -t228 * t139 - t166 * t140;
t188 = qJD(3) * t193;
t195 = qJD(3) * t150 + qJD(4) * t190;
t50 = -t166 * t188 + (qJD(5) - t110) * qJD(4) + t195;
t68 = pkin(4) * t128 - qJ(5) * t127 - qJD(5) * t134;
t18 = -t159 * t50 + t162 * t68;
t14 = pkin(5) * t128 - pkin(8) * t200 + t18;
t19 = t159 * t68 + t162 * t50;
t15 = -pkin(8) * t201 + t19;
t76 = t166 * t117 + t118 * t228;
t67 = qJD(4) * qJ(5) + t76;
t174 = -cos(pkin(9)) * pkin(1) - pkin(3) * t163 - pkin(2);
t132 = qJD(1) * t174 + qJD(3);
t85 = pkin(4) * t133 - qJ(5) * t134 + t132;
t29 = -t159 * t67 + t162 * t85;
t16 = pkin(5) * t133 - pkin(8) * t115 + t29;
t30 = t159 * t85 + t162 * t67;
t20 = pkin(8) * t186 + t30;
t5 = t16 * t167 - t165 * t20;
t1 = qJD(6) * t5 + t14 * t165 + t15 * t167;
t6 = t16 * t165 + t167 * t20;
t2 = -qJD(6) * t6 + t14 * t167 - t15 * t165;
t245 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t38 + Ifges(7,6) * t39;
t218 = Ifges(6,2) * t159;
t220 = Ifges(6,4) * t162;
t181 = -t218 + t220;
t221 = Ifges(6,4) * t159;
t182 = t162 * Ifges(6,1) - t221;
t183 = mrSges(6,1) * t159 + mrSges(6,2) * t162;
t229 = t162 / 0.2e1;
t230 = -t159 / 0.2e1;
t65 = -qJD(4) * pkin(4) + qJD(5) - t75;
t244 = t181 * t186 / 0.2e1 + t182 * t115 / 0.2e1 + t65 * t183 + t132 * mrSges(5,2) + (-t30 * t159 - t29 * t162) * mrSges(6,3) + t254 + Ifges(5,5) * qJD(4) + t252 + (t115 * Ifges(6,4) + Ifges(6,2) * t186 + t133 * Ifges(6,6)) * t230 + (t115 * Ifges(6,1) + Ifges(6,4) * t186 + t133 * Ifges(6,5)) * t229;
t243 = Ifges(7,4) * t241 + Ifges(7,2) * t240 + Ifges(7,6) * t233;
t242 = Ifges(7,1) * t241 + Ifges(7,4) * t240 + Ifges(7,5) * t233;
t239 = -t255 / 0.2e1;
t238 = t255 / 0.2e1;
t237 = -t73 / 0.2e1;
t236 = t73 / 0.2e1;
t131 = qJD(6) + t133;
t232 = -t131 / 0.2e1;
t231 = t131 / 0.2e1;
t227 = Ifges(7,4) * t73;
t77 = t172 * qJD(3) + t246 * qJD(4);
t86 = pkin(4) * t138 - qJ(5) * t137 - qJD(5) * t145;
t34 = t159 * t86 + t162 * t77;
t107 = -t166 * t139 + t140 * t228;
t99 = -pkin(4) * t172 - qJ(5) * t145 + t174;
t48 = t162 * t107 + t159 * t99;
t219 = Ifges(6,5) * t162;
t217 = Ifges(6,6) * t159;
t171 = t145 * qJD(3);
t54 = qJD(1) * t171 + qJD(4) * t76;
t215 = t246 * t54;
t212 = t128 * Ifges(6,5);
t211 = t128 * Ifges(6,6);
t208 = t134 * Ifges(5,4);
t207 = t172 * t54;
t206 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t186 + mrSges(6,2) * t115 + t134 * mrSges(5,3);
t198 = t137 * t159;
t197 = t145 * t159;
t31 = -mrSges(7,1) * t255 + mrSges(7,2) * t73;
t192 = -t31 - t206;
t191 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t33 = -t159 * t77 + t162 * t86;
t187 = t128 * mrSges(5,1) + t127 * mrSges(5,2);
t47 = -t107 * t159 + t162 * t99;
t180 = -t217 + t219;
t179 = t159 * t19 + t162 * t18;
t178 = -t159 * t18 + t162 * t19;
t177 = t159 * t29 - t162 * t30;
t32 = -pkin(5) * t172 - t145 * t226 + t47;
t40 = -pkin(8) * t197 + t48;
t11 = -t165 * t40 + t167 * t32;
t12 = t165 * t32 + t167 * t40;
t176 = -(-t148 * t160 + t156) * t160 + t126 * t163;
t121 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t133;
t90 = -mrSges(6,2) * t133 + mrSges(6,3) * t186;
t91 = mrSges(6,1) * t133 - mrSges(6,3) * t115;
t173 = -t159 * t91 + t162 * t90 + t121;
t78 = qJD(4) * t107 + t171;
t169 = t132 * mrSges(5,1) + t29 * mrSges(6,1) + t5 * mrSges(7,1) - Ifges(5,6) * qJD(4) - t208 / 0.2e1 + t131 * Ifges(7,3) + t73 * Ifges(7,5) + t255 * Ifges(7,6) + t115 * Ifges(6,5) + t186 * Ifges(6,6) - t30 * mrSges(6,2) - t6 * mrSges(7,2) + t257 * t253;
t153 = -pkin(5) * t162 - pkin(4);
t124 = Ifges(7,3) * t128;
t103 = t175 * t145;
t102 = t144 * t145;
t95 = mrSges(6,1) * t128 - mrSges(6,3) * t200;
t94 = -mrSges(6,2) * t128 - mrSges(6,3) * t201;
t81 = pkin(5) * t197 - t246;
t69 = Ifges(7,4) * t255;
t61 = t127 * t182 + t212;
t60 = t127 * t181 + t211;
t56 = pkin(5) * t198 + t78;
t55 = -pkin(5) * t199 + t76;
t53 = (-qJD(4) * t118 - t188) * t166 + t195;
t52 = mrSges(7,1) * t131 - mrSges(7,3) * t73;
t51 = -mrSges(7,2) * t131 + mrSges(7,3) * t255;
t46 = t145 * t135 - t137 * t144;
t45 = -t136 * t145 - t137 * t175;
t44 = -pkin(5) * t186 + t65;
t41 = pkin(5) * t201 + t54;
t27 = -mrSges(7,2) * t128 + mrSges(7,3) * t39;
t26 = mrSges(7,1) * t128 - mrSges(7,3) * t38;
t25 = -pkin(8) * t198 + t34;
t24 = Ifges(7,1) * t73 + Ifges(7,5) * t131 + t69;
t23 = Ifges(7,2) * t255 + Ifges(7,6) * t131 + t227;
t17 = pkin(5) * t138 - t137 * t226 + t33;
t4 = -qJD(6) * t12 - t165 * t25 + t167 * t17;
t3 = qJD(6) * t11 + t165 * t17 + t167 * t25;
t7 = [(t60 * t230 + t61 * t229 + t54 * t183 + t180 * t233 - t179 * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1) * t162 ^ 2 / 0.2e1 + (-t220 + t218 / 0.2e1) * t159) * t127) * t145 + t206 * t78 + (t133 * t191 + t169) * t138 + t174 * t187 + (-t107 * t128 - t127 * t246 - t137 * t75 - t138 * t76 + t145 * t54 + t172 * t53) * mrSges(5,3) - t246 * t87 + (Ifges(7,4) * t45 + Ifges(7,2) * t46) * t238 - t103 * t242 - t102 * t243 + (t137 * t258 - t134 * t138 / 0.2e1 + t172 * t127 - t145 * t128) * Ifges(5,4) + (-Ifges(7,4) * t103 - Ifges(7,2) * t102) * t240 + (-Ifges(7,1) * t103 - Ifges(7,4) * t102) * t241 + (-Ifges(7,5) * t103 - Ifges(7,6) * t102) * t233 + t41 * (mrSges(7,1) * t102 - mrSges(7,2) * t103) + (-t1 * t102 + t103 * t2 - t45 * t5 + t46 * t6) * mrSges(7,3) + m(7) * (t1 * t12 + t11 * t2 + t3 * t6 + t4 * t5 + t41 * t81 + t44 * t56) + m(6) * (t18 * t47 + t19 * t48 + t29 * t33 + t30 * t34 + t65 * t78 - t215) + m(5) * (t107 * t53 - t75 * t78 + t76 * t77 - t215) + (Ifges(7,5) * t45 + Ifges(7,6) * t46) * t231 + (Ifges(7,1) * t45 + Ifges(7,4) * t46) * t236 + (m(4) * (t151 * t185 + t176) + 0.2e1 * t256) * qJD(3) - (t124 / 0.2e1 - t19 * mrSges(6,2) + t18 * mrSges(6,1) + t180 * t127 + (Ifges(7,3) / 0.2e1 + t257) * t128 + t245) * t172 + t11 * t26 + t12 * t27 + t45 * t24 / 0.2e1 + t46 * t23 / 0.2e1 + t44 * (-mrSges(7,1) * t46 + mrSges(7,2) * t45) + t3 * t51 + t4 * t52 + t56 * t31 + t81 * t13 + t34 * t90 + t33 * t91 + t48 * t94 + t47 * t95 + t77 * t121 + (t180 * t253 + t244 + t252) * t137; -t102 * t26 - t103 * t27 + t45 * t51 + t46 * t52 + (-t128 * mrSges(5,3) - t159 * t95 + t162 * t94) * t145 - (mrSges(5,3) * t127 + t249) * t172 - t192 * t138 + t173 * t137 + m(5) * (t137 * t76 - t138 * t75 + t145 * t53 - t207) + m(6) * (-t137 * t177 + t138 * t65 + t145 * t178 - t207) + m(7) * (-t1 * t103 - t102 * t2 + t138 * t44 - t172 * t41 + t45 * t6 + t46 * t5); -t175 * t26 + t144 * t27 + t159 * t94 + t162 * t95 + t205 * t52 + t204 * t51 + t192 * t134 + t173 * t133 - m(5) * (-t133 * t76 - t134 * t75) + t187 + (-m(4) * t176 - t256) * qJD(1) + (t1 * t144 - t134 * t44 - t175 * t2 + t204 * t6 + t205 * t5) * m(7) + (-t133 * t177 - t134 * t65 + t179) * m(6); ((Ifges(6,1) * t159 + t220) * t229 + (Ifges(6,2) * t162 + t221) * t230 + Ifges(5,5)) * t127 + (t19 * mrSges(6,3) + qJD(5) * t90 + qJ(5) * t94 + t60 / 0.2e1 - t54 * mrSges(6,1) + t211 / 0.2e1) * t162 + (-t18 * mrSges(6,3) - qJD(5) * t91 - qJ(5) * t95 + t61 / 0.2e1 + t54 * mrSges(6,2) + t212 / 0.2e1) * t159 + (-t169 + t208 / 0.2e1 + t76 * mrSges(5,3)) * t134 + (-mrSges(7,1) * t205 + mrSges(7,2) * t204) * t44 - t206 * t76 - t175 * t243 + (-t1 * t175 - t144 * t2 - t204 * t5 + t205 * t6) * mrSges(7,3) + (Ifges(7,4) * t144 - Ifges(7,2) * t175) * t240 + (Ifges(7,1) * t144 - Ifges(7,4) * t175) * t241 + (Ifges(7,4) * t93 + Ifges(7,2) * t92) * t239 + t144 * t242 + (-t75 * mrSges(5,3) + t254 + (t219 / 0.2e1 - t217 / 0.2e1) * t133 + (t259 - t191) * t134 + t244) * t133 + (-t93 / 0.2e1 - t135 / 0.2e1) * t24 + (Ifges(7,5) * t144 - Ifges(7,6) * t175) * t233 + t41 * (mrSges(7,1) * t175 + mrSges(7,2) * t144) + (Ifges(7,5) * t93 + Ifges(7,6) * t92) * t232 + (Ifges(7,1) * t93 + Ifges(7,4) * t92) * t237 + (-pkin(4) * t54 + qJ(5) * t178 - qJD(5) * t177 - t29 * t42 - t30 * t43 - t65 * t76) * m(6) - t53 * mrSges(5,2) - t54 * mrSges(5,1) - t55 * t31 + (-Ifges(7,4) * t135 - Ifges(7,2) * t136) * t238 + (-Ifges(7,5) * t135 - Ifges(7,6) * t136) * t231 + (-Ifges(7,1) * t135 - Ifges(7,4) * t136) * t236 + (-t92 / 0.2e1 - t136 / 0.2e1) * t23 - pkin(4) * t87 - t43 * t90 - t42 * t91 + t111 * t26 + t112 * t27 - t75 * t121 - Ifges(5,6) * t128 + t153 * t13 + t250 * t51 + t251 * t52 + (t1 * t112 + t111 * t2 + t153 * t41 + t250 * t6 + t251 * t5 - t44 * t55) * m(7); t115 * t91 - t186 * t90 - t255 * t51 + t73 * t52 + (-t255 * t6 + t5 * t73 + t41) * m(7) + (t115 * t29 - t186 * t30 + t54) * m(6) + t249; t124 - t44 * (mrSges(7,1) * t73 + mrSges(7,2) * t255) + (Ifges(7,1) * t255 - t227) * t237 + t23 * t236 + (Ifges(7,5) * t255 - Ifges(7,6) * t73) * t232 - t5 * t51 + t6 * t52 + (t255 * t5 + t6 * t73) * mrSges(7,3) + (-Ifges(7,2) * t73 + t24 + t69) * t239 + t245;];
tauc  = t7(:);
