% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2018-11-23 15:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:55:27
% EndTime: 2018-11-23 15:55:31
% DurationCPUTime: 4.52s
% Computational Cost: add. (3838->423), mult. (8590->548), div. (0->0), fcn. (5413->6), ass. (0->209)
t253 = -mrSges(5,1) + mrSges(6,2);
t130 = sin(qJ(3));
t132 = cos(qJ(3));
t195 = sin(pkin(9));
t196 = cos(pkin(9));
t102 = -t195 * t130 + t196 * t132;
t131 = cos(qJ(6));
t129 = sin(qJ(6));
t206 = Ifges(7,4) * t129;
t151 = Ifges(7,2) * t131 + t206;
t205 = Ifges(7,4) * t131;
t153 = Ifges(7,1) * t129 + t205;
t156 = mrSges(7,1) * t131 - mrSges(7,2) * t129;
t101 = t130 * t196 + t132 * t195;
t94 = t101 * qJD(1);
t71 = qJD(3) * t131 + t129 * t94;
t224 = Ifges(7,4) * t71;
t70 = -qJD(3) * t129 + t131 * t94;
t96 = t102 * qJD(1);
t93 = qJD(6) + t96;
t19 = Ifges(7,2) * t70 + Ifges(7,6) * t93 + t224;
t69 = Ifges(7,4) * t70;
t20 = Ifges(7,1) * t71 + Ifges(7,5) * t93 + t69;
t202 = Ifges(7,6) * t131;
t203 = Ifges(7,5) * t129;
t227 = -t129 / 0.2e1;
t236 = -t93 / 0.2e1;
t238 = -t71 / 0.2e1;
t239 = -t70 / 0.2e1;
t229 = t94 * pkin(5);
t189 = qJD(1) * t130;
t133 = -pkin(1) - pkin(7);
t112 = qJD(1) * t133 + qJD(2);
t192 = t112 * t130;
t89 = -qJ(4) * t189 + t192;
t172 = t196 * t89;
t104 = t132 * t112;
t188 = qJD(1) * t132;
t90 = -qJ(4) * t188 + t104;
t81 = qJD(3) * pkin(3) + t90;
t45 = t195 * t81 + t172;
t41 = -qJD(3) * qJ(5) - t45;
t25 = -t41 - t229;
t256 = -t131 * t19 / 0.2e1 + t20 * t227 + (t202 + t203) * t236 + t153 * t238 + t151 * t239 + t25 * t156;
t106 = pkin(3) * t189 + qJD(1) * qJ(2) + qJD(4);
t146 = -qJ(5) * t96 + t106;
t46 = pkin(4) * t94 + t146;
t254 = -t106 * mrSges(5,1) + t46 * mrSges(6,2) + t256;
t252 = Ifges(6,4) - Ifges(5,5);
t251 = Ifges(6,5) - Ifges(5,6);
t218 = t94 * mrSges(5,3);
t225 = mrSges(6,1) * t94;
t76 = -qJD(3) * mrSges(6,3) + t225;
t210 = -qJD(3) * mrSges(5,2) - t218 - t76;
t36 = -mrSges(7,1) * t70 + mrSges(7,2) * t71;
t250 = -t76 + t36;
t216 = t96 * mrSges(5,3);
t217 = t96 * mrSges(6,1);
t209 = t253 * qJD(3) + t216 + t217;
t249 = qJ(2) * (m(3) + m(4));
t181 = qJD(1) * qJD(3);
t171 = t132 * t181;
t85 = t96 * qJD(3);
t42 = qJD(6) * t70 + t129 * t85;
t97 = t101 * qJD(3);
t86 = qJD(1) * t97;
t21 = -mrSges(7,1) * t86 - mrSges(7,3) * t42;
t43 = -qJD(6) * t71 + t131 * t85;
t22 = mrSges(7,2) * t86 + mrSges(7,3) * t43;
t247 = -t129 * t22 - t131 * t21;
t105 = pkin(3) * t171 + qJD(1) * qJD(2);
t141 = qJ(5) * t86 - qJD(5) * t96 + t105;
t230 = pkin(4) + pkin(8);
t14 = t230 * t85 + t141;
t187 = qJD(3) * t130;
t174 = t112 * t187;
t184 = qJD(4) * t132;
t135 = -t174 + (qJ(4) * t187 - t184) * qJD(1);
t186 = qJD(3) * t132;
t173 = t112 * t186;
t185 = qJD(4) * t130;
t68 = t173 + (-qJ(4) * t186 - t185) * qJD(1);
t30 = -t196 * t135 + t195 * t68;
t16 = -t86 * pkin(5) + t30;
t78 = t195 * t89;
t44 = t196 * t81 - t78;
t145 = qJD(5) - t44;
t228 = t96 * pkin(5);
t23 = -qJD(3) * t230 + t145 + t228;
t27 = t230 * t94 + t146;
t5 = -t129 * t27 + t131 * t23;
t1 = qJD(6) * t5 + t129 * t16 + t131 * t14;
t6 = t129 * t23 + t131 * t27;
t2 = -qJD(6) * t6 - t129 * t14 + t131 * t16;
t246 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t42 + Ifges(7,6) * t43;
t162 = t5 * t129 - t6 * t131;
t245 = qJD(6) * t162 - t1 * t129 - t131 * t2;
t243 = t5 * mrSges(7,1) + t106 * mrSges(5,2) - t6 * mrSges(7,2) - t46 * mrSges(6,3);
t242 = qJD(1) ^ 2;
t241 = t42 / 0.2e1;
t240 = t43 / 0.2e1;
t237 = t71 / 0.2e1;
t235 = -t94 / 0.2e1;
t233 = -t96 / 0.2e1;
t232 = t96 / 0.2e1;
t226 = t131 / 0.2e1;
t190 = qJ(4) - t133;
t107 = t190 * t130;
t108 = t190 * t132;
t66 = -t107 * t195 + t196 * t108;
t223 = t30 * t66;
t222 = t70 * Ifges(7,6);
t221 = t71 * Ifges(7,5);
t220 = t86 * mrSges(6,1);
t219 = t93 * Ifges(7,3);
t215 = t96 * Ifges(5,4);
t214 = t96 * Ifges(6,6);
t213 = -qJD(3) / 0.2e1;
t212 = mrSges(6,1) + mrSges(5,3);
t211 = Ifges(5,4) + Ifges(6,6);
t31 = t195 * t135 + t196 * t68;
t208 = Ifges(4,4) * t130;
t207 = Ifges(4,4) * t132;
t201 = t102 * t30;
t200 = t129 * mrSges(7,3);
t198 = t131 * mrSges(7,3);
t194 = Ifges(4,5) * qJD(3);
t193 = Ifges(4,6) * qJD(3);
t121 = t130 * pkin(3) + qJ(2);
t183 = qJD(6) * t129;
t182 = qJD(6) * t131;
t113 = pkin(3) * t186 + qJD(2);
t180 = -t36 - t210;
t124 = pkin(3) * t188;
t179 = Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1;
t176 = t196 * pkin(3);
t175 = t195 * pkin(3);
t168 = qJ(5) * t94 + t124;
t120 = -t176 - pkin(4);
t165 = t1 * t131 - t2 * t129;
t163 = t6 * t129 + t5 * t131;
t39 = -qJD(3) * pkin(4) + t145;
t95 = t102 * qJD(3);
t161 = t39 * t97 - t41 * t95;
t160 = -t44 * t97 + t45 * t95;
t67 = -t107 * t196 - t108 * t195;
t159 = -t66 * t86 - t67 * t85;
t158 = -qJ(5) * t102 + t121;
t26 = -qJD(3) * qJD(5) - t31;
t87 = t187 * t190 - t184;
t88 = -qJD(3) * t108 - t185;
t49 = t195 * t88 - t196 * t87;
t53 = t195 * t90 + t172;
t157 = mrSges(4,1) * t130 + mrSges(4,2) * t132;
t155 = mrSges(7,1) * t129 + mrSges(7,2) * t131;
t154 = Ifges(7,1) * t131 - t206;
t152 = -Ifges(7,2) * t129 + t205;
t150 = Ifges(7,5) * t131 - Ifges(7,6) * t129;
t40 = t101 * t230 + t158;
t51 = t102 * pkin(5) + t66;
t13 = t129 * t51 + t131 * t40;
t12 = -t129 * t40 + t131 * t51;
t47 = -mrSges(7,2) * t93 + mrSges(7,3) * t70;
t48 = mrSges(7,1) * t93 - mrSges(7,3) * t71;
t148 = t129 * t48 - t131 * t47;
t147 = -t129 * t47 - t131 * t48;
t143 = qJ(2) * (mrSges(4,1) * t132 - mrSges(4,2) * t130);
t142 = -t147 + t209;
t140 = qJ(5) * t97 - qJD(5) * t102 + t113;
t50 = t195 * t87 + t196 * t88;
t54 = t196 * t90 - t78;
t118 = t175 + qJ(5);
t110 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t188;
t109 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t189;
t103 = t157 * qJD(1);
t99 = t194 + (t132 * Ifges(4,1) - t208) * qJD(1);
t98 = t193 + (-t130 * Ifges(4,2) + t207) * qJD(1);
t92 = Ifges(5,4) * t94;
t91 = Ifges(6,6) * t94;
t84 = Ifges(7,3) * t86;
t83 = t86 * mrSges(5,2);
t82 = t86 * mrSges(6,3);
t62 = pkin(4) * t101 + t158;
t61 = -mrSges(6,2) * t94 - mrSges(6,3) * t96;
t60 = mrSges(5,1) * t94 + mrSges(5,2) * t96;
t59 = t96 * Ifges(5,1) + Ifges(5,5) * qJD(3) - t92;
t58 = -t94 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t215;
t57 = Ifges(6,4) * qJD(3) - t96 * Ifges(6,2) + t91;
t56 = Ifges(6,5) * qJD(3) + t94 * Ifges(6,3) - t214;
t55 = pkin(4) * t96 + t168;
t52 = -t101 * pkin(5) + t67;
t35 = pkin(4) * t95 + t140;
t34 = t230 * t96 + t168;
t33 = t54 - t228;
t32 = t53 - t229;
t29 = -t95 * pkin(5) + t50;
t28 = -t97 * pkin(5) + t49;
t24 = pkin(4) * t85 + t141;
t18 = t219 + t221 + t222;
t17 = t230 * t95 + t140;
t15 = -pkin(5) * t85 - t26;
t11 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t10 = Ifges(7,1) * t42 + Ifges(7,4) * t43 - Ifges(7,5) * t86;
t9 = Ifges(7,4) * t42 + Ifges(7,2) * t43 - Ifges(7,6) * t86;
t8 = t129 * t32 + t131 * t34;
t7 = -t129 * t34 + t131 * t32;
t4 = -qJD(6) * t13 - t129 * t17 + t131 * t28;
t3 = qJD(6) * t12 + t129 * t28 + t131 * t17;
t37 = [t35 * t61 + t3 * t47 + t4 * t48 + t52 * t11 + t29 * t36 + t12 * t21 + t13 * t22 + t62 * (-t85 * mrSges(6,2) + t82) + qJD(2) * t103 + t113 * t60 + t121 * (t85 * mrSges(5,1) - t83) + (-t222 / 0.2e1 - t221 / 0.2e1 - t219 / 0.2e1 + t57 / 0.2e1 - t59 / 0.2e1 - t18 / 0.2e1 + (-Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t96 + t179 * t94 - t243) * t97 - (t58 / 0.2e1 - t56 / 0.2e1 + t179 * t96 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1) * t94 + t162 * mrSges(7,3) + t254) * t95 + t210 * t50 + t209 * t49 + m(6) * (t24 * t62 - t26 * t67 + t35 * t46 + t39 * t49 - t41 * t50 + t223) + m(7) * (t1 * t13 + t12 * t2 + t15 * t52 + t25 * t29 + t3 * t6 + t4 * t5) + m(5) * (t105 * t121 + t106 * t113 + t31 * t67 - t44 * t49 + t45 * t50 + t223) + (t159 - t160) * mrSges(5,3) + (t159 - t161) * mrSges(6,1) + ((Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t97 - (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t95 + (-t98 / 0.2e1 + t133 * t109 - t193 / 0.2e1) * t132 + (-t99 / 0.2e1 - t133 * t110 - t194 / 0.2e1) * t130) * qJD(3) + (-t84 / 0.2e1 - t24 * mrSges(6,3) + t105 * mrSges(5,2) - t211 * t85 + t212 * t30 + (-Ifges(7,3) / 0.2e1 - Ifges(5,1) - Ifges(6,2)) * t86 + t246) * t102 + (-t24 * mrSges(6,2) + t105 * mrSges(5,1) + t26 * mrSges(6,1) - t31 * mrSges(5,3) + t129 * t10 / 0.2e1 + t9 * t226 - t15 * t156 + t153 * t241 + t151 * t240 + (-t203 / 0.2e1 - t202 / 0.2e1 + t211) * t86 - (-Ifges(6,3) - Ifges(5,2)) * t85 + t165 * mrSges(7,3) + (t25 * t155 + t70 * t152 / 0.2e1 + t154 * t237 + t93 * t150 / 0.2e1 + t19 * t227 + t20 * t226 - t163 * mrSges(7,3)) * qJD(6)) * t101 + (((2 * mrSges(3,3)) + t157 + 0.2e1 * t249) * qJD(2) + (0.2e1 * t143 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t130 * t132 + (-0.3e1 / 0.2e1 * t132 ^ 2 + 0.3e1 / 0.2e1 * t130 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1); (t109 * t132 - t110 * t130) * qJD(3) - t180 * t95 + (-mrSges(3,3) - t249) * t242 + (-t212 * t85 + t11) * t101 + t142 * t97 + (t148 * qJD(6) + t212 * t86 + t247) * t102 + m(6) * (-t101 * t26 + t161 - t201) + m(5) * (t101 * t31 + t160 - t201) + m(7) * (t101 * t15 + t102 * t245 + t163 * t97 + t25 * t95) + (-m(5) * t106 - m(6) * t46 + m(7) * t162 - t103 + t148 - t60 - t61) * qJD(1); (t56 - t215) * t233 + (t58 + t214) * t232 + (-t132 * (-Ifges(4,1) * t130 - t207) / 0.2e1 + t130 * (-Ifges(4,2) * t132 - t208) / 0.2e1 - t143) * t242 + (-t6 * t182 + t5 * t183) * mrSges(7,3) + ((t195 * t31 - t196 * t30) * pkin(3) - t106 * t124 + t44 * t53 - t45 * t54) * m(5) + (Ifges(6,3) * t235 - t6 * t198 + t5 * t200 + t251 * t213 + t254) * t96 - t181 * Ifges(4,5) * t130 / 0.2e1 + (-Ifges(5,2) * t96 + t18 + t59 - t92) * t94 / 0.2e1 + (-t118 * t26 + t120 * t30 - t39 * t53 - t46 * t55 + (t54 - qJD(5)) * t41) * m(6) + (t118 * t15 - t5 * t7 - t6 * t8 + (-t33 + qJD(5)) * t25) * m(7) + (-m(7) * t245 + t182 * t47 - t183 * t48 - t247) * (-pkin(8) + t120) + t256 * qJD(6) - t55 * t61 - t7 * t48 - t8 * t47 - t31 * mrSges(5,2) - t33 * t36 - t26 * mrSges(6,3) + (t57 + t91) * t235 - t210 * t54 - (mrSges(6,1) * t118 + mrSges(5,3) * t175 - t251) * t85 + (-t150 / 0.2e1 + mrSges(5,3) * t176 + t252) * t86 + (-Ifges(5,1) * t233 - Ifges(7,5) * t238 + Ifges(6,2) * t232 - Ifges(7,6) * t239 - Ifges(7,3) * t236 + t252 * t213 + t243) * t94 + t253 * t30 + t9 * t227 + t152 * t240 + t154 * t241 + t45 * t216 + t39 * t225 + t10 * t226 + t110 * t192 - t209 * t53 + t250 * qJD(5) + t15 * t155 - Ifges(4,6) * t171 / 0.2e1 - mrSges(4,2) * t173 - mrSges(4,1) * t174 - t60 * t124 + t98 * t188 / 0.2e1 + t99 * t189 / 0.2e1 - t109 * t104 - t2 * t198 - t1 * t200 + t118 * t11 - t41 * t217 - t44 * t218 - t120 * t220; -t129 * t21 + t131 * t22 + t82 - t83 - t253 * t85 + t147 * qJD(6) - t180 * t94 - t142 * t96 + (-t163 * t93 + t25 * t94 + t165) * m(7) + (-t39 * t96 - t41 * t94 + t24) * m(6) + (t44 * t96 + t45 * t94 + t105) * m(5); -t220 + t96 * t61 - t250 * qJD(3) + (t47 * t93 + t21) * t131 + (-t48 * t93 + t22) * t129 + (-qJD(3) * t25 - t162 * t96 - t245) * m(7) + (qJD(3) * t41 + t46 * t96 + t30) * m(6); -t84 - t25 * (mrSges(7,1) * t71 + mrSges(7,2) * t70) + (Ifges(7,1) * t70 - t224) * t238 + t19 * t237 + (Ifges(7,5) * t70 - Ifges(7,6) * t71) * t236 - t5 * t47 + t6 * t48 + (t5 * t70 + t6 * t71) * mrSges(7,3) + (-Ifges(7,2) * t71 + t20 + t69) * t239 + t246;];
tauc  = t37(:);
