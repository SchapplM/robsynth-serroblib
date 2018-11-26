% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2018-11-23 16:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:02:27
% EndTime: 2018-11-23 16:02:32
% DurationCPUTime: 4.56s
% Computational Cost: add. (2951->441), mult. (6148->545), div. (0->0), fcn. (2953->4), ass. (0->203)
t256 = -qJD(1) / 0.2e1;
t268 = Ifges(6,1) + Ifges(7,1);
t253 = Ifges(7,4) + Ifges(6,5);
t128 = sin(qJ(5));
t130 = cos(qJ(5));
t132 = -pkin(3) - pkin(8);
t131 = cos(qJ(3));
t186 = t131 * qJD(1);
t133 = -pkin(1) - pkin(7);
t109 = qJD(1) * t133 + qJD(2);
t99 = t131 * t109;
t75 = -pkin(4) * t186 + t99;
t259 = -t75 + qJD(4);
t53 = qJD(3) * t132 + t259;
t129 = sin(qJ(3));
t203 = qJ(4) * t131;
t153 = pkin(8) * t129 - t203;
t193 = qJD(1) * t129;
t194 = pkin(3) * t193 + qJD(1) * qJ(2);
t65 = qJD(1) * t153 + t194;
t17 = -t128 * t65 + t130 * t53;
t18 = t128 * t53 + t130 * t65;
t151 = t17 * t128 - t18 * t130;
t115 = qJD(5) + t186;
t249 = qJD(6) - t17;
t10 = -pkin(5) * t115 + t249;
t11 = qJ(6) * t115 + t18;
t152 = t10 * t128 + t11 * t130;
t167 = mrSges(7,1) * t130 + mrSges(7,3) * t128;
t169 = mrSges(6,1) * t130 - mrSges(6,2) * t128;
t126 = qJD(3) * qJ(4);
t98 = t129 * t109;
t74 = -pkin(4) * t193 + t98;
t67 = t126 + t74;
t192 = qJD(3) * t128;
t92 = -t130 * t193 + t192;
t190 = qJD(3) * t130;
t93 = t128 * t193 + t190;
t21 = pkin(5) * t92 - qJ(6) * t93 + t67;
t227 = t130 / 0.2e1;
t228 = -t130 / 0.2e1;
t230 = -t128 / 0.2e1;
t223 = Ifges(7,5) * t92;
t88 = Ifges(6,4) * t92;
t251 = t253 * t115 + t268 * t93 + t223 - t88;
t87 = Ifges(7,5) * t93;
t26 = t115 * Ifges(7,6) + t92 * Ifges(7,3) + t87;
t224 = Ifges(6,4) * t93;
t29 = -t92 * Ifges(6,2) + t115 * Ifges(6,6) + t224;
t269 = -t152 * mrSges(7,2) + t151 * mrSges(6,3) + t21 * t167 + t67 * t169 + t227 * t26 + t228 * t29 + t230 * t251;
t267 = qJD(4) - t99;
t225 = mrSges(6,3) * t93;
t61 = mrSges(6,1) * t115 - t225;
t62 = -mrSges(7,1) * t115 + mrSges(7,2) * t93;
t217 = t61 - t62;
t226 = mrSges(6,3) * t92;
t60 = -mrSges(6,2) * t115 - t226;
t63 = -mrSges(7,2) * t92 + mrSges(7,3) * t115;
t218 = t60 + t63;
t141 = t128 * t217 - t130 * t218;
t266 = -m(6) * t151 + m(7) * t152 - t141;
t210 = Ifges(7,5) * t128;
t156 = -Ifges(7,3) * t130 + t210;
t214 = Ifges(6,4) * t128;
t160 = Ifges(6,2) * t130 + t214;
t231 = t115 / 0.2e1;
t234 = t93 / 0.2e1;
t236 = t92 / 0.2e1;
t237 = -t92 / 0.2e1;
t209 = Ifges(7,5) * t130;
t213 = Ifges(6,4) * t130;
t246 = t128 * t268 - t209 + t213;
t207 = Ifges(7,6) * t130;
t208 = Ifges(6,6) * t130;
t211 = Ifges(6,5) * t128;
t212 = Ifges(7,4) * t128;
t247 = -t207 + t212 + t208 + t211;
t254 = qJD(3) / 0.2e1;
t255 = -qJD(3) / 0.2e1;
t79 = -qJ(4) * t186 + t194;
t84 = -t98 - t126;
t265 = -t84 * mrSges(5,1) + t79 * mrSges(5,2) - Ifges(5,5) * t254 - Ifges(4,6) * t255 - t156 * t236 - t160 * t237 - t247 * t231 - t246 * t234 + t269 + ((-Ifges(4,4) - Ifges(5,6)) * t131 + (Ifges(4,2) + Ifges(5,3)) * t129) * t256;
t185 = qJD(1) * qJD(3);
t174 = t129 * t185;
t173 = t131 * t185;
t57 = qJD(5) * t93 - t130 * t173;
t36 = -mrSges(7,2) * t57 - mrSges(7,3) * t174;
t187 = qJD(5) * t130;
t188 = qJD(5) * t128;
t189 = qJD(3) * t131;
t56 = -qJD(3) * t188 + (t128 * t189 + t129 * t187) * qJD(1);
t37 = -mrSges(6,1) * t174 - mrSges(6,3) * t56;
t38 = mrSges(7,1) * t174 + t56 * mrSges(7,2);
t39 = mrSges(6,2) * t174 - mrSges(6,3) * t57;
t142 = (t36 + t39) * t128 + (t37 - t38) * t130;
t264 = t141 * qJD(5) - t142;
t263 = Ifges(7,2) + Ifges(6,3);
t262 = (-Ifges(6,4) + Ifges(7,5)) * t57 + t268 * t56 - t253 * t174;
t155 = pkin(5) * t130 + qJ(6) * t128;
t146 = -pkin(4) - t155;
t143 = t146 * t131;
t261 = -qJD(1) * t143 + qJD(5) * t155 - qJD(6) * t130 + t267;
t260 = t130 * t268 + t210 - t214;
t258 = qJD(1) * (qJ(2) * (m(3) + m(4)) + mrSges(3,3));
t257 = (-mrSges(4,1) + mrSges(5,2)) * qJD(3);
t235 = -t93 / 0.2e1;
t232 = -t115 / 0.2e1;
t252 = -Ifges(6,6) + Ifges(7,6);
t250 = t129 * pkin(3) + qJ(2);
t147 = (qJD(3) * pkin(8) - qJD(4)) * t131;
t125 = qJD(1) * qJD(2);
t177 = pkin(3) * t173 + qJ(4) * t174 + t125;
t42 = qJD(1) * t147 + t177;
t191 = qJD(3) * t129;
t69 = (-pkin(4) * qJD(1) + t109) * t191;
t3 = t128 * t69 + t130 * t42 + t53 * t187 - t188 * t65;
t4 = -qJD(5) * t18 - t128 * t42 + t130 * t69;
t170 = t3 * t128 + t4 * t130;
t245 = -qJD(3) * t67 + t170;
t1 = -qJ(6) * t174 + qJD(6) * t115 + t3;
t2 = pkin(5) * t174 - t4;
t171 = t1 * t128 - t2 * t130;
t244 = -qJD(3) * t21 + t171;
t221 = pkin(4) - t133;
t103 = t221 * t131;
t86 = t153 + t250;
t216 = t128 * t103 + t130 * t86;
t175 = pkin(3) * t189 + qJ(4) * t191 + qJD(2);
t58 = t147 + t175;
t90 = t221 * t191;
t9 = -qJD(5) * t216 - t128 * t58 - t130 * t90;
t118 = Ifges(5,6) * t193;
t178 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t181 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t215 = Ifges(4,4) * t129;
t241 = -(Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t115 + t178 * t92 + t181 * t93 + t10 * mrSges(7,1) + t18 * mrSges(6,2) + t79 * mrSges(5,3) + Ifges(6,6) * t236 + Ifges(7,6) * t237 + Ifges(4,5) * t255 + (t131 * Ifges(4,1) - t215) * t256 + Ifges(5,4) * t254 - Ifges(5,2) * t186 / 0.2e1 + t118 / 0.2e1 - t11 * mrSges(7,3) - t17 * mrSges(6,1) + t253 * t235 + t263 * t232;
t240 = t56 / 0.2e1;
t239 = -t57 / 0.2e1;
t238 = t57 / 0.2e1;
t229 = t128 / 0.2e1;
t95 = pkin(3) * t186 + qJ(4) * t193;
t73 = pkin(8) * t186 + t95;
t25 = t128 * t74 + t130 * t73;
t206 = qJ(2) * mrSges(4,1);
t205 = qJ(2) * mrSges(4,2);
t106 = mrSges(5,1) * t193 - qJD(3) * mrSges(5,3);
t45 = mrSges(6,1) * t92 + mrSges(6,2) * t93;
t204 = -t106 + t45;
t76 = -qJD(3) * qJD(4) - t109 * t189;
t202 = qJD(3) * mrSges(4,2);
t199 = qJD(3) * t84;
t198 = t128 * t132;
t197 = t130 * t132;
t104 = -mrSges(4,3) * t193 - t202;
t196 = t104 - t106;
t195 = (mrSges(5,1) + mrSges(4,3)) * t186 + t257;
t44 = mrSges(7,1) * t92 - mrSges(7,3) * t93;
t184 = t44 + t204;
t183 = -Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t182 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t179 = -0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,4);
t176 = t109 * t191;
t78 = -qJD(3) * pkin(3) + t267;
t172 = -t78 + t99;
t168 = mrSges(6,1) * t128 + mrSges(6,2) * t130;
t166 = mrSges(7,1) * t128 - mrSges(7,3) * t130;
t161 = -Ifges(6,2) * t128 + t213;
t157 = Ifges(7,3) * t128 + t209;
t154 = pkin(5) * t128 - qJ(6) * t130;
t24 = -t128 * t73 + t130 * t74;
t40 = t103 * t130 - t128 * t86;
t8 = t103 * t187 - t128 * t90 + t130 * t58 - t188 * t86;
t59 = -pkin(4) * t173 - t76;
t139 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t119 = t129 * t133;
t113 = t133 * t189;
t102 = -pkin(4) * t129 + t119;
t101 = -t203 + t250;
t100 = qJ(4) + t154;
t97 = qJD(1) * (mrSges(4,1) * t129 + mrSges(4,2) * t131);
t96 = (-mrSges(5,2) * t129 - mrSges(5,3) * t131) * qJD(1);
t91 = -pkin(4) * t189 + t113;
t71 = -qJD(4) * t131 + t175;
t55 = -qJD(4) * t186 + t177;
t54 = t129 * t146 + t119;
t52 = Ifges(7,4) * t56;
t51 = Ifges(6,5) * t56;
t50 = Ifges(6,6) * t57;
t49 = Ifges(7,6) * t57;
t43 = pkin(5) * t93 + qJ(6) * t92;
t33 = -pkin(5) * t131 - t40;
t32 = qJ(6) * t131 + t216;
t23 = pkin(5) * t193 - t24;
t22 = -qJ(6) * t193 + t25;
t20 = mrSges(6,1) * t57 + mrSges(6,2) * t56;
t19 = mrSges(7,1) * t57 - mrSges(7,3) * t56;
t16 = t113 + (qJD(5) * t154 - qJD(6) * t128) * t129 + qJD(3) * t143;
t13 = t56 * Ifges(6,4) - t57 * Ifges(6,2) - Ifges(6,6) * t174;
t12 = t56 * Ifges(7,5) - Ifges(7,6) * t174 + t57 * Ifges(7,3);
t7 = pkin(5) * t191 - t9;
t6 = -qJ(6) * t191 + qJD(6) * t131 + t8;
t5 = pkin(5) * t57 - qJ(6) * t56 - qJD(6) * t93 + t59;
t14 = [t102 * t20 + t16 * t44 + t54 * t19 + t32 * t36 + t33 * t38 + t40 * t37 + t216 * t39 + t91 * t45 + t6 * t63 + t8 * t60 + t9 * t61 + t7 * t62 + t71 * t96 + m(5) * (t101 * t55 + t71 * t79) + m(6) * (t102 * t59 + t17 * t9 + t18 * t8 + t216 * t3 + t4 * t40 + t67 * t91) + m(7) * (t1 * t32 + t10 * t7 + t11 * t6 + t16 * t21 + t2 * t33 + t5 * t54) + (t97 + 0.2e1 * t258) * qJD(2) + (mrSges(4,2) * t125 - t55 * mrSges(5,3) + t52 / 0.2e1 - t50 / 0.2e1 + t51 / 0.2e1 + t49 / 0.2e1 - t178 * t57 - t181 * t56 + ((-m(5) * t84 + t196) * t133 + t182 * qJD(3) + (-t101 * mrSges(5,2) + t131 * t179 + 0.2e1 * t206) * qJD(1) - t265) * qJD(3) + t139) * t131 + (mrSges(4,1) * t125 + t12 * t228 + t13 * t227 - t59 * t169 - t5 * t167 + t156 * t238 + t160 * t239 - t55 * mrSges(5,2) + (-m(5) * t133 + mrSges(5,1)) * t76 + (-t128 * t4 + t130 * t3) * mrSges(6,3) + (t1 * t130 + t128 * t2) * mrSges(7,2) + (t29 * t230 + t67 * t168 + t21 * t166 + t157 * t236 + t161 * t237 + (-t128 * t18 - t130 * t17) * mrSges(6,3) + (t10 * t130 - t11 * t128) * mrSges(7,2) + t260 * t234 + (t128 * t252 + t130 * t253) * t231 + t251 * t227) * qJD(5) + ((-m(5) * t172 + t195) * t133 + t172 * mrSges(5,1) + t183 * qJD(3) + ((-t211 / 0.2e1 - t208 / 0.2e1 - t212 / 0.2e1 + t207 / 0.2e1 - t179) * t129 + t101 * mrSges(5,3) - 0.2e1 * t205 + (0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,1) - t263) * t131) * qJD(1) + t241) * qJD(3) + t246 * t240 + (qJD(5) * t26 + t262) * t229) * t129; (t19 + t20 + (t128 * t218 + t130 * t217 + t195) * qJD(3) + m(6) * (t17 * t190 + t18 * t192 + t59) + m(7) * (-t10 * t190 + t11 * t192 + t5) + m(5) * (qJD(3) * t78 - t76)) * t129 + ((t104 + t184) * qJD(3) + m(6) * (t17 * t188 - t18 * t187 - t245) + m(7) * (-t10 * t188 - t11 * t187 - t244) + m(5) * (-t176 - t199) + t264) * t131 + (-m(5) * t79 - t258 - t266 - t96 - t97) * qJD(1); (t132 * t266 + t156 * t237 + t160 * t236 + t247 * t232 + t246 * t235 + t269) * qJD(5) + t204 * qJD(4) + (-pkin(3) * t176 - qJ(4) * t76 - qJD(4) * t84 - t79 * t95 - (t129 * t78 - t131 * t84) * t109) * m(5) - t95 * t96 + t100 * t19 - t75 * t45 - t76 * mrSges(5,3) - t25 * t60 - t24 * t61 - t23 * t62 - t22 * t63 + qJ(4) * t20 + (((pkin(3) * mrSges(5,1) + t128 * t178 + t130 * t181 + t183) * qJD(3) - t118 / 0.2e1 + (t205 - t215 / 0.2e1) * qJD(1) + t78 * mrSges(5,1) - t241) * t129 + ((-t206 + (Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t131) * qJD(1) + (-qJ(4) * mrSges(5,1) + t182) * qJD(3) + (Ifges(5,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t193 + t265) * t131) * qJD(1) + (t1 * t198 - t10 * t23 + t100 * t5 - t11 * t22 - t2 * t197 + t21 * t261) * m(7) + t261 * t44 + t262 * t227 - t170 * mrSges(6,3) - t171 * mrSges(7,2) + t5 * t166 + t59 * t168 + t142 * t132 + (qJ(4) * t59 - t17 * t24 - t18 * t25 + t4 * t197 + t3 * t198 + t259 * t67) * m(6) + t260 * t240 + ((-t196 - t202) * t131 + (t257 - t195) * t129) * t109 + t157 * t238 + t161 * t239 + t12 * t229 + t13 * t230; (m(5) * t98 - t184) * qJD(3) + (-mrSges(5,1) * t191 + (t96 - t141) * t131) * qJD(1) - m(5) * (-t186 * t79 - t199) + (t115 * t152 + t244) * m(7) + (-t115 * t151 + t245) * m(6) - t264; (t217 + t225) * t18 + (-t218 - t226) * t17 - t263 * t174 + t52 - t50 + t51 + t49 - t67 * (mrSges(6,1) * t93 - mrSges(6,2) * t92) - t21 * (mrSges(7,1) * t93 + mrSges(7,3) * t92) + qJD(6) * t63 - pkin(5) * t38 - t43 * t44 + qJ(6) * t36 + t139 + (t10 * t92 + t11 * t93) * mrSges(7,2) + (Ifges(7,3) * t93 - t223) * t237 + t29 * t234 + (t252 * t93 - t253 * t92) * t232 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t18 + t11 * t249 - t21 * t43) * m(7) + (-Ifges(6,2) * t93 + t251 - t88) * t236 + (-t268 * t92 - t224 + t26 + t87) * t235; -t115 * t63 + t93 * t44 + 0.2e1 * (t2 / 0.2e1 + t11 * t232 + t21 * t234) * m(7) + t38;];
tauc  = t14(:);
