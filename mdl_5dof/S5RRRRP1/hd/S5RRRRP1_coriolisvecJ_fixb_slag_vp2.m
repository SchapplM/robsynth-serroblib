% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:46
% EndTime: 2019-12-05 18:45:01
% DurationCPUTime: 5.64s
% Computational Cost: add. (5504->372), mult. (14716->499), div. (0->0), fcn. (10321->6), ass. (0->177)
t191 = sin(qJ(3));
t192 = sin(qJ(2));
t194 = cos(qJ(3));
t195 = cos(qJ(2));
t166 = -t191 * t192 + t194 * t195;
t154 = t166 * qJD(1);
t167 = t191 * t195 + t194 * t192;
t155 = t167 * qJD(1);
t190 = sin(qJ(4));
t193 = cos(qJ(4));
t116 = t154 * t190 + t155 * t193;
t184 = -pkin(2) * t195 - pkin(1);
t177 = qJD(1) * t184;
t134 = -t154 * pkin(3) + t177;
t189 = qJD(2) + qJD(3);
t188 = qJD(4) + t189;
t199 = t193 * t154 - t155 * t190;
t132 = t189 * t166;
t122 = t132 * qJD(1);
t133 = t189 * t167;
t123 = t133 * qJD(1);
t38 = -qJD(4) * t116 - t122 * t190 - t123 * t193;
t204 = qJD(4) * t193;
t205 = qJD(4) * t190;
t249 = -pkin(7) - pkin(6);
t178 = t249 * t192;
t171 = qJD(1) * t178;
t162 = qJD(2) * pkin(2) + t171;
t203 = qJD(2) * t249;
t198 = qJD(1) * t203;
t163 = t192 * t198;
t164 = t195 * t198;
t179 = t249 * t195;
t172 = qJD(1) * t179;
t206 = qJD(3) * t194;
t207 = qJD(3) * t191;
t70 = t162 * t206 + t194 * t163 + t191 * t164 + t172 * t207;
t44 = -pkin(8) * t123 + t70;
t159 = t194 * t172;
t127 = t162 * t191 - t159;
t71 = -qJD(3) * t127 - t163 * t191 + t194 * t164;
t45 = -pkin(8) * t122 + t71;
t156 = t191 * t172;
t126 = t194 * t162 + t156;
t149 = t155 * pkin(8);
t92 = t126 - t149;
t81 = pkin(3) * t189 + t92;
t232 = pkin(8) * t154;
t93 = t127 + t232;
t7 = t190 * t45 + t193 * t44 + t81 * t204 - t205 * t93;
t2 = qJ(5) * t38 + qJD(5) * t199 + t7;
t241 = t116 / 0.2e1;
t260 = Ifges(5,5) + Ifges(6,5);
t261 = Ifges(5,1) + Ifges(6,1);
t231 = Ifges(5,4) + Ifges(6,4);
t264 = t231 * t199;
t256 = t261 * t116 + t260 * t188 + t264;
t258 = Ifges(5,6) + Ifges(6,6);
t259 = Ifges(5,2) + Ifges(6,2);
t266 = t231 * t116;
t257 = t258 * t188 + t259 * t199 + t266;
t262 = -t199 / 0.2e1;
t37 = t199 * qJD(4) + t122 * t193 - t123 * t190;
t84 = t193 * t93;
t40 = t190 * t81 + t84;
t8 = -qJD(4) * t40 - t190 * t44 + t193 * t45;
t3 = -qJ(5) * t37 - qJD(5) * t116 + t8;
t72 = -pkin(4) * t199 + qJD(5) + t134;
t268 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t7 * mrSges(5,2) - t2 * mrSges(6,2) + t258 * t38 + t260 * t37 - t134 * (mrSges(5,1) * t116 + mrSges(5,2) * t199) - t72 * (mrSges(6,1) * t116 + mrSges(6,2) * t199) + t257 * t241 - (-t258 * t116 + t260 * t199) * t188 / 0.2e1 + (-t259 * t116 + t256 + t264) * t262 - (t261 * t199 - t266) * t116 / 0.2e1;
t233 = pkin(4) * t116;
t216 = qJ(5) * t199;
t15 = t40 + t216;
t226 = t116 * t15;
t225 = t116 * t40;
t108 = qJ(5) * t116;
t82 = t190 * t93;
t39 = t193 * t81 - t82;
t14 = t39 - t108;
t13 = pkin(4) * t188 + t14;
t253 = t13 * t199;
t252 = t199 * t39;
t135 = t194 * t178 + t179 * t191;
t111 = -pkin(8) * t167 + t135;
t136 = t191 * t178 - t194 * t179;
t112 = pkin(8) * t166 + t136;
t64 = t190 * t111 + t193 * t112;
t248 = pkin(1) * mrSges(3,1);
t247 = pkin(1) * mrSges(3,2);
t239 = t154 / 0.2e1;
t238 = t155 / 0.2e1;
t235 = m(4) * t177;
t234 = pkin(3) * t155;
t47 = t193 * t92 - t82;
t130 = -t171 * t191 + t159;
t94 = t130 - t232;
t131 = t194 * t171 + t156;
t95 = -t149 + t131;
t49 = t190 * t94 + t193 * t95;
t96 = -mrSges(6,2) * t188 + mrSges(6,3) * t199;
t97 = -mrSges(5,2) * t188 + mrSges(5,3) * t199;
t230 = t96 + t97;
t98 = mrSges(6,1) * t188 - mrSges(6,3) * t116;
t99 = mrSges(5,1) * t188 - mrSges(5,3) * t116;
t229 = t98 + t99;
t228 = mrSges(4,3) * t154;
t227 = Ifges(3,4) * t192;
t183 = pkin(2) * t194 + pkin(3);
t211 = t191 * t193;
t151 = pkin(2) * t211 + t183 * t190;
t222 = t151 * t38;
t221 = t155 * mrSges(4,3);
t220 = t155 * Ifges(4,4);
t219 = t190 * t38;
t218 = Ifges(3,5) * qJD(2);
t217 = Ifges(3,6) * qJD(2);
t215 = qJD(2) * mrSges(3,1);
t214 = qJD(2) * mrSges(3,2);
t212 = t190 * t191;
t210 = qJD(1) * t192;
t209 = qJD(1) * t195;
t208 = qJD(2) * t192;
t186 = pkin(2) * t210;
t202 = t218 / 0.2e1;
t201 = -t217 / 0.2e1;
t200 = -t38 * mrSges(6,1) + t37 * mrSges(6,2);
t100 = pkin(3) * t123 + qJD(2) * t186;
t117 = pkin(2) * t208 + pkin(3) * t133;
t46 = -t190 * t92 - t84;
t48 = -t190 * t95 + t193 * t94;
t63 = t193 * t111 - t112 * t190;
t150 = -pkin(2) * t212 + t193 * t183;
t80 = t234 + t233;
t128 = t166 * t193 - t167 * t190;
t129 = t166 * t190 + t167 * t193;
t140 = -t166 * pkin(3) + t184;
t173 = t192 * t203;
t174 = t195 * t203;
t76 = t194 * t173 + t191 * t174 + t178 * t206 + t179 * t207;
t55 = -pkin(8) * t133 + t76;
t77 = -t136 * qJD(3) - t173 * t191 + t194 * t174;
t56 = -pkin(8) * t132 + t77;
t9 = t111 * t204 - t112 * t205 + t190 * t56 + t193 * t55;
t10 = -t64 * qJD(4) - t190 * t55 + t193 * t56;
t106 = t154 * Ifges(4,2) + t189 * Ifges(4,6) + t220;
t148 = Ifges(4,4) * t154;
t107 = t155 * Ifges(4,1) + t189 * Ifges(4,5) + t148;
t196 = -t177 * (mrSges(4,1) * t155 + mrSges(4,2) * t154) + mrSges(5,3) * t252 + mrSges(6,3) * t253 - t155 * (Ifges(4,1) * t154 - t220) / 0.2e1 - t189 * (Ifges(4,5) * t154 - Ifges(4,6) * t155) / 0.2e1 - Ifges(4,6) * t123 + t106 * t238 + t126 * t228 - t70 * mrSges(4,2) + t71 * mrSges(4,1) + Ifges(4,5) * t122 - (-Ifges(4,2) * t155 + t107 + t148) * t154 / 0.2e1 + t268;
t185 = Ifges(3,4) * t209;
t182 = pkin(3) * t193 + pkin(4);
t176 = mrSges(3,3) * t209 - t214;
t175 = -mrSges(3,3) * t210 + t215;
t153 = Ifges(3,1) * t210 + t185 + t218;
t152 = t217 + (t195 * Ifges(3,2) + t227) * qJD(1);
t147 = pkin(4) + t150;
t139 = mrSges(4,1) * t189 - t221;
t138 = -mrSges(4,2) * t189 + t228;
t137 = t186 + t234;
t125 = -mrSges(4,1) * t154 + mrSges(4,2) * t155;
t121 = -t183 * t205 + (-t191 * t204 + (-t190 * t194 - t211) * qJD(3)) * pkin(2);
t120 = t183 * t204 + (-t191 * t205 + (t193 * t194 - t212) * qJD(3)) * pkin(2);
t85 = -t128 * pkin(4) + t140;
t73 = t186 + t80;
t68 = -mrSges(5,1) * t199 + mrSges(5,2) * t116;
t67 = -mrSges(6,1) * t199 + mrSges(6,2) * t116;
t51 = -qJD(4) * t129 - t132 * t190 - t133 * t193;
t50 = qJD(4) * t128 + t132 * t193 - t133 * t190;
t27 = qJ(5) * t128 + t64;
t26 = -qJ(5) * t129 + t63;
t21 = -pkin(4) * t51 + t117;
t20 = -t108 + t49;
t19 = t48 - t216;
t18 = -t108 + t47;
t17 = t46 - t216;
t12 = -pkin(4) * t38 + t100;
t5 = -qJ(5) * t50 - qJD(5) * t129 + t10;
t4 = qJ(5) * t51 + qJD(5) * t128 + t9;
t1 = [(-t122 * t135 - t123 * t136 - t126 * t132 - t127 * t133 + t166 * t70 - t167 * t71) * mrSges(4,3) + (-t123 * t166 - t133 * t239) * Ifges(4,2) + (t122 * t166 - t123 * t167 + t132 * t239 - t133 * t238) * Ifges(4,4) + t184 * (mrSges(4,1) * t123 + mrSges(4,2) * t122) + (-mrSges(5,1) * t100 - mrSges(6,1) * t12 + mrSges(5,3) * t7 + mrSges(6,3) * t2 + t231 * t37 + t259 * t38) * t128 + (t258 * t51 + t260 * t50) * t188 / 0.2e1 + (t231 * t50 + t259 * t51) * t199 / 0.2e1 + (t231 * t51 + t261 * t50) * t241 + (mrSges(5,2) * t100 + mrSges(6,2) * t12 - mrSges(5,3) * t8 - mrSges(6,3) * t3 + t231 * t38 + t261 * t37) * t129 + t256 * t50 / 0.2e1 + t257 * t51 / 0.2e1 + (-t13 * t50 + t15 * t51 - t26 * t37 + t27 * t38) * mrSges(6,3) + (-t37 * t63 + t38 * t64 - t39 * t50 + t40 * t51) * mrSges(5,3) + m(6) * (t12 * t85 + t13 * t5 + t15 * t4 + t2 * t27 + t21 * t72 + t26 * t3) + m(5) * (t10 * t39 + t100 * t140 + t117 * t134 + t40 * t9 + t63 * t8 + t64 * t7) + t189 * (Ifges(4,5) * t132 - Ifges(4,6) * t133) / 0.2e1 + t177 * (mrSges(4,1) * t133 + mrSges(4,2) * t132) + t85 * t200 + t132 * t107 / 0.2e1 - t133 * t106 / 0.2e1 + t134 * (-mrSges(5,1) * t51 + mrSges(5,2) * t50) + t76 * t138 + t77 * t139 + t140 * (-mrSges(5,1) * t38 + mrSges(5,2) * t37) + (-pkin(6) * t175 + t153 / 0.2e1 + t202 + (-0.2e1 * t247 + 0.3e1 / 0.2e1 * Ifges(3,4) * t195) * qJD(1)) * t195 * qJD(2) + (t122 * t167 + t132 * t238) * Ifges(4,1) + (-pkin(6) * t176 - t152 / 0.2e1 + t201 + (-0.2e1 * t248 - 0.3e1 / 0.2e1 * t227 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t195) * qJD(1) + (t125 + qJD(1) * (-mrSges(4,1) * t166 + mrSges(4,2) * t167) + 0.2e1 * t235) * pkin(2)) * t208 + m(4) * (t126 * t77 + t127 * t76 + t135 * t71 + t136 * t70) + t21 * t67 + t72 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + t4 * t96 + t9 * t97 + t5 * t98 + t10 * t99 + t117 * t68; (-t147 * t37 + t222 + t226) * mrSges(6,3) + (-t150 * t37 + t222 + t225) * mrSges(5,3) + ((t202 - t153 / 0.2e1 - t185 / 0.2e1 + qJD(1) * t247 + (t175 - t215) * pkin(6)) * t195 + (t201 + t152 / 0.2e1 + (t248 + t227 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t195) * qJD(1) + (t176 + t214) * pkin(6) + (-t125 - t235) * pkin(2)) * t192) * qJD(1) + t196 + t229 * t121 + t230 * t120 - t137 * t68 - t131 * t138 - t130 * t139 + t127 * t221 + ((t138 * t194 - t139 * t191) * qJD(3) + (-t122 * t194 - t123 * t191) * mrSges(4,3)) * pkin(2) - t73 * t67 - t20 * t96 - t49 * t97 - t19 * t98 - t48 * t99 + (t147 * t3 + t151 * t2 - t72 * t73 + (-t20 + t120) * t15 + (-t19 + t121) * t13) * m(6) + (-t134 * t137 + t150 * t8 + t151 * t7 + (-t49 + t120) * t40 + (t121 - t48) * t39) * m(5) + ((t191 * t70 + t194 * t71 + (-t126 * t191 + t127 * t194) * qJD(3)) * pkin(2) - t126 * t130 - t127 * t131) * m(4); -m(5) * (t39 * t46 + t40 * t47) + (mrSges(6,3) * t219 - t155 * t68 + (-t193 * t37 + t219) * mrSges(5,3) + (-t190 * t229 + t193 * t230) * qJD(4) + (-t134 * t155 + t190 * t7 + t193 * t8 + t204 * t40 - t205 * t39) * m(5)) * pkin(3) + (-t182 * t37 + t226) * mrSges(6,3) + t196 + (t139 + t221) * t127 - t126 * t138 + mrSges(5,3) * t225 - t80 * t67 - t18 * t96 - t47 * t97 - t17 * t98 - t46 * t99 + ((-t13 * t205 + t15 * t204 + t190 * t2) * pkin(3) - t13 * t17 - t15 * t18 - t72 * t80 + t182 * t3) * m(6); (t3 * pkin(4) - t72 * t233 - (-t13 + t14) * t15) * m(6) + (-pkin(4) * t37 + t226 + t253) * mrSges(6,3) + (t252 + t225) * mrSges(5,3) - t67 * t233 - t14 * t96 - t39 * t97 + t15 * t98 + t40 * t99 + t268; -t199 * t96 + t116 * t98 + 0.2e1 * (t12 / 0.2e1 + t15 * t262 + t13 * t241) * m(6) + t200;];
tauc = t1(:);
