% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRRP1
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:27
% EndTime: 2019-03-08 19:55:36
% DurationCPUTime: 5.80s
% Computational Cost: add. (3443->424), mult. (8892->563), div. (0->0), fcn. (6188->10), ass. (0->201)
t273 = Ifges(6,4) + Ifges(7,4);
t139 = cos(qJ(4));
t190 = qJD(2) * t139;
t130 = Ifges(5,4) * t190;
t275 = -t130 / 0.2e1;
t274 = Ifges(6,1) + Ifges(7,1);
t263 = Ifges(6,5) + Ifges(7,5);
t272 = Ifges(6,2) + Ifges(7,2);
t262 = Ifges(6,6) + Ifges(7,6);
t136 = sin(qJ(4));
t134 = cos(pkin(6));
t122 = qJD(1) * t134 + qJD(3);
t199 = t122 * t136;
t140 = cos(qJ(2));
t132 = sin(pkin(6));
t192 = qJD(1) * t132;
t174 = t140 * t192;
t116 = qJD(2) * pkin(2) + t174;
t131 = sin(pkin(11));
t133 = cos(pkin(11));
t137 = sin(qJ(2));
t175 = t137 * t192;
t72 = t131 * t116 + t133 * t175;
t69 = qJD(2) * pkin(8) + t72;
t43 = t139 * t69 + t199;
t194 = t43 * qJD(4);
t149 = t131 * t137 - t133 * t140;
t245 = qJD(2) * t132;
t88 = t149 * t245;
t79 = qJD(1) * t88;
t18 = -t136 * t79 + t194;
t271 = m(5) * (-t18 + t194);
t270 = qJD(4) / 0.2e1;
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t188 = qJD(4) * t138;
t191 = qJD(2) * t136;
t110 = -t135 * t191 + t188;
t269 = t273 * t110;
t189 = qJD(4) * t135;
t111 = t138 * t191 + t189;
t268 = t273 * t111;
t267 = t273 * t138;
t266 = t273 * t135;
t173 = Ifges(5,5) * t270;
t40 = qJD(4) * pkin(9) + t43;
t148 = -pkin(4) * t139 - pkin(9) * t136 - pkin(3);
t117 = t131 * t175;
t71 = t116 * t133 - t117;
t48 = qJD(2) * t148 - t71;
t11 = -t135 * t40 + t138 * t48;
t12 = t135 * t48 + t138 * t40;
t151 = t11 * t138 + t12 * t135;
t162 = mrSges(7,1) * t135 + mrSges(7,2) * t138;
t164 = mrSges(6,1) * t135 + mrSges(6,2) * t138;
t42 = t122 * t139 - t136 * t69;
t39 = -qJD(4) * pkin(4) - t42;
t22 = -pkin(5) * t110 + qJD(6) + t39;
t226 = t138 / 0.2e1;
t229 = -t135 / 0.2e1;
t232 = t111 / 0.2e1;
t125 = qJD(5) - t190;
t247 = t274 * t111 + t263 * t125 + t269;
t248 = t272 * t110 + t262 * t125 + t268;
t249 = t125 / 0.2e1;
t250 = t110 / 0.2e1;
t252 = t274 * t138 - t266;
t254 = -t272 * t135 + t267;
t8 = -qJ(6) * t111 + t11;
t7 = pkin(5) * t125 + t8;
t9 = qJ(6) * t110 + t12;
t243 = t151 * mrSges(6,3) + (t9 * t135 + t7 * t138) * mrSges(7,3) - t22 * t162 - t39 * t164 - t254 * t250 - t252 * t232 - (-t135 * t262 + t138 * t263) * t249 - t248 * t229 - t247 * t226;
t264 = t191 / 0.2e1;
t68 = -qJD(2) * pkin(3) - t71;
t265 = -t68 * mrSges(5,2) + t42 * mrSges(5,3) - Ifges(5,1) * t264 - t173 + t243 + t275;
t183 = qJD(2) * qJD(4);
t171 = t136 * t183;
t182 = qJD(4) * qJD(5);
t186 = qJD(5) * t135;
t187 = qJD(4) * t139;
t80 = t138 * t182 + (-t136 * t186 + t138 * t187) * qJD(2);
t185 = qJD(5) * t138;
t146 = t135 * t187 + t136 * t185;
t81 = -qJD(2) * t146 - t135 * t182;
t261 = t262 * t171 + t272 * t81 + t273 * t80;
t260 = t263 * t171 + t273 * t81 + t274 * t80;
t216 = -qJ(6) - pkin(9);
t170 = qJD(5) * t216;
t184 = qJD(6) * t138;
t169 = pkin(4) * t136 - pkin(9) * t139;
t114 = t169 * qJD(2);
t21 = t135 * t114 + t138 * t42;
t259 = t184 - t21 + (qJ(6) * t190 + t170) * t135;
t195 = t138 * t139;
t147 = pkin(5) * t136 - qJ(6) * t195;
t20 = t138 * t114 - t135 * t42;
t258 = -qJD(2) * t147 - qJD(6) * t135 + t138 * t170 - t20;
t256 = t42 * qJD(4);
t255 = t272 * t138 + t266;
t253 = t274 * t135 + t267;
t251 = -m(5) * t42 + m(6) * t39;
t172 = -Ifges(5,6) * qJD(4) / 0.2e1;
t246 = t263 * t135 + t262 * t138;
t17 = -t139 * t79 + t256;
t115 = t169 * qJD(4);
t91 = (t131 * t140 + t133 * t137) * t132;
t86 = qJD(1) * t91;
t49 = (t115 + t86) * qJD(2);
t3 = t135 * t49 + t138 * t17 + t48 * t185 - t186 * t40;
t4 = -qJD(5) * t12 - t135 * t17 + t138 * t49;
t167 = -t135 * t4 + t138 * t3;
t178 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t179 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t180 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t242 = -t179 * t110 - t180 * t111 - t178 * t125 - t11 * mrSges(6,1) - t68 * mrSges(5,1) - t7 * mrSges(7,1) - t172 + (Ifges(5,4) * t136 + t139 * Ifges(5,2)) * qJD(2) / 0.2e1 + t12 * mrSges(6,2) + t43 * mrSges(5,3) + t9 * mrSges(7,2) - t262 * t250 - (Ifges(7,3) + Ifges(6,3)) * t249 - t263 * t232;
t241 = 0.2e1 * m(5);
t240 = t80 / 0.2e1;
t239 = t81 / 0.2e1;
t236 = m(7) * t22;
t235 = -t110 / 0.2e1;
t233 = -t111 / 0.2e1;
t231 = -t125 / 0.2e1;
t225 = pkin(2) * t133;
t224 = pkin(5) * t135;
t150 = t134 * t139 - t136 * t91;
t221 = t18 * t150;
t87 = qJD(2) * t91;
t78 = qJD(1) * t87;
t90 = t149 * t132;
t218 = t78 * t90;
t35 = -t81 * mrSges(7,1) + t80 * mrSges(7,2);
t36 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t215 = t35 + t36;
t57 = mrSges(7,1) * t171 - mrSges(7,3) * t80;
t58 = mrSges(6,1) * t171 - mrSges(6,3) * t80;
t214 = t57 + t58;
t59 = -mrSges(7,2) * t171 + mrSges(7,3) * t81;
t60 = -mrSges(6,2) * t171 + mrSges(6,3) * t81;
t213 = t59 + t60;
t82 = -mrSges(7,2) * t125 + mrSges(7,3) * t110;
t83 = -mrSges(6,2) * t125 + mrSges(6,3) * t110;
t212 = t82 + t83;
t84 = mrSges(7,1) * t125 - mrSges(7,3) * t111;
t85 = mrSges(6,1) * t125 - mrSges(6,3) * t111;
t211 = t84 + t85;
t106 = t148 - t225;
t204 = t106 * t185 + t135 * t115;
t127 = pkin(2) * t131 + pkin(8);
t113 = t127 * t195;
t63 = t135 * t106 + t113;
t177 = mrSges(5,3) * t191;
t203 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t110 + mrSges(6,2) * t111 + t177;
t200 = qJ(6) * t136;
t197 = t135 * t139;
t196 = t136 * t138;
t193 = t136 * t127 * t189 + t138 * t115;
t64 = -mrSges(7,1) * t110 + mrSges(7,2) * t111;
t181 = -t64 - t203;
t176 = mrSges(5,3) * t190;
t1 = pkin(5) * t171 - qJ(6) * t80 - qJD(6) * t111 + t4;
t2 = qJ(6) * t81 + qJD(6) * t110 + t3;
t168 = -t1 * t135 + t138 * t2;
t165 = mrSges(6,1) * t138 - mrSges(6,2) * t135;
t163 = mrSges(7,1) * t138 - t135 * mrSges(7,2);
t67 = t134 * t136 + t139 * t91;
t32 = t135 * t90 + t138 * t67;
t31 = -t135 * t67 + t138 * t90;
t145 = -t4 * mrSges(6,1) - t1 * mrSges(7,1) + t3 * mrSges(6,2) + t2 * mrSges(7,2);
t129 = -pkin(5) * t138 - pkin(4);
t128 = -pkin(3) - t225;
t124 = Ifges(6,3) * t171;
t123 = Ifges(7,3) * t171;
t121 = t216 * t138;
t120 = t216 * t135;
t119 = -qJD(4) * mrSges(5,2) + t176;
t112 = (-mrSges(5,1) * t139 + mrSges(5,2) * t136) * qJD(2);
t104 = (mrSges(5,1) * t136 + mrSges(5,2) * t139) * t183;
t97 = (t127 + t224) * t136;
t94 = t138 * t106;
t89 = t133 * t174 - t117;
t77 = Ifges(6,5) * t80;
t76 = Ifges(7,5) * t80;
t75 = Ifges(6,6) * t81;
t74 = Ifges(7,6) * t81;
t70 = pkin(5) * t146 + t127 * t187;
t62 = -t127 * t197 + t94;
t56 = -t135 * t200 + t63;
t47 = -qJ(6) * t196 + t94 + (-t127 * t135 - pkin(5)) * t139;
t37 = t199 + (qJD(2) * t224 + t69) * t139;
t34 = t135 * t86 + t195 * t89;
t33 = t138 * t86 - t197 * t89;
t30 = qJD(4) * t150 - t139 * t88;
t29 = qJD(4) * t67 - t136 * t88;
t24 = -qJD(5) * t63 + t193;
t23 = (-t136 * t188 - t139 * t186) * t127 + t204;
t14 = (-qJ(6) * qJD(5) - qJD(4) * t127) * t196 + (-qJD(6) * t136 + (-qJ(6) * qJD(4) - qJD(5) * t127) * t139) * t135 + t204;
t13 = -t136 * t184 + t147 * qJD(4) + (-t113 + (-t106 + t200) * t135) * qJD(5) + t193;
t10 = -pkin(5) * t81 + t18;
t6 = qJD(5) * t31 + t135 * t87 + t138 * t30;
t5 = -qJD(5) * t32 - t135 * t30 + t138 * t87;
t15 = [t90 * t104 + t87 * t112 + t30 * t119 - t215 * t150 + t212 * t6 + t211 * t5 + t213 * t32 + t214 * t31 - t181 * t29 + m(5) * (t17 * t67 - t29 * t42 + t30 * t43 + t68 * t87 + t218 - t221) + m(4) * (-t71 * t87 - t72 * t88 - t79 * t91 + t218) + m(7) * (t1 * t31 - t10 * t150 + t2 * t32 + t22 * t29 + t5 * t7 + t6 * t9) + m(6) * (t11 * t5 + t12 * t6 + t29 * t39 + t3 * t32 + t31 * t4 - t221) + (-t87 * mrSges(4,1) + t88 * mrSges(4,2) + (-t136 * t67 - t139 * t150) * qJD(4) * mrSges(5,3) + (-mrSges(3,1) * t137 - mrSges(3,2) * t140) * t245) * qJD(2); (qJD(2) * t89 + t79) * mrSges(4,2) + t128 * t104 - t86 * t112 + t97 * t35 + t14 * t82 + t23 * t83 + t13 * t84 + t24 * t85 + t70 * t64 + t47 * t57 + t56 * t59 + t62 * t58 + t63 * t60 - m(6) * (t11 * t33 + t12 * t34) + (qJD(2) * t86 - t78) * mrSges(4,1) + (t78 * mrSges(5,2) + t10 * t162 + (mrSges(5,3) + t164) * t18 + (-t1 * t138 - t2 * t135) * mrSges(7,3) + (-t3 * t135 - t4 * t138) * mrSges(6,3) + (((-0.3e1 / 0.2e1 * Ifges(5,4) + t180 * t138 - t179 * t135) * t136 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - t178) * t139) * qJD(2) + t172 - t242) * qJD(4) + (t39 * t165 + t22 * t163 + (t7 * t135 - t9 * t138) * mrSges(7,3) + (t11 * t135 - t12 * t138) * mrSges(6,3) + t255 * t235 + t253 * t233 + t246 * t231 - t248 * t138 / 0.2e1) * qJD(5) + t252 * t240 + t254 * t239 + (qJD(5) * t247 + t261) * t229 + t260 * t226 + (m(6) * t18 - qJD(4) * t119 - t271 + t36) * t127 + (t181 - t236 - t251) * t89) * t136 + m(6) * (t11 * t24 + t12 * t23 + t3 * t63 + t4 * t62) + (-t68 * t86 / 0.2e1 + t128 * t78 / 0.2e1) * t241 - t211 * t33 - t212 * t34 + (-t123 / 0.2e1 - t124 / 0.2e1 - t77 / 0.2e1 - t75 / 0.2e1 - t76 / 0.2e1 - t74 / 0.2e1 - t78 * mrSges(5,1) + t17 * mrSges(5,3) - t89 * t119 - t179 * t81 - t180 * t80 + (t127 * t17 / 0.2e1 - t43 * t89 / 0.2e1) * t241 + ((t203 + t251) * t127 + t173 + 0.3e1 / 0.2e1 * t130 - t265) * qJD(4) + t145) * t139 + (t1 * t47 + t10 * t97 + t2 * t56 + t22 * t70 + (-t34 + t14) * t9 + (-t33 + t13) * t7) * m(7) + (t71 * t86 - t72 * t89 + (-t131 * t79 - t133 * t78) * pkin(2)) * m(4); ((-t135 * t211 + t138 * t212 + t119 - t176) * qJD(4) + m(7) * (t188 * t9 - t189 * t7 - t10) + t271 + m(6) * (-t11 * t189 + t12 * t188 - t18) - t215) * t139 + (t213 * t138 - t214 * t135 + (-t135 * t212 - t138 * t211) * qJD(5) + (-t177 - t181) * qJD(4) + m(7) * (qJD(4) * t22 - t185 * t7 - t186 * t9 + t168) + m(5) * (t17 - t256) + m(6) * (qJD(4) * t39 - t11 * t185 - t12 * t186 + t167)) * t136; ((-m(6) * t151 - t135 * t83 - t138 * t85) * pkin(9) - t243 + (t64 + t236) * t224) * qJD(5) + t167 * mrSges(6,3) + t168 * mrSges(7,3) - t10 * t163 + (-mrSges(5,1) - t165) * t18 + t129 * t35 - t42 * t119 + t120 * t57 - t121 * t59 + m(6) * (-pkin(4) * t18 + pkin(9) * t167) - t21 * t83 - t20 * t85 - t37 * t64 - pkin(4) * t36 - t17 * mrSges(5,2) - m(6) * (t11 * t20 + t12 * t21 + t39 * t43) + (-t135 * t58 + t138 * t60) * pkin(9) + t259 * t82 + (t1 * t120 + t10 * t129 - t121 * t2 - t22 * t37 + t258 * t7 + t259 * t9) * m(7) + t260 * t135 / 0.2e1 + t261 * t226 - t203 * t43 + ((Ifges(5,4) * t264 + t246 * t270 + t172 + t242) * t136 + (t275 + t173 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t191 + t265) * t139) * qJD(2) + t258 * t84 + t253 * t240 + t255 * t239; t123 + t124 - t39 * (mrSges(6,1) * t111 + mrSges(6,2) * t110) - t22 * (mrSges(7,1) * t111 + mrSges(7,2) * t110) + t77 + t75 + t76 + t74 - t8 * t82 - t11 * t83 + t9 * t84 + t12 * t85 - t145 + (-t111 * t64 + t57) * pkin(5) + (t110 * t7 + t111 * t9) * mrSges(7,3) + (t11 * t110 + t111 * t12) * mrSges(6,3) + (-(-t7 + t8) * t9 + (-t111 * t22 + t1) * pkin(5)) * m(7) + (t274 * t110 - t268) * t233 + t248 * t232 + (t110 * t263 - t111 * t262) * t231 + (-t272 * t111 + t247 + t269) * t235; -t110 * t82 + t111 * t84 + 0.2e1 * (t10 / 0.2e1 + t9 * t235 + t7 * t232) * m(7) + t35;];
tauc  = t15(:);
