% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:30
% EndTime: 2020-01-03 12:05:36
% DurationCPUTime: 2.56s
% Computational Cost: add. (3996->298), mult. (6931->435), div. (0->0), fcn. (4319->8), ass. (0->186)
t163 = sin(pkin(9));
t164 = cos(pkin(9));
t145 = -t164 * pkin(3) - t163 * pkin(7) - pkin(2);
t166 = sin(qJ(4));
t169 = cos(qJ(4));
t232 = t164 * t169;
t111 = qJ(3) * t232 + t166 * t145;
t167 = sin(qJ(2));
t225 = qJD(3) * t164;
t170 = cos(qJ(2));
t231 = t164 * t170;
t247 = qJD(1) * pkin(1);
t264 = -t111 * qJD(4) - t166 * t225 - (-t166 * t231 + t167 * t169) * t247;
t223 = qJD(4) * t169;
t273 = (t166 * t167 + t169 * t231) * t247 - t145 * t223 - t169 * t225;
t224 = qJD(4) * t166;
t203 = t163 * t224;
t150 = pkin(8) * t203;
t272 = t150 + t264;
t233 = t164 * t166;
t211 = qJ(3) * t233;
t234 = t163 * t169;
t218 = pkin(8) * t234;
t271 = -(-t211 - t218) * qJD(4) + t273;
t162 = qJD(1) + qJD(2);
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t140 = t165 * t169 + t166 * t168;
t259 = qJD(4) + qJD(5);
t91 = t259 * t140;
t67 = t91 * t163;
t56 = t162 * t67;
t270 = Ifges(6,5) * t56;
t183 = t165 * t166 - t168 * t169;
t120 = t183 * t163;
t68 = t259 * t120;
t57 = t162 * t68;
t269 = Ifges(6,6) * t57;
t235 = t163 * t166;
t207 = t162 * t235;
t143 = qJ(3) * t162 + t167 * t247;
t209 = t143 * t232;
t217 = t170 * t247;
t191 = qJD(3) - t217;
t97 = t145 * t162 + t191;
t61 = t166 * t97 + t209;
t50 = -pkin(8) * t207 + t61;
t245 = t165 * t50;
t236 = t162 * t164;
t149 = qJD(4) - t236;
t206 = t162 * t234;
t210 = t143 * t233;
t175 = -pkin(8) * t206 - t210;
t86 = t169 * t97;
t49 = t175 + t86;
t40 = pkin(4) * t149 + t49;
t12 = t168 * t40 - t245;
t212 = qJD(2) * t247;
t196 = t170 * t212;
t220 = t162 * qJD(3);
t136 = t196 + t220;
t195 = t167 * t212;
t213 = t136 * t232 + t166 * t195 + t97 * t223;
t32 = qJD(4) * t175 + t213;
t190 = -t136 * t233 + t169 * t195;
t237 = t162 * t163;
t33 = (-t209 + (pkin(8) * t237 - t97) * t166) * qJD(4) + t190;
t4 = qJD(5) * t12 + t165 * t33 + t168 * t32;
t244 = t168 * t50;
t13 = t165 * t40 + t244;
t5 = -qJD(5) * t13 - t165 * t32 + t168 * t33;
t268 = -t5 * mrSges(6,1) + t4 * mrSges(6,2);
t201 = t164 * t224;
t36 = -t143 * t201 + t213;
t37 = -qJD(4) * t61 + t190;
t267 = -t37 * mrSges(5,1) + t36 * mrSges(5,2);
t135 = t169 * t145;
t74 = -t218 + t135 + (-qJ(3) * t166 - pkin(4)) * t164;
t219 = pkin(8) * t235;
t94 = -t219 + t111;
t43 = t165 * t74 + t168 * t94;
t266 = -qJD(5) * t43 + t271 * t165 + t272 * t168;
t42 = -t165 * t94 + t168 * t74;
t265 = qJD(5) * t42 + t272 * t165 - t271 * t168;
t160 = t163 ^ 2;
t161 = t164 ^ 2;
t198 = (t160 + t161) * mrSges(4,3);
t263 = -qJ(3) * t201 - t273;
t262 = (t236 - t259) * t183;
t261 = t140 * t236 - t91;
t242 = t136 * t160;
t241 = t143 * t160;
t256 = pkin(1) * t170;
t131 = t145 - t256;
t158 = pkin(1) * t167 + qJ(3);
t88 = t166 * t131 + t158 * t232;
t226 = qJD(3) * t143;
t260 = qJ(3) * t242 + t160 * t226 - t217 * t241;
t102 = -t165 * t207 + t168 * t206;
t257 = t102 / 0.2e1;
t255 = pkin(4) * t162;
t254 = t12 * t67;
t253 = t269 - t270;
t119 = t140 * t163;
t100 = t162 * t119;
t252 = mrSges(6,3) * t100;
t251 = Ifges(5,4) * t166;
t250 = Ifges(5,4) * t169;
t249 = Ifges(6,4) * t102;
t248 = pkin(1) * qJD(2);
t246 = t102 * mrSges(6,3);
t243 = t169 * t61;
t240 = t143 * t162;
t215 = t170 * t248;
t154 = qJD(3) + t215;
t239 = t154 * t163;
t238 = t158 * t163;
t230 = t154 * t241 + t158 * t242;
t222 = qJD(5) * t165;
t221 = qJD(5) * t168;
t216 = t167 * t248;
t208 = t158 * t233;
t205 = t131 * t223 + t154 * t232 + t166 * t216;
t202 = t163 * t223;
t200 = -qJD(4) * t163 / 0.2e1;
t187 = mrSges(5,1) * t166 + mrSges(5,2) * t169;
t109 = t187 * t237;
t59 = mrSges(6,1) * t100 + mrSges(6,2) * t102;
t199 = (-t109 - t59) * t163;
t197 = mrSges(5,3) * t203;
t194 = mrSges(5,3) * t202;
t193 = t166 * t200;
t192 = t169 * t200;
t189 = -mrSges(4,1) * t164 + mrSges(4,2) * t163;
t188 = mrSges(5,1) * t169 - mrSges(5,2) * t166;
t186 = -Ifges(5,5) * t166 - Ifges(5,6) * t169;
t123 = t169 * t131;
t64 = -t218 + t123 + (-t158 * t166 - pkin(4)) * t164;
t73 = -t219 + t88;
t34 = -t165 * t73 + t168 * t64;
t35 = t165 * t64 + t168 * t73;
t60 = t86 - t210;
t185 = -t166 * t60 + t243;
t107 = -mrSges(5,2) * t149 - mrSges(5,3) * t207;
t108 = mrSges(5,1) * t149 - mrSges(5,3) * t206;
t184 = t169 * t107 - t166 * t108;
t182 = t166 * (-Ifges(5,2) * t169 - t251);
t181 = t169 * (-Ifges(5,1) * t166 - t250);
t180 = (t169 * Ifges(5,1) - t251) * t163;
t179 = (-t166 * Ifges(5,2) + t250) * t163;
t178 = qJD(4) * t188;
t177 = t186 * qJD(4);
t144 = qJD(5) + t149;
t45 = -Ifges(6,2) * t100 + Ifges(6,6) * t144 + t249;
t98 = Ifges(6,4) * t100;
t46 = Ifges(6,1) * t102 + Ifges(6,5) * t144 - t98;
t99 = (t166 * t255 + t143) * t163;
t172 = -t12 * t252 + t45 * t257 - t99 * (mrSges(6,1) * t102 - mrSges(6,2) * t100) + t253 - t102 * (-Ifges(6,1) * t100 - t249) / 0.2e1 - t144 * (-Ifges(6,5) * t100 - Ifges(6,6) * t102) / 0.2e1 + (-Ifges(6,2) * t102 + t46 - t98) * t100 / 0.2e1 - t268;
t55 = -t88 * qJD(4) - t154 * t233 + t169 * t216;
t77 = t149 * Ifges(5,6) + t162 * t179;
t78 = t149 * Ifges(5,5) + t162 * t180;
t92 = (t223 * t255 + t136) * t163;
t171 = t144 * (-Ifges(6,5) * t67 + Ifges(6,6) * t68) / 0.2e1 + t99 * (-mrSges(6,1) * t68 - mrSges(6,2) * t67) - t100 * (-Ifges(6,4) * t67 + Ifges(6,2) * t68) / 0.2e1 - t67 * t46 / 0.2e1 + t68 * t45 / 0.2e1 + t13 * t68 * mrSges(6,3) + t189 * t195 + t60 * t197 + t178 * t241 + t187 * t242 + (-Ifges(6,1) * t67 + Ifges(6,4) * t68) * t257 + t78 * t193 + t77 * t192 + (t149 / 0.2e1 - t236 / 0.2e1) * t163 * t177 + (-t37 * t234 - t36 * t235) * mrSges(5,3) + (t136 * t161 + t242) * mrSges(4,3) + (-t269 / 0.2e1 + t270 / 0.2e1 - t253 / 0.2e1 + t267 + t268) * t164 + (-t92 * mrSges(6,2) + t5 * mrSges(6,3) + Ifges(6,1) * t56 - Ifges(6,4) * t57) * t120 + (t92 * mrSges(6,1) - t4 * mrSges(6,3) + Ifges(6,4) * t56 - Ifges(6,2) * t57) * t119 + ((-t182 + t181) * qJD(4) * t160 + t179 * t192 + t180 * t193 + (-Ifges(5,5) * t193 - Ifges(5,6) * t192) * t164) * t162;
t157 = pkin(4) * t235;
t151 = pkin(4) * t202;
t138 = qJ(3) * t163 + t157;
t137 = -t162 * pkin(2) + t191;
t132 = qJD(3) * t163 + t151;
t127 = t189 * t162;
t126 = t157 + t238;
t118 = t160 * t240;
t112 = t151 + t239;
t110 = t135 - t211;
t104 = t178 * t237;
t87 = t123 - t208;
t72 = mrSges(6,1) * t144 - t246;
t71 = -mrSges(6,2) * t144 - t252;
t54 = -t158 * t201 + t205;
t48 = t150 + t55;
t47 = (-t208 - t218) * qJD(4) + t205;
t19 = -mrSges(6,1) * t57 - mrSges(6,2) * t56;
t17 = t168 * t49 - t245;
t16 = -t165 * t49 - t244;
t7 = -qJD(5) * t35 - t165 * t47 + t168 * t48;
t6 = qJD(5) * t34 + t165 * t48 + t168 * t47;
t1 = [t127 * t216 + t112 * t59 + t126 * t19 + t54 * t107 + t55 * t108 + t6 * t71 + t7 * t72 + m(6) * (t112 * t99 + t12 * t7 + t126 * t92 + t13 * t6 + t34 * t5 + t35 * t4) + t104 * t238 + t109 * t239 + t171 - t61 * t194 - mrSges(3,1) * t195 - mrSges(3,2) * t196 + m(5) * (t36 * t88 + t37 * t87 + t54 * t61 + t55 * t60 + t230) + m(4) * ((t136 * t158 + t143 * t154) * t161 + (t137 + (-pkin(2) - t256) * qJD(1)) * t216 + t230) + (t34 * t56 + t35 * t57 + t254) * mrSges(6,3) + (-mrSges(3,1) * t216 - mrSges(3,2) * t215 + t154 * t198 - t88 * t194 + t87 * t197) * t162; t264 * t108 + t263 * t107 + (t42 * t56 + t43 * t57 + t254) * mrSges(6,3) + t132 * t59 + t138 * t19 + ((-t127 + (-qJD(2) + t162) * mrSges(3,1)) * t167 + (-mrSges(3,2) * qJD(2) + t199 + (mrSges(3,2) - t198) * t162) * t170) * t247 + t266 * t72 + t265 * t71 + t220 * t198 + t171 + (qJ(3) * t104 + qJD(3) * t109 + (-t243 + (t110 * t166 - t111 * t169) * t162) * qJD(4) * mrSges(5,3)) * t163 + (t138 * t92 + t4 * t43 + t42 * t5 + (-t163 * t217 + t132) * t99 + t265 * t13 + t266 * t12) * m(6) + (t110 * t37 + t111 * t36 + t263 * t61 + t264 * t60 + t260) * m(5) + (-pkin(2) * t195 + (qJ(3) * t136 + t226) * t161 - (t143 * t161 * t170 + t137 * t167) * t247 + t260) * m(4); t261 * t72 + t262 * t71 + t184 * qJD(4) + (t140 * t57 - t183 * t56) * mrSges(6,3) + (-t198 * t162 - t164 * t184 + t199) * t162 + (t261 * t12 + t262 * t13 + t140 * t4 - t183 * t5 - t99 * t237) * m(6) + (t149 * t185 + t166 * t36 + t169 * t37 - t118) * m(5) + (-t161 * t240 - t118 + t195) * m(4); (t71 * t221 - t72 * t222 + m(6) * (-t12 * t222 + t13 * t221 + t165 * t4 + t168 * t5) + (t165 * t57 + t168 * t56) * mrSges(6,3)) * pkin(4) - m(6) * (t12 * t16 + t13 * t17) - t60 * t107 + t61 * t108 - t17 * t71 - t16 * t72 + t13 * t246 + (t166 * t78 / 0.2e1 - t149 * t186 / 0.2e1 + (-t143 * t188 + (-t181 / 0.2e1 + t182 / 0.2e1) * t162) * t163 + t177 + t185 * mrSges(5,3) + (t77 / 0.2e1 + (-m(6) * t99 - t59) * pkin(4)) * t169) * t237 + t172 - t267; -t12 * t71 + (t72 + t246) * t13 + t172;];
tauc = t1(:);
