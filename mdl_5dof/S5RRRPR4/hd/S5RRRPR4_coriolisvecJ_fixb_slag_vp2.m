% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:06
% DurationCPUTime: 2.55s
% Computational Cost: add. (2487->298), mult. (4302->408), div. (0->0), fcn. (2224->6), ass. (0->160)
t141 = qJD(1) + qJD(2);
t148 = sin(qJ(2));
t216 = pkin(1) * qJD(1);
t118 = pkin(7) * t141 + t148 * t216;
t147 = sin(qJ(3));
t252 = t147 * t118 + qJD(4);
t233 = qJD(3) - qJD(5);
t251 = -t233 / 0.2e1;
t250 = mrSges(5,1) + mrSges(4,1);
t207 = t141 * t147;
t249 = -pkin(8) * t207 + t252;
t150 = cos(qJ(3));
t151 = cos(qJ(2));
t215 = pkin(1) * qJD(2);
t187 = qJD(1) * t215;
t180 = t151 * t187;
t167 = t147 * t180;
t196 = qJD(3) * t150;
t52 = t118 * t196 + t167;
t214 = t147 * t52;
t123 = t150 * t180;
t197 = qJD(3) * t147;
t51 = -t118 * t197 + t123;
t248 = t150 * t51 + t214;
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t165 = t146 * t147 + t149 * t150;
t85 = t165 * t141;
t230 = t85 / 0.2e1;
t247 = t118 * (t147 ^ 2 + t150 ^ 2);
t142 = qJD(3) * qJD(4);
t45 = t142 + t51;
t246 = t150 * t45 + t214;
t206 = t141 * t150;
t188 = mrSges(5,2) * t206;
t117 = qJD(3) * mrSges(5,3) + t188;
t199 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t206 + t117;
t200 = (-mrSges(5,2) - mrSges(4,3)) * t207 + t250 * qJD(3);
t245 = t147 * t199 + t150 * t200;
t227 = pkin(3) + pkin(4);
t42 = -qJD(3) * t227 + t249;
t143 = qJD(3) * qJ(4);
t105 = t150 * t118;
t69 = -pkin(8) * t206 + t105;
t50 = t143 + t69;
t19 = -t146 * t50 + t149 * t42;
t182 = pkin(8) * t141 - t118;
t33 = t182 * t197 + t123 + t142;
t34 = -t182 * t196 + t167;
t2 = qJD(5) * t19 + t146 * t34 + t149 * t33;
t234 = -t146 * t150 + t147 * t149;
t86 = t234 * t141;
t228 = t86 / 0.2e1;
t75 = Ifges(6,4) * t86;
t25 = -Ifges(6,2) * t85 - Ifges(6,6) * t233 + t75;
t74 = Ifges(6,4) * t85;
t26 = Ifges(6,1) * t86 - Ifges(6,5) * t233 - t74;
t44 = t233 * t165;
t28 = t44 * t141;
t43 = t233 * t234;
t29 = t43 * t141;
t20 = t146 * t42 + t149 * t50;
t3 = -qJD(5) * t20 - t146 * t33 + t149 * t34;
t137 = t147 * qJ(4);
t183 = pkin(2) + t137;
t191 = t151 * t216;
t37 = t191 + (t150 * t227 + t183) * t141;
t244 = -(Ifges(6,1) * t85 + t75) * t86 / 0.2e1 - t3 * mrSges(6,1) + t2 * mrSges(6,2) - Ifges(6,5) * t28 - Ifges(6,6) * t29 - t25 * t228 + (-Ifges(6,5) * t85 - Ifges(6,6) * t86) * t251 + t37 * (t86 * mrSges(6,1) - t85 * mrSges(6,2)) + (Ifges(6,2) * t86 - t26 + t74) * t230;
t120 = -t146 * qJ(4) - t149 * t227;
t242 = qJD(5) * t120 - t146 * t69 + t249 * t149;
t121 = t149 * qJ(4) - t146 * t227;
t241 = -qJD(5) * t121 - t249 * t146 - t149 * t69;
t139 = t150 * pkin(3);
t236 = t137 + t139;
t80 = -qJD(3) * pkin(3) + t252;
t94 = t105 + t143;
t231 = t200 * t147 - m(5) * (t147 * t80 + t150 * t94) - t199 * t150;
t226 = pkin(7) - pkin(8);
t223 = t80 * mrSges(5,2);
t132 = pkin(1) * t148 + pkin(7);
t222 = -pkin(8) + t132;
t31 = mrSges(6,1) * t85 + mrSges(6,2) * t86;
t172 = -mrSges(5,1) * t150 - mrSges(5,3) * t147;
t98 = t172 * t141;
t221 = -t31 + t98;
t220 = Ifges(4,4) * t147;
t218 = Ifges(5,5) * t147;
t217 = Ifges(5,5) * t150;
t209 = qJ(4) * t150;
t201 = t94 * qJD(3);
t198 = qJD(2) * t148;
t136 = t147 * qJD(4);
t195 = -qJD(1) - t141;
t194 = -qJD(2) + t141;
t192 = pkin(2) + t236;
t190 = pkin(1) * t198;
t189 = t151 * t215;
t128 = t226 * t150;
t97 = pkin(3) * t197 - qJ(4) * t196 - t136;
t133 = -pkin(1) * t151 - pkin(2);
t109 = t222 * t150;
t185 = t197 / 0.2e1;
t181 = t148 * t187;
t179 = t141 * t185;
t175 = -t19 * t44 - t234 * t3;
t174 = -t19 * t85 + t20 * t86;
t173 = -mrSges(4,1) * t150 + mrSges(4,2) * t147;
t171 = pkin(3) * t147 - t209;
t47 = mrSges(6,2) * t233 - t85 * mrSges(6,3);
t48 = -mrSges(6,1) * t233 - t86 * mrSges(6,3);
t169 = -t146 * t48 + t149 * t47;
t168 = -t147 * t94 + t150 * t80;
t106 = t133 - t236;
t108 = t222 * t147;
t35 = t108 * t149 - t109 * t146;
t36 = t108 * t146 + t109 * t149;
t127 = t226 * t147;
t57 = t127 * t149 - t128 * t146;
t58 = t127 * t146 + t128 * t149;
t164 = -t147 * t227 + t209;
t163 = (t147 * Ifges(5,1) - t217) * t141;
t162 = (t150 * Ifges(4,2) + t220) * t141;
t70 = -pkin(4) * t197 - t97;
t161 = (mrSges(4,1) * t147 + mrSges(4,2) * t150) * qJD(3);
t160 = (mrSges(5,1) * t147 - mrSges(5,3) * t150) * qJD(3);
t154 = m(5) * t246 + m(4) * t248 + (m(5) * t168 - t245) * qJD(3);
t119 = -pkin(2) * t141 - t191;
t30 = -t181 + (qJD(3) * t164 + t136) * t141;
t38 = t181 + (qJD(3) * t171 - t136) * t141;
t49 = -t191 + (-t183 - t139) * t141;
t129 = Ifges(5,5) * t207;
t81 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t206 + t129;
t82 = Ifges(4,6) * qJD(3) + t162;
t83 = Ifges(5,4) * qJD(3) + t163;
t130 = Ifges(4,4) * t206;
t84 = Ifges(4,1) * t207 + Ifges(4,5) * qJD(3) + t130;
t153 = t30 * (mrSges(6,1) * t165 + mrSges(6,2) * t234) + (Ifges(6,5) * t44 + Ifges(6,6) * t43) * t251 + t43 * t25 / 0.2e1 + t44 * t26 / 0.2e1 + t37 * (-mrSges(6,1) * t43 + mrSges(6,2) * t44) + t173 * t181 + t49 * t160 + t119 * t161 + t38 * t172 + t81 * t185 + (-t150 * Ifges(5,3) + t218) * t179 + 0.2e1 * (t218 - t220 + (Ifges(4,1) + Ifges(5,1)) * t150) * t179 + (-t141 * (Ifges(5,3) * t147 + t217) + t223) * t196 + ((Ifges(5,4) + Ifges(4,5)) * t150 + (-Ifges(4,6) + Ifges(5,6)) * t147) * qJD(3) ^ 2 / 0.2e1 - (t162 + t82) * t197 / 0.2e1 + (-t29 * t165 - t43 * t230) * Ifges(6,2) + (t44 * t228 + t234 * t28) * Ifges(6,1) + (-t165 * t2 + t20 * t43) * mrSges(6,3) + t248 * mrSges(4,3) + (-t201 * t147 + t246) * mrSges(5,2) + (t29 * t234 - t28 * t165 - t85 * t44 / 0.2e1 + t43 * t228) * Ifges(6,4) + ((0.3e1 * Ifges(4,4) * t150 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t147) * t141 + t163 + t84 + t83) * t196 / 0.2e1;
t138 = t150 * pkin(4);
t134 = pkin(8) * t197;
t113 = qJD(3) * t128;
t112 = -pkin(7) * t197 + t134;
t107 = t138 + t192;
t100 = t171 * t141;
t99 = t173 * t141;
t93 = -t106 + t138;
t88 = t141 * t161;
t87 = t141 * t160;
t73 = t165 * t191;
t72 = t234 * t191;
t71 = t97 + t190;
t63 = t164 * t141;
t61 = qJD(3) * t109 + t147 * t189;
t60 = -t132 * t197 + t150 * t189 + t134;
t46 = t70 - t190;
t22 = -qJD(5) * t58 - t112 * t146 + t113 * t149;
t21 = qJD(5) * t57 + t112 * t149 + t113 * t146;
t10 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t9 = -qJD(5) * t36 - t146 * t60 + t149 * t61;
t8 = qJD(5) * t35 + t146 * t61 + t149 * t60;
t1 = [(-t28 * t35 + t29 * t36 + t175) * mrSges(6,3) + m(5) * (t106 * t38 + t49 * t71) + t133 * t88 + t106 * t87 + t93 * t10 + t71 * t98 + m(6) * (t19 * t9 + t2 * t36 + t20 * t8 + t3 * t35 + t30 * t93 + t37 * t46) + t153 + t46 * t31 + t8 * t47 + t9 * t48 + t154 * t132 + ((m(4) * (qJD(1) * t133 + t119) + t99 + t195 * mrSges(3,1)) * t148 + (m(4) * t247 + t195 * mrSges(3,2) - t231) * t151) * t215; (-t28 * t57 + t29 * t58 + t175) * mrSges(6,3) - m(6) * (t19 * t72 + t20 * t73) - t192 * t87 + t107 * t10 - pkin(2) * t88 + t97 * t98 + m(6) * (t107 * t30 + t19 * t22 + t2 * t58 + t20 * t21 + t3 * t57 + t37 * t70) + m(5) * (-t192 * t38 + t49 * t97) + t153 + t70 * t31 + t154 * pkin(7) + (-t72 + t22) * t48 + (-t73 + t21) * t47 + ((t194 * mrSges(3,2) + t231) * t151 + (-m(5) * t49 + m(6) * t37 + mrSges(3,1) * t194 - t221 - t99) * t148 + (-pkin(2) * t198 - t119 * t148 - t151 * t247) * m(4)) * t216; t241 * t48 + t242 * t47 + ((-t119 * mrSges(4,2) + t49 * mrSges(5,3) - t130 / 0.2e1 - t83 / 0.2e1 - t84 / 0.2e1 - t223 + Ifges(5,5) * t206 / 0.2e1) * t150 + (-t119 * mrSges(4,1) - t49 * mrSges(5,1) - t129 / 0.2e1 - t81 / 0.2e1 + t82 / 0.2e1 + t94 * mrSges(5,2) + Ifges(4,4) * t207 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1) * t206) * t147 + ((-pkin(3) * mrSges(5,2) + Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t150 + (-qJ(4) * mrSges(5,2) + Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t147) * qJD(3)) * t141 + (-t120 * t28 + t121 * t29 - t174) * mrSges(6,3) + t245 * t118 + qJD(4) * t117 - t100 * t98 - t63 * t31 + t45 * mrSges(5,3) - t51 * mrSges(4,2) - t250 * t52 + (t120 * t3 + t121 * t2 + t19 * t241 + t20 * t242 - t37 * t63) * m(6) + (-pkin(3) * t52 + qJ(4) * t45 + qJD(4) * t94 - t100 * t49 - t118 * t168) * m(5) + t244; t221 * t207 + t169 * qJD(5) + (t146 * t29 - t149 * t28) * mrSges(6,3) + (-t117 - t169 + t188) * qJD(3) + (t2 * t146 + t3 * t149 - t37 * t207 - t233 * (-t146 * t19 + t149 * t20)) * m(6) + (t207 * t49 - t201 + t52) * m(5); t174 * mrSges(6,3) - t19 * t47 + t20 * t48 - t244;];
tauc = t1(:);
