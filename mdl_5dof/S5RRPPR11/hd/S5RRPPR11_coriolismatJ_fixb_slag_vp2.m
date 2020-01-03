% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:29
% EndTime: 2019-12-31 19:46:34
% DurationCPUTime: 2.36s
% Computational Cost: add. (5374->318), mult. (10407->436), div. (0->0), fcn. (10003->6), ass. (0->167)
t264 = pkin(3) + pkin(6);
t280 = -Ifges(3,4) - Ifges(4,6);
t182 = cos(qJ(2));
t177 = sin(pkin(8));
t178 = cos(pkin(8));
t180 = sin(qJ(5));
t249 = cos(qJ(5));
t277 = -t177 * t180 + t178 * t249;
t119 = t277 * t182;
t279 = t119 * mrSges(6,3);
t196 = t177 * t249 + t178 * t180;
t137 = Ifges(6,4) * t196;
t96 = Ifges(6,1) * t277 - t137;
t278 = -Ifges(6,2) * t277 - t137 + t96;
t181 = sin(qJ(2));
t208 = pkin(2) * t181 - qJ(3) * t182;
t139 = qJ(4) * t181 + t208;
t159 = t264 * t182;
t148 = t178 * t159;
t66 = pkin(4) * t182 + t148 + (-pkin(7) * t181 - t139) * t177;
t218 = t178 * t181;
t90 = t139 * t178 + t159 * t177;
t78 = pkin(7) * t218 + t90;
t31 = -t180 * t78 + t249 * t66;
t32 = t180 * t66 + t249 * t78;
t276 = t196 * t32 + t277 * t31;
t270 = t178 ^ 2;
t275 = -t177 ^ 2 - t270;
t274 = -t196 ^ 2 - t277 ^ 2;
t269 = -m(5) / 0.2e1;
t268 = m(5) / 0.2e1;
t267 = -m(6) / 0.2e1;
t266 = m(6) / 0.2e1;
t265 = -mrSges(6,1) / 0.2e1;
t118 = t277 * t181;
t263 = t118 / 0.2e1;
t262 = -t119 / 0.2e1;
t120 = t196 * t181;
t260 = t120 / 0.2e1;
t121 = t196 * t182;
t259 = -t121 / 0.2e1;
t258 = t277 / 0.2e1;
t257 = -t277 / 0.2e1;
t256 = t196 / 0.2e1;
t255 = -t196 / 0.2e1;
t253 = -t177 / 0.2e1;
t252 = t177 / 0.2e1;
t251 = -t178 / 0.2e1;
t250 = t178 / 0.2e1;
t179 = -pkin(2) - qJ(4);
t248 = -pkin(7) + t179;
t247 = Ifges(5,4) * t177;
t246 = Ifges(5,4) * t178;
t245 = Ifges(6,4) * t277;
t100 = -mrSges(6,2) * t181 - t279;
t101 = mrSges(6,1) * t182 - mrSges(6,3) * t120;
t102 = mrSges(6,1) * t181 + t121 * mrSges(6,3);
t116 = t182 * Ifges(5,6) + (Ifges(5,2) * t178 + t247) * t181;
t237 = t177 * Ifges(5,1);
t117 = t182 * Ifges(5,5) + (t237 + t246) * t181;
t129 = (-pkin(4) * t178 - t264) * t181;
t217 = t178 * t182;
t130 = pkin(4) * t217 + t159;
t235 = t178 * mrSges(5,1);
t238 = t177 * mrSges(5,2);
t206 = t235 - t238;
t131 = t206 * t181;
t219 = t177 * t181;
t231 = t182 * mrSges(5,1);
t149 = -mrSges(5,3) * t219 + t231;
t233 = t181 * mrSges(5,1);
t150 = mrSges(5,3) * t177 * t182 + t233;
t230 = t182 * mrSges(5,2);
t151 = mrSges(5,3) * t218 - t230;
t232 = t181 * mrSges(5,2);
t152 = -mrSges(5,3) * t217 - t232;
t226 = qJ(3) * t181;
t204 = -pkin(2) * t182 - t226;
t155 = -pkin(1) + t204;
t158 = t264 * t181;
t198 = Ifges(6,5) * t260 + Ifges(6,6) * t263;
t156 = t182 * mrSges(4,2) - t181 * mrSges(4,3);
t210 = m(4) * t155 + t156;
t234 = t178 * Ifges(5,6);
t236 = t177 * Ifges(5,5);
t240 = t121 * mrSges(6,2);
t242 = t119 * mrSges(6,1);
t138 = t179 * t182 - pkin(1) - t226;
t147 = t178 * t158;
t65 = pkin(4) * t181 + t147 + (pkin(7) * t182 - t138) * t177;
t88 = t138 * t178 + t158 * t177;
t77 = -pkin(7) * t217 + t88;
t29 = -t180 * t77 + t249 * t65;
t30 = t180 * t65 + t249 * t77;
t61 = Ifges(6,4) * t120 + Ifges(6,2) * t118 + t182 * Ifges(6,6);
t239 = t121 * Ifges(6,4);
t62 = -t119 * Ifges(6,2) + t181 * Ifges(6,6) - t239;
t63 = Ifges(6,1) * t120 + Ifges(6,4) * t118 + t182 * Ifges(6,5);
t115 = Ifges(6,4) * t119;
t64 = -t121 * Ifges(6,1) + Ifges(6,5) * t181 - t115;
t241 = t120 * mrSges(6,2);
t243 = t118 * mrSges(6,1);
t68 = t241 - t243;
t87 = -t138 * t177 + t147;
t89 = -t139 * t177 + t148;
t99 = -mrSges(6,2) * t182 + t118 * mrSges(6,3);
t1 = t90 * t152 - t159 * t131 + t87 * t149 + t89 * t150 + t88 * t151 + t130 * t68 + t64 * t260 + t63 * t259 + t129 * (-t240 + t242) + t62 * t263 + t61 * t262 + t30 * t99 + t32 * t100 + t29 * t101 + t31 * t102 + t210 * t208 + m(6) * (t129 * t130 + t29 * t31 + t30 * t32) + m(5) * (-t158 * t159 + t87 * t89 + t88 * t90) + (-pkin(1) * mrSges(3,2) - t155 * mrSges(4,3) + Ifges(6,5) * t259 + Ifges(6,6) * t262 + t116 * t251 + t117 * t253 - t158 * t206 + (-t236 / 0.2e1 - t234 / 0.2e1 - t280) * t182) * t182 + (-pkin(1) * mrSges(3,1) - t155 * mrSges(4,2) + (t234 + t236 + t280) * t181 + (-Ifges(3,2) - Ifges(4,3) + Ifges(3,1) + Ifges(5,3) + Ifges(6,3) + Ifges(4,2) - t270 * Ifges(5,2) / 0.2e1 + (-t246 - t237 / 0.2e1) * t177) * t182 + t198) * t181;
t244 = t1 * qJD(1);
t215 = -Ifges(6,5) * t119 + Ifges(6,6) * t121;
t67 = -mrSges(6,1) * t121 - mrSges(6,2) * t119;
t69 = Ifges(6,2) * t121 - t115;
t70 = -Ifges(6,1) * t119 + t239;
t4 = t130 * t67 + t181 * t215 / 0.2e1 + t29 * t100 - t30 * t102 + (t30 * mrSges(6,3) - t70 / 0.2e1 + t62 / 0.2e1) * t121 - (-t29 * mrSges(6,3) + t64 / 0.2e1 + t69 / 0.2e1) * t119;
t227 = t4 * qJD(1);
t190 = (-t119 * t257 + t121 * t256) * mrSges(6,3) + t100 * t258 + t102 * t255;
t197 = -t118 * mrSges(6,2) / 0.2e1 + t120 * t265;
t10 = t190 + t197;
t225 = t10 * qJD(1);
t224 = t118 * t196;
t223 = t119 * t196;
t191 = m(5) * (t177 * t87 - t178 * t88) - t178 * t152 + t177 * t150;
t12 = m(6) * (-t118 * t30 + t120 * t29) - t118 * t100 + t120 * t102 + (t191 - t210) * t181;
t222 = t12 * qJD(1);
t221 = t120 * t277;
t220 = t121 * t277;
t214 = -Ifges(6,5) * t196 - Ifges(6,6) * t277;
t213 = t269 + t267;
t212 = (t220 - t223) * t266;
t94 = mrSges(6,1) * t196 + mrSges(6,2) * t277;
t207 = m(5) * t275;
t205 = -Ifges(6,1) * t196 - t245;
t203 = t177 * t90 + t178 * t89;
t166 = pkin(4) * t177 + qJ(3);
t93 = mrSges(6,1) * t277 - mrSges(6,2) * t196;
t95 = -Ifges(6,2) * t196 + t245;
t13 = -t166 * t93 + t205 * t257 + t256 * t278 + t258 * t95;
t153 = t248 * t177;
t154 = t248 * t178;
t91 = -t153 * t180 + t154 * t249;
t92 = t153 * t249 + t154 * t180;
t185 = t130 * t93 / 0.2e1 + t166 * t67 / 0.2e1 + t181 * t214 / 0.4e1 - t92 * t102 / 0.2e1 + (t100 / 0.2e1 + t279 / 0.2e1) * t91 - t278 * t119 / 0.4e1 - (t64 + t69) * t196 / 0.4e1 - (t62 / 0.4e1 - t70 / 0.4e1) * t277 + (t92 * mrSges(6,3) / 0.2e1 + t95 / 0.4e1 - t205 / 0.4e1) * t121;
t189 = -Ifges(6,3) * t182 / 0.2e1 + t31 * t265 + t32 * mrSges(6,2) / 0.2e1 - t198;
t3 = t185 + t189;
t202 = qJD(1) * t3 - qJD(2) * t13;
t21 = -m(6) * (-t196 * t92 - t277 * t91) - t179 * t207 + t274 * mrSges(6,3) + t275 * mrSges(5,3);
t186 = (-t220 / 0.2e1 + t223 / 0.2e1) * mrSges(6,3) + (-t177 * t88 - t178 * t87) * t268 + (-t119 * t92 + t121 * t91 - t196 * t30 - t277 * t29) * t266 + t102 * t257 + t100 * t255;
t192 = t158 * t268 + t129 * t267 + t243 / 0.2e1 - t241 / 0.2e1;
t7 = (-t150 / 0.2e1 + t233 / 0.2e1) * t178 + (-t152 / 0.2e1 - t232 / 0.2e1) * t177 + t186 + t192;
t201 = qJD(1) * t7 - qJD(2) * t21;
t44 = m(6) * t166 + mrSges(4,3) + t178 * mrSges(5,2) + t177 * mrSges(5,1) + (m(5) + m(4)) * qJ(3) + t94;
t187 = (-t221 / 0.2e1 + t224 / 0.2e1) * mrSges(6,3) + t159 * t268 + (-t118 * t92 + t120 * t91 + t130) * t266 + t242 / 0.2e1 - t240 / 0.2e1;
t188 = t101 * t257 + t203 * t269 + t255 * t99 + t267 * t276;
t8 = (t231 / 0.2e1 - t149 / 0.2e1) * t178 + (-t230 / 0.2e1 - t151 / 0.2e1) * t177 + t187 + t188;
t200 = qJD(1) * t8 + qJD(2) * t44;
t199 = qJD(1) * t67 + qJD(2) * t93;
t14 = m(6) * (-t119 * t30 + t121 * t29) - t119 * t100 + t121 * t102 + t191 * t182;
t195 = qJD(1) * t14 + qJD(3) * t212;
t193 = t274 * t266 + t207 / 0.2e1;
t35 = t193 + t213;
t194 = qJD(1) * t212 + qJD(2) * t35;
t37 = qJD(4) * t212;
t34 = t193 - t213;
t9 = t190 - t197;
t6 = t152 * t253 + t150 * t251 + mrSges(5,2) * t219 / 0.2e1 - mrSges(5,1) * t218 / 0.2e1 + t186 - t192;
t5 = t151 * t252 + t149 * t250 + (mrSges(4,1) - t238 / 0.2e1 + t235 / 0.2e1 + m(4) * pkin(6)) * t182 + t187 - t188;
t2 = t185 - t189;
t11 = [qJD(2) * t1 + qJD(3) * t12 + qJD(4) * t14 + qJD(5) * t4, t5 * qJD(3) + t6 * qJD(4) + t2 * qJD(5) + t244 + (-qJ(3) * t131 + t91 * t101 + t129 * t94 + t166 * t68 + t61 * t255 + t63 * t258 + t96 * t260 + t95 * t263 + t92 * t99 + (-t158 * mrSges(5,2) + t117 / 0.2e1 - t89 * mrSges(5,3) + t179 * t149) * t178 + (-t158 * mrSges(5,1) - t116 / 0.2e1 - t90 * mrSges(5,3) + t179 * t151) * t177 + 0.2e1 * (t129 * t166 + t31 * t91 + t32 * t92) * t266 + 0.2e1 * (-qJ(3) * t158 + t179 * t203) * t268 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t250 + Ifges(6,5) * t258 + Ifges(5,6) * t253 + Ifges(6,6) * t255 - Ifges(4,4) + Ifges(3,5)) * t182 + (Ifges(4,5) - Ifges(3,6) - qJ(3) * mrSges(4,1) + (Ifges(5,1) * t178 - t247) * t252 + (-Ifges(5,2) * t177 + t246) * t250) * t181 + (m(4) * t204 - t182 * mrSges(3,1) + t181 * mrSges(3,2) + t156) * pkin(6) - t276 * mrSges(6,3)) * qJD(2), t222 + t5 * qJD(2) + m(6) * (t221 - t224) * qJD(3) + t37 + t9 * qJD(5), qJD(2) * t6 + t195, t227 + t2 * qJD(2) + t9 * qJD(3) + (-mrSges(6,1) * t30 - mrSges(6,2) * t29 + t215) * qJD(5); qJD(3) * t8 + qJD(4) * t7 + qJD(5) * t3 - t244, qJD(3) * t44 - qJD(4) * t21 - qJD(5) * t13, qJD(4) * t34 + t200, qJD(3) * t34 + t201, (-mrSges(6,1) * t92 - mrSges(6,2) * t91 + t214) * qJD(5) + t202; -qJD(2) * t8 + qJD(5) * t10 - t222 + t37, qJD(4) * t35 - t200, 0, t194, -qJD(5) * t94 + t225; -qJD(2) * t7 + qJD(5) * t67 - t195, -qJD(3) * t35 + qJD(5) * t93 - t201, -t194, 0, t199; -qJD(2) * t3 - qJD(3) * t10 - qJD(4) * t67 - t227, -qJD(4) * t93 - t202, -t225, -t199, 0;];
Cq = t11;
