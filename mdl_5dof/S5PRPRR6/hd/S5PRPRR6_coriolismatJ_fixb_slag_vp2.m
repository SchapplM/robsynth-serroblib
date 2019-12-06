% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:33
% EndTime: 2019-12-05 15:56:38
% DurationCPUTime: 2.05s
% Computational Cost: add. (4882->254), mult. (11618->376), div. (0->0), fcn. (12776->10), ass. (0->152)
t258 = m(6) / 0.2e1;
t163 = sin(pkin(10));
t165 = cos(pkin(10));
t212 = t163 ^ 2 + t165 ^ 2;
t269 = t212 * mrSges(4,3);
t248 = sin(qJ(4));
t250 = cos(qJ(4));
t143 = -t163 * t250 - t165 * t248;
t167 = cos(qJ(5));
t233 = t167 * mrSges(6,2);
t166 = sin(qJ(5));
t237 = t166 * mrSges(6,1);
t146 = t233 + t237;
t102 = t146 * t143;
t268 = t143 * mrSges(5,3) + t102;
t164 = sin(pkin(5));
t168 = cos(qJ(2));
t218 = t164 * t168;
t121 = t143 * t218;
t141 = t163 * t248 - t165 * t250;
t122 = t141 * t218;
t249 = sin(qJ(2));
t203 = t249 * t164;
t99 = -t122 * t167 + t166 * t203;
t224 = t99 * t167;
t98 = t122 * t166 + t167 * t203;
t225 = t98 * t166;
t267 = (t224 / 0.2e1 - t225 / 0.2e1) * mrSges(6,3) + (pkin(4) * t121 + (t224 - t225) * pkin(8)) * t258 + t122 * mrSges(5,2) / 0.2e1;
t266 = -m(6) * (t166 * t99 + t167 * t98) / 0.2e1 - (m(4) + m(5)) * t203 / 0.2e1;
t161 = t166 ^ 2;
t162 = t167 ^ 2;
t264 = mrSges(6,3) * (t162 / 0.2e1 + t161 / 0.2e1);
t156 = Ifges(6,5) * t167;
t241 = Ifges(6,6) * t166;
t186 = -t156 / 0.2e1 + t241 / 0.2e1;
t263 = Ifges(5,4) + t186;
t157 = Ifges(6,4) * t167;
t148 = Ifges(6,1) * t166 + t157;
t223 = cos(pkin(5));
t132 = t163 * t223 + t165 * t203;
t178 = t163 * t203 - t165 * t223;
t262 = t165 * t132 + t163 * t178;
t246 = t143 * pkin(4);
t247 = t141 * pkin(8);
t110 = -t246 + t247;
t245 = pkin(7) + qJ(3);
t144 = t245 * t165;
t201 = t245 * t163;
t111 = t144 * t248 + t201 * t250;
t58 = t110 * t167 + t111 * t166;
t59 = t110 * t166 - t111 * t167;
t261 = -t58 * t166 + t59 * t167;
t195 = Ifges(6,2) * t166 - t157;
t145 = -mrSges(6,1) * t167 + mrSges(6,2) * t166;
t259 = -m(6) * pkin(4) - mrSges(5,1) + t145;
t86 = t132 * t250 - t178 * t248;
t70 = -t166 * t86 - t167 * t218;
t257 = t70 / 0.2e1;
t256 = -t145 / 0.2e1;
t255 = t146 / 0.2e1;
t254 = -t166 / 0.2e1;
t253 = t166 / 0.2e1;
t252 = -t167 / 0.2e1;
t251 = t167 / 0.2e1;
t243 = Ifges(6,4) * t166;
t240 = t141 * mrSges(5,3);
t85 = t132 * t248 + t178 * t250;
t238 = t143 * t85;
t236 = t166 * mrSges(6,3);
t235 = t166 * t70;
t75 = Ifges(6,6) * t141 + t143 * t195;
t234 = t166 * t75;
t232 = t167 * mrSges(6,3);
t155 = -pkin(3) * t165 - pkin(2);
t104 = pkin(4) * t141 + pkin(8) * t143 + t155;
t112 = t144 * t250 - t201 * t248;
t52 = t104 * t166 + t112 * t167;
t231 = t167 * t52;
t71 = -t166 * t218 + t167 * t86;
t230 = t167 * t71;
t196 = Ifges(6,1) * t167 - t243;
t77 = Ifges(6,5) * t141 - t143 * t196;
t229 = t167 * t77;
t226 = t85 * t121;
t222 = t111 * t121;
t221 = t111 * t143;
t204 = t249 * t164 ^ 2;
t12 = m(6) * (t70 * t98 + t71 * t99 - t226) + m(4) * (t164 * t262 - t204) * t168 + (-t122 * t86 - t168 * t204 - t226) * m(5);
t220 = t12 * qJD(1);
t190 = -t230 + t235;
t15 = m(6) * (t190 + t86) * t85;
t219 = t15 * qJD(1);
t108 = mrSges(6,1) * t141 + t143 * t232;
t216 = t166 * t108;
t208 = t143 * t236;
t106 = -mrSges(6,2) * t141 + t208;
t215 = t167 * t106;
t134 = t141 * mrSges(5,2);
t211 = -t161 - t162;
t174 = -t141 * t264 + (t211 * t247 + t246) * t258;
t105 = mrSges(6,2) * t143 + t141 * t236;
t107 = -mrSges(6,1) * t143 + t141 * t232;
t176 = (t166 * t59 + t167 * t58) * t258 + t105 * t253 + t107 * t251;
t18 = t134 + (t256 + mrSges(5,1)) * t143 + t174 - t176;
t214 = t18 * qJD(2);
t185 = t233 / 0.2e1 + t237 / 0.2e1;
t180 = t185 * t141;
t182 = -t215 / 0.2e1 + t216 / 0.2e1;
t22 = t180 + t182;
t213 = t22 * qJD(2);
t207 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t199 = t156 - t241;
t198 = t212 * qJ(3);
t147 = Ifges(6,2) * t167 + t243;
t194 = Ifges(6,5) * t166 + Ifges(6,6) * t167;
t101 = t146 * t141;
t109 = -t143 * mrSges(5,1) - t134;
t51 = t104 * t167 - t112 * t166;
t191 = t166 * t51 - t231;
t169 = ((t112 + t191) * t258 - t101 / 0.2e1 + t182) * t85 + (t111 * t86 + t58 * t70 + t59 * t71) * t258 + t107 * t257 + t71 * t105 / 0.2e1 - t86 * t102 / 0.2e1 - t109 * t218 / 0.2e1;
t2 = -(t256 + mrSges(5,1) / 0.2e1) * t121 + t169 - t267;
t74 = -Ifges(6,6) * t143 + t141 * t195;
t76 = -Ifges(6,5) * t143 - t141 * t196;
t3 = -t112 * t102 - t111 * t101 + t59 * t106 + t52 * t105 + t58 * t108 + t51 * t107 + m(6) * (t111 * t112 + t51 * t58 + t52 * t59) + t155 * t109 + (-t229 / 0.2e1 + t234 / 0.2e1 + t263 * t141) * t141 + (t76 * t252 + t74 * t253 - t263 * t143 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2)) * t141) * t143;
t193 = t2 * qJD(1) + t3 * qJD(2);
t100 = t143 * t145;
t4 = -t111 * t100 + t52 * t108 + (t77 * t254 + t75 * t252 - mrSges(6,3) * t231 - t141 * t194 / 0.2e1 + (t147 * t254 + t148 * t251) * t143) * t143 + (-t106 + t208) * t51;
t173 = (t230 / 0.2e1 - t235 / 0.2e1) * t143 * mrSges(6,3) + t106 * t257 - t71 * t108 / 0.2e1 + t85 * t100 / 0.2e1;
t187 = t98 * mrSges(6,1) / 0.2e1 - t99 * mrSges(6,2) / 0.2e1;
t7 = t173 - t187;
t192 = t7 * qJD(1) - t4 * qJD(2);
t10 = t268 * t143 + t269 + (-t215 + t216 + t240) * t141 + m(6) * (t141 * t191 - t221) + m(5) * (-t112 * t141 - t221) + m(4) * t198;
t170 = (t141 * t190 - t238) * t258 + m(5) * (-t141 * t86 - t238) / 0.2e1 + m(4) * t262 / 0.2e1;
t14 = t170 + t266;
t189 = -qJD(1) * t14 - qJD(2) * t10;
t188 = -t58 * mrSges(6,1) / 0.2e1 + t59 * mrSges(6,2) / 0.2e1;
t183 = t147 * t253 + t148 * t252;
t181 = t143 * t194;
t16 = (-t146 / 0.2e1 + t185) * t85;
t57 = pkin(4) * t146 - t195 * t252 + t196 * t254 + t183;
t171 = pkin(8) * t264 + (t148 / 0.4e1 + t157 / 0.4e1 + t207 * t166) * t166 + (0.3e1 / 0.4e1 * t243 + t147 / 0.4e1 - t207 * t167) * t167;
t172 = (t106 * t254 + t108 * t252) * pkin(8) - pkin(4) * t100 / 0.2e1 + t111 * t255 - t234 / 0.4e1 + t229 / 0.4e1;
t6 = (-0.3e1 / 0.4e1 * t241 + 0.3e1 / 0.4e1 * t156) * t141 + (Ifges(6,3) / 0.2e1 + t171) * t143 + t172 + t188;
t179 = t16 * qJD(1) - t6 * qJD(2) + t57 * qJD(4);
t23 = t180 - t182;
t19 = t143 * t256 + t174 + t176;
t17 = (t185 + t255) * t85;
t13 = t170 - t266;
t8 = t173 + t187;
t5 = t172 - t188 + (-Ifges(6,3) / 0.2e1 + t171) * t143 + (t199 / 0.4e1 + t186) * t141;
t1 = t169 - (t145 / 0.2e1 - mrSges(5,1) / 0.2e1) * t121 + t267;
t9 = [qJD(2) * t12 + qJD(4) * t15, t13 * qJD(3) + t1 * qJD(4) + t8 * qJD(5) + t220 + (m(6) * (t51 * t98 + t52 * t99 - t222) + t99 * t106 + t98 * t108 + t122 * t240 + m(5) * (-t112 * t122 - t222) + m(4) * (-pkin(2) * t249 + t168 * t198) * t164 + (m(5) * t155 - mrSges(4,1) * t165 + mrSges(5,1) * t141 + mrSges(4,2) * t163 - mrSges(5,2) * t143 - mrSges(3,1)) * t203 + t268 * t121 + (-mrSges(3,2) + t269) * t218) * qJD(2), qJD(2) * t13, t219 + t1 * qJD(2) + (t259 * t86 + (mrSges(5,2) + (m(6) * pkin(8) + mrSges(6,3)) * t211) * t85) * qJD(4) + t17 * qJD(5), t8 * qJD(2) + t17 * qJD(4) + (-mrSges(6,1) * t71 - mrSges(6,2) * t70) * qJD(5); qJD(3) * t14 + qJD(4) * t2 + qJD(5) * t7 - t220, qJD(3) * t10 + qJD(4) * t3 - qJD(5) * t4, qJD(4) * t19 + qJD(5) * t23 - t189, t19 * qJD(3) + t5 * qJD(5) + t193 + (pkin(4) * t101 - t181 / 0.2e1 + t76 * t253 + t74 * t251 + Ifges(5,6) * t143 + t111 * mrSges(5,2) + (-Ifges(5,5) + t183) * t141 + t259 * t112 + (m(6) * t261 + t167 * t105 - t166 * t107) * pkin(8) + t261 * mrSges(6,3)) * qJD(4), t23 * qJD(3) + t5 * qJD(4) + (-t52 * mrSges(6,1) - t51 * mrSges(6,2) + t181) * qJD(5) + t192; -qJD(2) * t14, -qJD(4) * t18 - qJD(5) * t22 + t189, 0, -t214, -qJD(5) * t146 - t213; -qJD(2) * t2 - qJD(5) * t16 - t219, qJD(3) * t18 + qJD(5) * t6 - t193, t214, -t57 * qJD(5), (pkin(8) * t145 + t199) * qJD(5) - t179; -t7 * qJD(2) + t16 * qJD(4), qJD(3) * t22 - qJD(4) * t6 - t192, t213, t179, 0;];
Cq = t9;
