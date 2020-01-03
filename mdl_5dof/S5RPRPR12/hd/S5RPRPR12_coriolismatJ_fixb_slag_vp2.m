% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:30
% EndTime: 2019-12-31 18:29:36
% DurationCPUTime: 2.90s
% Computational Cost: add. (8695->303), mult. (17827->443), div. (0->0), fcn. (20106->8), ass. (0->162)
t179 = cos(pkin(8));
t245 = pkin(6) + qJ(2);
t169 = t245 * t179;
t181 = sin(qJ(3));
t177 = sin(pkin(8));
t201 = t245 * t177;
t247 = cos(qJ(3));
t144 = t169 * t247 - t181 * t201;
t263 = -m(5) / 0.2e1;
t162 = t177 * t181 - t179 * t247;
t176 = sin(pkin(9));
t178 = cos(pkin(9));
t180 = sin(qJ(5));
t182 = cos(qJ(5));
t197 = t176 * t180 - t178 * t182;
t113 = t197 * t162;
t234 = t113 * mrSges(6,2);
t164 = t176 * t182 + t178 * t180;
t110 = t164 * t162;
t236 = t110 * mrSges(6,1);
t266 = -t236 / 0.2e1 + t234 / 0.2e1;
t220 = t162 * t176;
t93 = -pkin(4) * t220 + t144;
t269 = t144 * t263 - m(6) * t93 / 0.2e1 - t266;
t165 = t177 * t247 + t179 * t181;
t109 = t197 * t165;
t256 = -t109 / 0.2e1;
t252 = -t164 / 0.2e1;
t268 = t109 * t252;
t174 = t178 ^ 2;
t267 = -t174 / 0.2e1;
t191 = t164 * t165;
t258 = -t191 / 0.2e1;
t223 = qJ(4) * t162;
t246 = pkin(3) * t165;
t135 = t223 + t246;
t142 = t169 * t181 + t201 * t247;
t72 = t135 * t178 + t142 * t176;
t73 = t135 * t176 - t142 * t178;
t264 = -t176 * t72 + t178 * t73;
t262 = m(5) / 0.2e1;
t261 = m(6) / 0.2e1;
t259 = t110 / 0.2e1;
t257 = t113 / 0.2e1;
t158 = Ifges(6,4) * t197;
t140 = Ifges(6,1) * t164 - t158;
t255 = t140 / 0.2e1;
t244 = pkin(7) + qJ(4);
t166 = t244 * t176;
t168 = t244 * t178;
t143 = -t166 * t180 + t168 * t182;
t254 = -t143 / 0.2e1;
t253 = -t197 / 0.2e1;
t251 = t164 / 0.2e1;
t250 = t165 / 0.2e1;
t249 = t176 / 0.2e1;
t248 = t178 / 0.2e1;
t242 = Ifges(5,4) * t176;
t241 = Ifges(5,4) * t178;
t240 = Ifges(6,4) * t109;
t239 = Ifges(6,4) * t164;
t154 = t162 * mrSges(4,2);
t141 = -t166 * t182 - t168 * t180;
t171 = -pkin(4) * t178 - pkin(3);
t172 = t176 ^ 2;
t208 = t172 + t174;
t185 = (t110 * t252 + t113 * t253) * mrSges(6,3) + (t267 - t172 / 0.2e1) * t162 * mrSges(5,3) + (-t208 * t223 - t246) * t262 + (t141 * t110 + t143 * t113 + t171 * t165) * t261;
t127 = -mrSges(5,2) * t165 + mrSges(5,3) * t220;
t219 = t162 * t178;
t129 = mrSges(5,1) * t165 + mrSges(5,3) * t219;
t44 = pkin(4) * t165 + pkin(7) * t219 + t72;
t53 = pkin(7) * t220 + t73;
t26 = -t180 * t53 + t182 * t44;
t27 = t180 * t44 + t182 * t53;
t75 = -mrSges(6,2) * t165 + mrSges(6,3) * t110;
t77 = mrSges(6,1) * t165 - t113 * mrSges(6,3);
t186 = (t176 * t73 + t178 * t72) * t262 + (t164 * t27 - t197 * t26) * t261 + t77 * t253 + t75 * t251 + t127 * t249 + t129 * t248;
t136 = mrSges(6,1) * t197 + mrSges(6,2) * t164;
t167 = -mrSges(5,1) * t178 + mrSges(5,2) * t176;
t203 = t136 / 0.2e1 + t167 / 0.2e1;
t8 = t154 + (-mrSges(4,1) + t203) * t165 + t185 - t186;
t238 = qJD(1) * t8;
t228 = t178 * mrSges(5,2);
t232 = t176 * mrSges(5,1);
t200 = t228 + t232;
t122 = t200 * t162;
t123 = t200 * t165;
t218 = t165 * t176;
t128 = -t162 * mrSges(5,2) - mrSges(5,3) * t218;
t217 = t165 * t178;
t130 = t162 * mrSges(5,1) - mrSges(5,3) * t217;
t194 = Ifges(6,5) * t257 + Ifges(6,6) * t259;
t204 = -pkin(2) * t179 - pkin(1);
t124 = pkin(3) * t162 - qJ(4) * t165 + t204;
t68 = t124 * t178 - t144 * t176;
t41 = pkin(4) * t162 - pkin(7) * t217 + t68;
t69 = t124 * t176 + t144 * t178;
t50 = -pkin(7) * t218 + t69;
t22 = -t180 * t50 + t182 * t41;
t227 = t178 * Ifges(5,5);
t23 = t180 * t41 + t182 * t50;
t230 = t176 * Ifges(5,6);
t231 = t176 * Ifges(5,2);
t45 = Ifges(6,4) * t113 + Ifges(6,2) * t110 + Ifges(6,6) * t165;
t46 = -Ifges(6,2) * t191 + Ifges(6,6) * t162 - t240;
t47 = Ifges(6,1) * t113 + Ifges(6,4) * t110 + Ifges(6,5) * t165;
t108 = Ifges(6,4) * t191;
t48 = -Ifges(6,1) * t109 + Ifges(6,5) * t162 - t108;
t55 = t234 - t236;
t56 = mrSges(6,1) * t191 - mrSges(6,2) * t109;
t76 = -mrSges(6,2) * t162 - mrSges(6,3) * t191;
t78 = mrSges(6,1) * t162 + t109 * mrSges(6,3);
t87 = Ifges(5,6) * t165 + (t231 - t241) * t162;
t88 = Ifges(5,5) * t165 + (-Ifges(5,1) * t178 + t242) * t162;
t92 = pkin(4) * t218 + t142;
t1 = t23 * t75 + t27 * t76 + t22 * t77 + t26 * t78 + t92 * t55 + t93 * t56 + t46 * t259 + t45 * t258 + t48 * t257 + t47 * t256 + t69 * t127 + t73 * t128 + t68 * t129 + t72 * t130 - t142 * t122 + t144 * t123 - t204 * t154 + m(5) * (t142 * t144 + t68 * t72 + t69 * t73) + m(6) * (t22 * t26 + t23 * t27 + t92 * t93) + (Ifges(6,5) * t256 + Ifges(6,6) * t258 + t88 * t248 - t176 * t87 / 0.2e1 + t204 * mrSges(4,1) + (-Ifges(4,4) + t227 / 0.2e1 - t230 / 0.2e1) * t165) * t165 + ((Ifges(5,1) * t267 + Ifges(6,3) - Ifges(4,1) + Ifges(4,2) + Ifges(5,3) + (t241 - t231 / 0.2e1) * t176) * t165 + t194 + (Ifges(4,4) - t227 + t230) * t162) * t162;
t237 = t1 * qJD(1);
t210 = -Ifges(6,5) * t191 + Ifges(6,6) * t109;
t54 = -t109 * mrSges(6,1) - mrSges(6,2) * t191;
t57 = Ifges(6,2) * t109 - t108;
t58 = -Ifges(6,1) * t191 + t240;
t4 = t92 * t54 + t162 * t210 / 0.2e1 + t22 * t76 - t23 * t78 - (-t23 * mrSges(6,3) + t58 / 0.2e1 - t46 / 0.2e1) * t109 - (-t22 * mrSges(6,3) + t48 / 0.2e1 + t57 / 0.2e1) * t191;
t225 = t4 * qJD(1);
t199 = t176 * t68 - t178 * t69;
t215 = t178 * t128;
t216 = t176 * t130;
t221 = t142 * t165;
t5 = t110 * t78 + t113 * t76 + (mrSges(4,3) * t165 + t123 + t56) * t165 + (mrSges(4,3) * t162 - t215 + t216) * t162 + m(6) * (t110 * t22 + t113 * t23 + t165 * t92) + m(5) * (t162 * t199 + t221) + m(4) * (-t144 * t162 + t221) + (m(3) * qJ(2) + mrSges(3,3)) * (t177 ^ 2 + t179 ^ 2);
t224 = t5 * qJD(1);
t12 = m(6) * (t109 * t22 - t191 * t23) - t191 * t76 + t109 * t78 + (-t176 * t128 - t178 * t130 + m(5) * (-t176 * t69 - t178 * t68)) * t165;
t222 = qJD(1) * t12;
t190 = (-t109 * t197 - t164 * t191) * t261;
t202 = m(5) * t208;
t29 = t190 + 0.2e1 * (-t202 / 0.4e1 - m(5) / 0.4e1 - m(6) / 0.4e1) * t165;
t212 = t29 * qJD(1);
t209 = -Ifges(6,5) * t197 - Ifges(6,6) * t164;
t134 = mrSges(6,1) * t164 - mrSges(6,2) * t197;
t206 = t134 * qJD(5);
t36 = (t164 ^ 2 + t197 ^ 2) * mrSges(6,3) + m(6) * (-t141 * t164 - t143 * t197) + (m(5) * qJ(4) + mrSges(5,3)) * t208;
t193 = t252 * t78 + t253 * t76;
t184 = (-t191 * t253 + t268) * mrSges(6,3) + t199 * t263 + (t109 * t141 - t143 * t191 - t164 * t22 - t197 * t23) * t261 - t216 / 0.2e1 + t215 / 0.2e1 + t193;
t7 = (t228 / 0.2e1 + t232 / 0.2e1) * t162 + t184 + t269;
t198 = qJD(1) * t7 + qJD(3) * t36;
t31 = 0.2e1 * mrSges(6,1) * t256 + 0.2e1 * t258 * mrSges(6,2);
t196 = qJD(1) * t31 + qJD(3) * t134;
t187 = (t197 * t258 - t268) * mrSges(6,3) + t193;
t10 = t187 + t266;
t192 = t10 * qJD(1);
t137 = -Ifges(6,2) * t164 - t158;
t138 = -Ifges(6,2) * t197 + t239;
t139 = -Ifges(6,1) * t197 - t239;
t16 = t171 * t134 + (t139 / 0.2e1 - t138 / 0.2e1) * t164 - (t255 + t137 / 0.2e1) * t197;
t183 = -(t48 / 0.4e1 + t57 / 0.4e1) * t197 + (t58 / 0.4e1 - t46 / 0.4e1) * t164 - (-t141 * mrSges(6,3) / 0.2e1 + t140 / 0.4e1 + t137 / 0.4e1) * t191 - (mrSges(6,3) * t254 + t139 / 0.4e1 - t138 / 0.4e1) * t109 + t141 * t76 / 0.2e1 + t78 * t254 + t162 * t209 / 0.4e1 + t171 * t54 / 0.2e1 + t92 * t134 / 0.2e1;
t188 = -Ifges(6,3) * t165 / 0.2e1 - t26 * mrSges(6,1) / 0.2e1 + t27 * mrSges(6,2) / 0.2e1 - t194;
t3 = t183 + t188;
t189 = -qJD(1) * t3 - qJD(3) * t16;
t28 = t190 + (m(5) + m(6) - t202) * t250;
t11 = t187 - t266;
t9 = t165 * t203 + t185 + t186;
t6 = -mrSges(5,2) * t219 / 0.2e1 - mrSges(5,1) * t220 / 0.2e1 + t184 - t269;
t2 = t183 - t188;
t13 = [qJD(2) * t5 + qJD(3) * t1 + qJD(4) * t12 + qJD(5) * t4, t224 + m(6) * (-t110 * t197 + t113 * t164) * qJD(2) + t9 * qJD(3) + t28 * qJD(4) + t11 * qJD(5), t9 * qJD(2) + t6 * qJD(4) + t2 * qJD(5) + t237 + (t178 * qJ(4) * t127 - t176 * qJ(4) * t129 + (-t178 * (Ifges(5,1) * t176 + t241) / 0.2e1 + (Ifges(5,2) * t178 + t242) * t249 - Ifges(4,5)) * t162 + t87 * t248 + t88 * t249 + t47 * t251 + t45 * t253 + t113 * t255 + t138 * t259 + 0.2e1 * (t141 * t26 + t143 * t27 + t171 * t93) * t261 + 0.2e1 * (-pkin(3) * t144 + qJ(4) * t264) * t262 + pkin(3) * t122 + t93 * t136 + t141 * t77 + t142 * mrSges(4,2) + t143 * t75 - t144 * mrSges(4,1) - Ifges(4,6) * t165 + t144 * t167 + t171 * t55 + (Ifges(5,5) * t176 + Ifges(6,5) * t164 + Ifges(5,6) * t178 - Ifges(6,6) * t197) * t250 + (-t164 * t26 - t197 * t27) * mrSges(6,3) + t264 * mrSges(5,3)) * qJD(3), qJD(2) * t28 + qJD(3) * t6 + t222, t225 + t11 * qJD(2) + t2 * qJD(3) + (-mrSges(6,1) * t23 - mrSges(6,2) * t22 + t210) * qJD(5); -qJD(3) * t8 + qJD(4) * t29 + qJD(5) * t10 - t224, 0, -t238, t212, t192 - t206; qJD(2) * t8 + qJD(4) * t7 + qJD(5) * t3 - t237, t238, qJD(4) * t36 + qJD(5) * t16, t198, (-mrSges(6,1) * t143 - mrSges(6,2) * t141 + t209) * qJD(5) - t189; -qJD(2) * t29 - qJD(3) * t7 + qJD(5) * t31 - t222, -t212, -t198 + t206, 0, t196; -qJD(2) * t10 - qJD(3) * t3 - qJD(4) * t31 - t225, -t192, -qJD(4) * t134 + t189, -t196, 0;];
Cq = t13;
