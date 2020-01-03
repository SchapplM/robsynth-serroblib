% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP3
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:09
% EndTime: 2019-12-31 21:49:13
% DurationCPUTime: 1.79s
% Computational Cost: add. (3044->242), mult. (6711->308), div. (0->0), fcn. (4516->6), ass. (0->157)
t136 = sin(qJ(4));
t134 = t136 ^ 2;
t139 = cos(qJ(4));
t135 = t139 ^ 2;
t201 = t134 + t135;
t186 = (-Ifges(5,4) + Ifges(6,5)) * t139;
t260 = mrSges(5,3) + mrSges(6,2);
t266 = Ifges(5,1) + Ifges(6,1);
t265 = t260 * t201;
t140 = cos(qJ(3));
t229 = t140 * pkin(2);
t264 = t201 * t229;
t137 = sin(qJ(3));
t230 = t137 * pkin(2);
t124 = pkin(8) + t230;
t263 = t201 * t124;
t138 = sin(qJ(2));
t121 = t137 * t138 * pkin(1);
t141 = cos(qJ(2));
t231 = pkin(1) * t141;
t126 = pkin(2) + t231;
t86 = t126 * t140 - t121;
t262 = t201 * t86;
t261 = -Ifges(6,3) + t266;
t215 = t139 * mrSges(5,2);
t217 = t136 * mrSges(5,1);
t173 = t215 + t217;
t81 = -pkin(3) - t86;
t160 = t81 * t173;
t214 = t139 * mrSges(6,3);
t216 = t136 * mrSges(6,1);
t172 = -t214 + t216;
t209 = qJ(5) * t136;
t171 = -pkin(4) * t139 - t209;
t105 = -pkin(3) + t171;
t46 = t105 - t86;
t246 = t46 / 0.2e1;
t259 = t172 * t246 + t160 / 0.2e1;
t125 = -pkin(3) - t229;
t159 = t125 * t173;
t91 = t105 - t229;
t242 = t91 / 0.2e1;
t258 = t172 * t242 + t159 / 0.2e1;
t107 = -t139 * mrSges(6,1) - t136 * mrSges(6,3);
t108 = -t139 * mrSges(5,1) + t136 * mrSges(5,2);
t257 = t108 + t107;
t154 = -mrSges(4,2) + t265;
t191 = -mrSges(4,1) + t257;
t256 = t154 * t229 + t191 * t230;
t204 = t138 * t140;
t87 = pkin(1) * t204 + t126 * t137;
t255 = t154 * t86 + t191 * t87;
t92 = (t137 * t141 + t204) * pkin(1);
t93 = t140 * t231 - t121;
t253 = t154 * t93 + (-mrSges(3,1) * t138 - mrSges(3,2) * t141) * pkin(1) + t191 * t92;
t252 = 2 * qJD(3);
t251 = m(5) / 0.2e1;
t250 = -m(6) / 0.2e1;
t249 = m(6) / 0.2e1;
t248 = pkin(3) / 0.2e1;
t247 = -mrSges(4,2) / 0.2e1;
t245 = -t81 / 0.2e1;
t244 = -t86 / 0.2e1;
t243 = t86 / 0.2e1;
t240 = -t93 / 0.2e1;
t239 = t93 / 0.2e1;
t238 = t105 / 0.2e1;
t237 = -t125 / 0.2e1;
t236 = t136 / 0.2e1;
t234 = m(6) * qJ(5);
t109 = pkin(4) * t136 - qJ(5) * t139;
t233 = m(6) * t109;
t232 = m(6) * t139;
t82 = pkin(8) + t87;
t228 = t262 * t82;
t218 = t135 * t93;
t220 = t134 * t93;
t227 = (t218 + t220) * t82;
t226 = t46 + t91;
t225 = t263 * t93;
t224 = t262 * pkin(8);
t223 = Ifges(5,4) * t136;
t128 = t139 * mrSges(6,2);
t3 = m(5) * (t81 * t87 + t228) + m(6) * (t46 * t87 + t228) + t255;
t213 = t3 * qJD(1);
t4 = m(5) * (t81 * t92 + t227) + m(6) * (t46 * t92 + t227) + m(4) * (-t86 * t92 + t87 * t93) + t253;
t212 = t4 * qJD(1);
t211 = t105 + t46;
t210 = t105 + t91;
t32 = (m(6) * t46 + t107) * t136;
t208 = qJD(1) * t32;
t130 = Ifges(6,5) * t136;
t80 = t109 * t107;
t155 = t80 - t136 * t223 / 0.2e1 + (0.2e1 * t130 - t223) * t236 + (t261 * t236 + (-Ifges(5,2) - Ifges(6,3) / 0.2e1 + t266 / 0.2e1) * t136 - t186) * t139;
t158 = t172 + t233;
t11 = t158 * t46 + t155 + t160;
t207 = t11 * qJD(1);
t203 = t264 * t124;
t202 = t264 * pkin(8);
t127 = mrSges(6,3) + t234;
t200 = t127 * qJD(4);
t199 = pkin(4) * t250;
t196 = t234 / 0.2e1;
t195 = m(6) * t236;
t194 = -t229 / 0.2e1;
t193 = t229 / 0.2e1;
t192 = t244 + t239;
t187 = t263 * t86 + t264 * t82;
t185 = m(6) * t193;
t184 = t226 * t233;
t183 = t211 * t233;
t182 = t210 * t233;
t181 = t201 * t93 * pkin(8);
t180 = t244 + t245 + t248;
t177 = t239 + t246 + t242;
t176 = t243 + t246 + t238;
t175 = t240 + t245 + t237;
t174 = mrSges(4,1) / 0.2e1 - t107 / 0.2e1 - t108 / 0.2e1;
t12 = t158 * t91 + t155 + t159;
t153 = t223 + (Ifges(5,2) - t261) * t139 - t130;
t7 = -t80 - t184 / 0.2e1 + (mrSges(5,2) * t175 + mrSges(6,3) * t177 + t196 * t93 + t186) * t139 + (mrSges(5,1) * t175 - mrSges(6,1) * t177 + t199 * t93 + t153) * t136;
t170 = -t7 * qJD(1) + t12 * qJD(2);
t15 = m(6) * (t230 * t91 + t203) + m(5) * (t125 * t230 + t203) + t256;
t142 = -t174 * t87 + (t125 * t87 + t230 * t81 + t187) * t251 + (t230 * t46 + t87 * t91 + t187) * t249;
t144 = -m(5) * (-pkin(3) * t92 + t181) / 0.2e1 + (t105 * t92 + t181) * t250;
t162 = t174 * t137;
t2 = t192 * mrSges(4,2) + t174 * t92 + (t140 * t247 - t162) * pkin(2) + t142 + t144 + (t193 - t192) * t265;
t169 = t2 * qJD(1) + t15 * qJD(2);
t17 = (m(6) * t177 + t107) * t136;
t40 = (m(6) * t91 + t107) * t136;
t167 = qJD(1) * t17 + qJD(2) * t40;
t166 = t194 + t237 + t248;
t164 = t193 + t242 + t238;
t163 = pkin(3) * t173;
t13 = t105 * t158 + t155 - t163;
t5 = -t80 - t183 / 0.2e1 + (mrSges(5,2) * t180 + mrSges(6,3) * t176 + t196 * t86 + t186) * t139 + (mrSges(5,1) * t180 - mrSges(6,1) * t176 + t199 * t86 + t153) * t136;
t9 = -t80 - t182 / 0.2e1 + (mrSges(5,2) * t166 + mrSges(6,3) * t164 + qJ(5) * t185 + t186) * t139 + (m(6) * pkin(4) * t194 + mrSges(5,1) * t166 - mrSges(6,1) * t164 + t153) * t136;
t157 = t5 * qJD(1) + t9 * qJD(2) - t13 * qJD(3);
t19 = (m(6) * t176 + t107) * t136;
t25 = (m(6) * t164 + t107) * t136;
t44 = (m(6) * t105 + t107) * t136;
t156 = qJD(1) * t19 + qJD(2) * t25 + qJD(3) * t44;
t152 = t172 * t238 + t155 - t163 / 0.2e1;
t147 = (-mrSges(6,2) * t209 - pkin(4) * t128 + (Ifges(6,4) + Ifges(5,5)) * t139 + (-Ifges(5,6) + Ifges(6,6)) * t136) * qJD(4);
t145 = qJD(4) * (m(6) * t171 + t257);
t143 = -t215 / 0.2e1 - t217 / 0.2e1 - t216 / 0.2e1 + t214 / 0.2e1 - t233 / 0.2e1;
t106 = pkin(8) * t232 + t128;
t88 = t124 * t232 + t128;
t45 = t232 * t82 + t128;
t26 = -t210 * t195 + (-t107 + t185) * t136;
t20 = -t211 * t195 + (m(6) * t243 - t107) * t136;
t18 = -t226 * t195 + (m(6) * t239 - t107) * t136;
t10 = t143 * t229 + t182 / 0.2e1 + t152 + t258;
t8 = t143 * t93 + t184 / 0.2e1 + t155 + t258 + t259;
t6 = t143 * t86 + t183 / 0.2e1 + t152 + t259;
t1 = mrSges(4,2) * t240 - pkin(2) * t162 + t142 - t144 + (t86 + t229) * (t247 + t260 * (t134 / 0.2e1 + t135 / 0.2e1)) + (-mrSges(4,1) / 0.2e1 + t257 / 0.2e1) * t92 + t260 * (t218 / 0.2e1 + t220 / 0.2e1);
t14 = [qJD(2) * t4 + qJD(3) * t3 + qJD(4) * t11 - qJD(5) * t32, t1 * qJD(3) + t8 * qJD(4) + t18 * qJD(5) + t212 + (0.2e1 * (t125 * t92 + t225) * t251 + 0.2e1 * (t91 * t92 + t225) * t249 + m(4) * (t137 * t93 - t140 * t92) * pkin(2) + t253) * qJD(2), t213 + t1 * qJD(2) + t6 * qJD(4) + t20 * qJD(5) + ((-pkin(3) * t87 + t224) * t251 + (t105 * t87 + t224) * t249) * t252 + t255 * qJD(3), t8 * qJD(2) + t6 * qJD(3) + t45 * qJD(5) + t145 * t82 + t147 + t207, qJD(2) * t18 + qJD(3) * t20 + qJD(4) * t45 - t208; qJD(3) * t2 - qJD(4) * t7 - qJD(5) * t17 - t212, qJD(3) * t15 + qJD(4) * t12 - qJD(5) * t40, t10 * qJD(4) + t26 * qJD(5) + ((t105 * t230 + t202) * t249 + (-pkin(3) * t230 + t202) * t251) * t252 + t169 + t256 * qJD(3), t10 * qJD(3) + t88 * qJD(5) + t124 * t145 + t147 + t170, qJD(3) * t26 + qJD(4) * t88 - t167; -qJD(2) * t2 - qJD(4) * t5 - qJD(5) * t19 - t213, -qJD(4) * t9 - qJD(5) * t25 - t169, qJD(4) * t13 - qJD(5) * t44, pkin(8) * t145 + t106 * qJD(5) + t147 - t157, qJD(4) * t106 - t156; qJD(2) * t7 + qJD(3) * t5 - t207, qJD(3) * t9 - t170, t157, t127 * qJD(5), t200; qJD(2) * t17 + qJD(3) * t19 + t208, qJD(3) * t25 + t167, t156, -t200, 0;];
Cq = t14;
