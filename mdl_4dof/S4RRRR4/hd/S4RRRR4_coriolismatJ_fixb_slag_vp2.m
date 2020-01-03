% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:43
% EndTime: 2019-12-31 17:25:47
% DurationCPUTime: 1.81s
% Computational Cost: add. (4635->226), mult. (9714->318), div. (0->0), fcn. (9571->6), ass. (0->153)
t243 = sin(qJ(3));
t244 = sin(qJ(2));
t245 = cos(qJ(3));
t246 = cos(qJ(2));
t128 = t243 * t244 - t245 * t246;
t129 = -t243 * t246 - t244 * t245;
t152 = sin(qJ(4));
t153 = cos(qJ(4));
t186 = Ifges(5,5) * t152 + Ifges(5,6) * t153;
t175 = t129 * t186;
t247 = t153 / 0.2e1;
t249 = t152 / 0.2e1;
t146 = Ifges(5,4) * t153;
t136 = Ifges(5,1) * t152 + t146;
t214 = t153 * t136;
t240 = Ifges(5,4) * t152;
t135 = Ifges(5,2) * t153 + t240;
t215 = t152 * t135;
t261 = -t214 / 0.2e1 + t215 / 0.2e1;
t187 = Ifges(5,2) * t152 - t146;
t56 = -Ifges(5,6) * t129 + t128 * t187;
t188 = Ifges(5,1) * t153 - t240;
t58 = -Ifges(5,5) * t129 - t128 * t188;
t173 = t56 * t247 + t58 * t249 + Ifges(4,6) * t129 - t175 / 0.2e1 + (-Ifges(4,5) + t261) * t128;
t133 = -mrSges(5,1) * t153 + mrSges(5,2) * t152;
t206 = t244 * pkin(5);
t172 = -pkin(6) * t244 - t206;
t149 = t246 * pkin(5);
t213 = pkin(6) * t246 + t149;
t257 = t172 * t243 + t213 * t245;
t262 = t257 * t133;
t265 = t257 * mrSges(4,1);
t106 = -t172 * t245 + t213 * t243;
t273 = t106 * mrSges(4,2);
t275 = t173 + t262 - t265 + t273;
t274 = -t262 / 0.2e1 + t265 / 0.2e1 - t273 / 0.2e1;
t272 = t106 * t153;
t218 = t106 * t257;
t271 = t152 * t106;
t270 = t243 * t106;
t145 = Ifges(5,5) * t153;
t238 = Ifges(5,6) * t152;
t269 = -t145 / 0.2e1 + t238 / 0.2e1;
t248 = -t153 / 0.2e1;
t266 = pkin(3) * t257;
t264 = Ifges(4,4) + t269;
t208 = t245 * pkin(2);
t143 = -t208 - pkin(3);
t263 = t143 * t257;
t144 = -pkin(2) * t246 - pkin(1);
t87 = pkin(3) * t128 + pkin(7) * t129 + t144;
t30 = -t152 * t257 + t153 * t87;
t31 = t152 * t87 + t153 * t257;
t134 = mrSges(5,1) * t152 + mrSges(5,2) * t153;
t85 = t134 * t128;
t86 = t134 * t129;
t231 = t152 * mrSges(5,3);
t89 = mrSges(5,2) * t129 + t128 * t231;
t228 = t153 * mrSges(5,3);
t91 = -mrSges(5,1) * t129 + t128 * t228;
t260 = -t257 * t86 - t106 * t85 + t144 * (-mrSges(4,1) * t129 - mrSges(4,2) * t128) + t30 * t91 + t31 * t89;
t150 = t152 ^ 2;
t151 = t153 ^ 2;
t259 = t150 + t151;
t258 = t187 * t248 + t188 * t249;
t194 = t145 - t238;
t256 = -Ifges(5,3) * t129 / 0.2e1 + (t194 / 0.4e1 + t269) * t128;
t255 = m(5) * pkin(2);
t254 = -mrSges(5,1) / 0.2e1;
t253 = mrSges(5,2) / 0.2e1;
t252 = Ifges(5,3) / 0.2e1;
t251 = pkin(3) * t85;
t250 = -t152 / 0.2e1;
t242 = pkin(3) * t134;
t163 = t58 * t248 + t56 * t249 - t264 * t129 + (-Ifges(5,3) + Ifges(4,1) - Ifges(4,2)) * t128;
t59 = Ifges(5,5) * t128 - t129 * t188;
t227 = t153 * t59;
t57 = Ifges(5,6) * t128 + t129 * t187;
t230 = t152 * t57;
t164 = -t227 / 0.2e1 + t230 / 0.2e1 + t264 * t128;
t207 = t244 * pkin(2);
t94 = -pkin(3) * t129 + pkin(7) * t128;
t88 = t207 + t94;
t38 = t153 * t88 + t271;
t39 = t152 * t88 - t272;
t204 = t129 * t231;
t90 = -mrSges(5,2) * t128 + t204;
t92 = mrSges(5,1) * t128 + t129 * t228;
t1 = -pkin(1) * (mrSges(3,1) * t244 + mrSges(3,2) * t246) + t38 * t92 + t39 * t90 + m(5) * (t30 * t38 + t31 * t39 + t218) + m(4) * t144 * t207 + (mrSges(4,1) * t207 + t164) * t128 + (-mrSges(4,2) * t207 + t163) * t129 + (-Ifges(3,2) + Ifges(3,1)) * t246 * t244 + (-t244 ^ 2 + t246 ^ 2) * Ifges(3,4) + t260;
t237 = t1 * qJD(1);
t232 = t143 * t85;
t229 = t152 * t91;
t226 = t153 * t89;
t46 = t153 * t94 + t271;
t47 = t152 * t94 - t272;
t3 = m(5) * (t30 * t46 + t31 * t47 + t218) + t47 * t90 + t46 * t92 + t164 * t128 + t163 * t129 + t260;
t225 = t3 * qJD(1);
t224 = t38 * t152;
t223 = t39 * t153;
t222 = t46 * t152;
t221 = t47 * t153;
t124 = -t215 / 0.2e1;
t84 = t133 * t129;
t7 = -t106 * t84 + t31 * t92 + (t59 * t250 + t57 * t248 - t128 * t186 / 0.2e1 - t31 * t228 + (t136 * t247 + t124) * t129) * t129 + (-t90 + t204) * t30;
t220 = t7 * qJD(1);
t216 = t143 * t134;
t212 = pkin(7) * t229;
t211 = pkin(7) * t226;
t210 = mrSges(5,3) * t224;
t209 = mrSges(5,3) * t223;
t205 = t243 * pkin(2);
t142 = t205 + pkin(7);
t203 = t142 * t229;
t202 = t142 * t226;
t201 = -Ifges(5,2) / 0.4e1 + Ifges(5,1) / 0.4e1;
t198 = t152 * t245;
t197 = t153 * t245;
t192 = t124 + t214 / 0.2e1 + t258;
t190 = -t198 / 0.2e1;
t189 = mrSges(5,3) * (t151 / 0.2e1 + t150 / 0.2e1);
t185 = t223 - t224;
t184 = t221 - t222;
t158 = (-mrSges(4,1) + t133) * t205 + (mrSges(5,3) * t259 - mrSges(4,2)) * t208;
t171 = t259 * t245;
t20 = (t142 * t171 + t143 * t243) * t255 + t158;
t154 = m(5) * (t263 + t184 * t142 + (t197 * t31 - t198 * t30 + t270) * pkin(2)) / 0.2e1 - t232 / 0.2e1 - t203 / 0.2e1 + t202 / 0.2e1 - t86 * t205 / 0.2e1 + pkin(2) * t92 * t190 + pkin(2) * t90 * t197 / 0.2e1 + (-t222 / 0.2e1 + t221 / 0.2e1) * mrSges(5,3) - t274;
t155 = -m(5) * (pkin(7) * t185 - t266) / 0.2e1 - t251 / 0.2e1 + t212 / 0.2e1 - t211 / 0.2e1 + t210 / 0.2e1 - t209 / 0.2e1 + t274;
t4 = t154 + t155;
t182 = t4 * qJD(1) + t20 * qJD(2);
t25 = t192 + t216;
t159 = (t136 / 0.4e1 + t146 / 0.4e1 + t201 * t152) * t152 + (0.3e1 / 0.4e1 * t240 + t135 / 0.4e1 - t201 * t153) * t153;
t156 = t142 * t189 + t159;
t168 = t106 * t134 / 0.2e1 - t230 / 0.4e1 + t227 / 0.4e1;
t177 = t248 * t92 + t250 * t90;
t160 = t177 * t142 + t143 * t84 / 0.2e1 + t168;
t169 = (-0.3e1 / 0.4e1 * t238 + 0.3e1 / 0.4e1 * t145) * t128;
t180 = t253 * t39 + t254 * t38;
t6 = t169 + (t252 + t156) * t129 + t160 + t180;
t181 = t6 * qJD(1) + t25 * qJD(2);
t179 = t253 * t47 + t254 * t46;
t161 = -t258 + t261;
t165 = (-mrSges(5,2) * t197 / 0.2e1 + mrSges(5,1) * t190) * pkin(2);
t15 = (pkin(3) / 0.2e1 - t143 / 0.2e1) * t134 + t165 + t161;
t40 = t161 + t242;
t157 = pkin(7) * t189 + t159;
t162 = t177 * pkin(7) - pkin(3) * t84 / 0.2e1 + t168;
t9 = t169 + (t252 + t157) * t129 + t162 + t179;
t170 = t9 * qJD(1) - t15 * qJD(2) - t40 * qJD(3);
t16 = -t242 / 0.2e1 + t216 / 0.2e1 + t165 + t192;
t8 = t129 * t157 + t162 - t179 + t256;
t5 = t129 * t156 + t160 - t180 + t256;
t2 = -t155 + t154 + t173;
t10 = [qJD(2) * t1 + qJD(3) * t3 - qJD(4) * t7, t237 + (t202 + t209 - t232 + m(5) * (t142 * t185 + t263) - t203 - t210 + mrSges(3,2) * t206 - mrSges(3,1) * t149 + Ifges(3,5) * t246 - Ifges(3,6) * t244 + m(4) * (-t245 * t257 - t270) * pkin(2) + (t128 * t208 + t129 * t205) * mrSges(4,3) + t275) * qJD(2) + t2 * qJD(3) + t5 * qJD(4), t225 + t2 * qJD(2) + (t251 + m(5) * (pkin(7) * t184 - t266) + t211 - t212 + t184 * mrSges(5,3) + t275) * qJD(3) + t8 * qJD(4), -t220 + t5 * qJD(2) + t8 * qJD(3) + (-t31 * mrSges(5,1) - t30 * mrSges(5,2) + t175) * qJD(4); qJD(3) * t4 + qJD(4) * t6 - t237, qJD(3) * t20 + qJD(4) * t25, ((-pkin(3) * t243 + pkin(7) * t171) * t255 + t158) * qJD(3) + t16 * qJD(4) + t182, t16 * qJD(3) + (t133 * t142 + t194) * qJD(4) + t181; -qJD(2) * t4 + qJD(4) * t9 - t225, -qJD(4) * t15 - t182, -t40 * qJD(4), (pkin(7) * t133 + t194) * qJD(4) + t170; -qJD(2) * t6 - qJD(3) * t9 + t220, qJD(3) * t15 - t181, -t170, 0;];
Cq = t10;
