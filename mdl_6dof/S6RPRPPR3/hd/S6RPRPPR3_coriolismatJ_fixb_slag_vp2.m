% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:09
% EndTime: 2019-03-09 02:44:12
% DurationCPUTime: 1.86s
% Computational Cost: add. (3207->297), mult. (5891->392), div. (0->0), fcn. (4534->6), ass. (0->176)
t150 = sin(qJ(6));
t146 = t150 ^ 2;
t152 = cos(qJ(6));
t147 = t152 ^ 2;
t201 = t146 + t147;
t276 = t201 * mrSges(7,3);
t225 = t152 * mrSges(7,2);
t232 = t150 * mrSges(7,1);
t275 = -t225 / 0.2e1 - t232 / 0.2e1;
t181 = t225 + t232;
t153 = cos(qJ(3));
t129 = sin(pkin(9)) * pkin(1) + pkin(7);
t268 = -qJ(5) + t129;
t89 = t268 * t153;
t273 = t89 * t181;
t138 = t152 * mrSges(7,1);
t231 = t150 * mrSges(7,2);
t267 = -t138 + t231;
t272 = t267 - mrSges(6,1);
t149 = qJ(4) + pkin(5);
t271 = t149 * t153;
t270 = t149 * t181;
t239 = Ifges(7,4) * t152;
t179 = Ifges(7,1) * t150 + t239;
t248 = -t152 / 0.2e1;
t269 = t179 * t248;
t246 = t153 / 0.2e1;
t266 = t246 * t276;
t209 = t152 * t153;
t151 = sin(qJ(3));
t227 = t151 * mrSges(7,1);
t102 = mrSges(7,3) * t209 + t227;
t210 = t152 * t102;
t212 = t150 * t153;
t226 = t151 * mrSges(7,2);
t100 = mrSges(7,3) * t212 - t226;
t214 = t150 * t100;
t265 = t214 + t210;
t235 = Ifges(7,6) * t150;
t237 = Ifges(7,5) * t152;
t167 = t237 / 0.2e1 - t235 / 0.2e1;
t264 = Ifges(5,5) - Ifges(4,4) - Ifges(6,4) + t167;
t263 = 2 * m(7);
t262 = 2 * qJD(3);
t261 = m(6) / 0.2e1;
t260 = -m(7) / 0.2e1;
t259 = m(7) / 0.2e1;
t258 = t89 / 0.2e1;
t257 = m(5) + m(6);
t256 = -pkin(4) - pkin(8);
t142 = t151 * pkin(3);
t208 = t153 * qJ(4);
t108 = -t142 + t208;
t97 = -pkin(4) * t151 + t108;
t255 = m(6) * t97;
t188 = t201 * t151;
t254 = m(7) * (-t151 + t188) * t153;
t145 = -pkin(3) + t256;
t252 = t145 / 0.2e1;
t251 = -t150 / 0.2e1;
t250 = t150 / 0.2e1;
t249 = t151 / 0.2e1;
t247 = t152 / 0.2e1;
t245 = m(7) * t151;
t244 = -mrSges(4,1) - mrSges(5,1);
t243 = mrSges(4,2) - mrSges(5,3);
t242 = -cos(pkin(9)) * pkin(1) - pkin(2);
t241 = mrSges(7,3) * t151;
t240 = Ifges(7,4) * t150;
t238 = Ifges(7,5) * t151;
t236 = Ifges(7,2) * t152;
t234 = Ifges(7,6) * t151;
t87 = t267 * t153;
t233 = t149 * t87;
t49 = t256 * t151 - t142 + t271;
t32 = -t150 * t89 + t152 * t49;
t230 = t150 * t32;
t33 = t150 * t49 + t152 * t89;
t229 = t150 * t33;
t178 = -Ifges(7,2) * t150 + t239;
t159 = t178 * t153;
t67 = -t159 + t234;
t228 = t150 * t67;
t224 = t152 * t32;
t223 = t152 * t33;
t128 = Ifges(7,4) * t212;
t69 = -Ifges(7,1) * t209 + t128 + t238;
t222 = t152 * t69;
t221 = t153 * mrSges(7,1);
t220 = t153 * mrSges(7,2);
t219 = t153 * t89;
t154 = -pkin(3) - pkin(4);
t105 = t151 * t154 + t208;
t190 = t146 / 0.2e1 + t147 / 0.2e1;
t155 = t190 * t241 + t246 * t267;
t101 = -t152 * t241 + t221;
t99 = -t150 * t241 - t220;
t164 = t101 * t248 + t99 * t251;
t202 = t153 * mrSges(6,1) + t151 * mrSges(6,2);
t43 = t145 * t188 + t271;
t9 = (-t43 / 0.4e1 - t229 / 0.4e1 - t224 / 0.4e1) * t263 + 0.2e1 * (-t105 / 0.4e1 - t97 / 0.4e1) * m(6) + t155 + t164 - t202;
t218 = t9 * qJD(1);
t109 = mrSges(6,1) * t151 - mrSges(6,2) * t153;
t86 = -t153 * pkin(3) - t151 * qJ(4) + t242;
t189 = -m(5) * t86 + t153 * mrSges(5,1) + t151 * mrSges(5,3);
t65 = t153 * pkin(4) - t86;
t41 = pkin(5) * t151 + pkin(8) * t153 + t65;
t88 = t268 * t151;
t28 = -t150 * t88 + t152 * t41;
t29 = t150 * t41 + t152 * t88;
t11 = (m(7) * (t150 * t29 + t152 * t28) + m(6) * t65 + t109 + t189 + t265) * t151;
t217 = qJD(1) * t11;
t163 = t214 / 0.2e1 + t210 / 0.2e1;
t183 = mrSges(7,3) * t190;
t171 = t153 * t183;
t13 = t87 * t249 + (-t171 + t163) * t153;
t216 = qJD(1) * t13;
t193 = t227 / 0.2e1;
t14 = (t193 + t102 / 0.2e1) * t152 + (-t226 / 0.2e1 + t100 / 0.2e1) * t150 - t171;
t215 = t14 * qJD(1);
t213 = t150 * t102;
t211 = t152 * t100;
t162 = t213 / 0.2e1 - t211 / 0.2e1;
t191 = t225 / 0.2e1;
t205 = t150 * t193 + t151 * t191;
t17 = t162 + t205;
t207 = t17 * qJD(1);
t35 = (t190 * t263 + m(6)) * t151;
t206 = t35 * qJD(1);
t204 = t275 * t151;
t203 = Ifges(7,5) * t212 + Ifges(7,6) * t209;
t200 = qJD(3) * t151;
t199 = qJD(3) * t153;
t198 = -t254 / 0.2e1;
t197 = t254 / 0.2e1;
t196 = t245 / 0.2e1;
t185 = m(5) * t129 + mrSges(5,2) - mrSges(6,3);
t184 = t201 * t259;
t180 = Ifges(7,1) * t152 - t240;
t177 = -t236 - t240;
t176 = t235 - t237;
t66 = Ifges(7,6) * t153 + t178 * t151;
t68 = Ifges(7,5) * t153 + t180 * t151;
t1 = t29 * t99 + t33 * t100 + t28 * t101 + t32 * t102 + t97 * t109 + t189 * t108 + m(7) * (t28 * t32 + t29 * t33 - t88 * t89) + (t242 * mrSges(4,2) - t86 * mrSges(5,3) - t264 * t153 + t88 * t181 + t68 * t248 + t66 * t250) * t153 + (-t228 / 0.2e1 + t222 / 0.2e1 + t273 + t242 * mrSges(4,1) + t86 * mrSges(5,1) + (Ifges(6,2) + Ifges(5,1) - Ifges(5,3) + Ifges(4,1) - Ifges(4,2) - Ifges(6,1) + Ifges(7,3)) * t153 + t264 * t151) * t151 + (t202 + t255) * t65;
t165 = t101 * t251 + t99 * t247;
t166 = t191 + t232 / 0.2e1;
t173 = t223 - t230;
t174 = t150 * t28 - t152 * t29;
t6 = ((-t173 + t89) * t260 + t166 * t153 + t165) * t153 + ((-t174 - t88) * t260 - t166 * t151 + t162) * t151;
t175 = t1 * qJD(1) - t6 * qJD(2);
t5 = t28 * t100 - t29 * t102 + t89 * t87 + t203 * t249 + (t153 * t269 - t174 * mrSges(7,3) + t67 * t247 + (t236 * t153 + t128 + t69) * t250) * t153;
t172 = t5 * qJD(1) + t13 * qJD(2);
t170 = t105 * t261 + t43 * t259;
t169 = t32 * mrSges(7,1) / 0.2e1 - t33 * mrSges(7,2) / 0.2e1;
t168 = -t6 * qJD(1) - qJD(2) * t254;
t8 = (t151 * mrSges(6,3) - t211 + t213) * t151 + (mrSges(6,3) + t181) * t153 ^ 2 + m(7) * (t174 * t151 - t219) + m(6) * (-t151 * t88 - t219);
t161 = -t8 * qJD(1) + qJD(2) * t198;
t160 = t151 * t181;
t25 = -t147 * Ifges(7,4) + t270 + (t240 + (-Ifges(7,1) + Ifges(7,2)) * t152) * t150;
t3 = -t233 / 0.2e1 + (Ifges(7,3) / 0.2e1 - t145 * t183) * t153 + (-0.3e1 / 0.4e1 * t234 + mrSges(7,1) * t258 - t67 / 0.4e1 + t100 * t252 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.4e1) * t212) * t150 + (0.3e1 / 0.4e1 * t238 + mrSges(7,2) * t258 + t69 / 0.4e1 + t128 / 0.4e1 + t102 * t252 + (0.5e1 / 0.4e1 * t240 + (Ifges(7,2) / 0.2e1 - Ifges(7,1) / 0.4e1) * t152) * t153) * t152 + t169;
t39 = t160 / 0.2e1 + t204;
t157 = t3 * qJD(1) + t39 * qJD(2) + t25 * qJD(3);
t16 = (t99 / 0.2e1 + t220 / 0.2e1) * t152 + (-t101 / 0.2e1 + t221 / 0.2e1) * t150 + (-t230 / 0.4e1 + t223 / 0.4e1 - t89 / 0.4e1) * t263;
t44 = m(7) * t149 + t257 * qJ(4) + mrSges(5,3) - t272;
t48 = (-0.1e1 / 0.2e1 + t190) * t245;
t156 = qJD(1) * t16 + qJD(2) * t48 - qJD(3) * t44;
t40 = -t160 / 0.2e1 + t204;
t37 = t196 + (t184 + t257) * t151;
t36 = t151 * t184 - t201 * t196;
t18 = -t162 + t205;
t15 = (-t231 / 0.2e1 + t138 / 0.2e1) * t151 - t163 + t266;
t12 = m(6) * t89 + m(7) * t258 + t173 * t259 + t165 + (t185 + t275) * t153;
t10 = (t224 + t229) * t259 + t255 / 0.2e1 + t155 - t164 - t170;
t4 = -t273 / 0.2e1 - t222 / 0.4e1 - t179 * t212 / 0.2e1 + t228 / 0.4e1 - t152 * (Ifges(7,2) * t209 + t128) / 0.4e1 + t233 / 0.2e1 - t150 * t159 / 0.4e1 + Ifges(7,3) * t246 + t169 + (t176 / 0.4e1 + t167) * t151 + (t180 + t177) * t209 / 0.4e1 + (-t265 / 0.2e1 + t266) * t145;
t2 = -t6 * qJD(3) + qJD(5) * t197 + t13 * qJD(6);
t7 = [qJD(3) * t1 + qJD(4) * t11 + qJD(5) * t8 + qJD(6) * t5, t2, t12 * qJD(4) + t10 * qJD(5) + t4 * qJD(6) + ((t173 * t145 - t149 * t88) * t259 + (-qJ(4) * t88 + t154 * t89) * t261) * t262 + (-pkin(3) * mrSges(5,2) - t154 * mrSges(6,3) + Ifges(6,6) + Ifges(4,5) + Ifges(5,4) + Ifges(7,5) * t251 + Ifges(7,6) * t248 + (-m(5) * pkin(3) + t244) * t129) * t199 + (-t185 * qJ(4) + t243 * t129 + t177 * t251 - Ifges(6,5) - Ifges(4,6) + Ifges(5,6) + t269 + t270) * t200 + t175 + (t89 * mrSges(6,2) + (-t33 * mrSges(7,3) + t145 * t99 - t66 / 0.2e1) * t152 + (t32 * mrSges(7,3) - t145 * t101 - t68 / 0.2e1) * t150 + t272 * t88) * qJD(3), qJD(3) * t12 + qJD(5) * t36 + qJD(6) * t15 + t217, t10 * qJD(3) + t36 * qJD(4) + t18 * qJD(6) - t161, t4 * qJD(3) + t15 * qJD(4) + t18 * qJD(5) + (-mrSges(7,1) * t29 - mrSges(7,2) * t28 + t203) * qJD(6) + t172; t2, -qJD(3) * t254, t202 * qJD(3) + t37 * qJD(4) + t40 * qJD(6) + (-t267 - t243) * t199 + (m(5) * t108 / 0.2e1 + t170) * t262 + (t244 - t276) * t200 + t168, t37 * qJD(3), qJD(1) * t197, qJD(3) * t40 - qJD(6) * t87 + t216; -qJD(4) * t16 + qJD(5) * t9 - qJD(6) * t3 - t175, -t48 * qJD(4) - t39 * qJD(6) - t168, qJD(4) * t44 - qJD(6) * t25, -t156, t218 (t145 * t267 + t176) * qJD(6) - t157; qJD(3) * t16 - qJD(5) * t35 - qJD(6) * t14 - t217, t48 * qJD(3), t156, 0, -t206, qJD(6) * t267 - t215; -t9 * qJD(3) + t35 * qJD(4) - t17 * qJD(6) + t161, qJD(1) * t198, -t218, t206, 0, -qJD(6) * t181 - t207; qJD(3) * t3 + qJD(4) * t14 + qJD(5) * t17 - t172, qJD(3) * t39 - t216, t157, t215, t207, 0;];
Cq  = t7;
