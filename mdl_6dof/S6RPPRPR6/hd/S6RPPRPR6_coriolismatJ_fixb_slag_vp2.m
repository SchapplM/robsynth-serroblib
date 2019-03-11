% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:44
% EndTime: 2019-03-09 01:50:48
% DurationCPUTime: 1.87s
% Computational Cost: add. (2722->295), mult. (5103->395), div. (0->0), fcn. (3666->4), ass. (0->168)
t152 = sin(qJ(6));
t146 = t152 ^ 2;
t154 = cos(qJ(6));
t148 = t154 ^ 2;
t201 = -t148 - t146;
t268 = t201 * mrSges(7,3);
t137 = t152 * mrSges(7,1);
t140 = t154 * mrSges(7,2);
t184 = t137 + t140;
t267 = qJD(6) * t184;
t266 = -t201 / 0.2e1;
t153 = sin(qJ(4));
t243 = m(7) * t153;
t150 = qJ(2) - pkin(7);
t264 = -pkin(5) + t150;
t195 = t146 / 0.2e1 + t148 / 0.2e1;
t189 = mrSges(7,3) * t195;
t238 = Ifges(7,6) * t154;
t240 = Ifges(7,5) * t152;
t171 = t238 / 0.2e1 + t240 / 0.2e1;
t263 = Ifges(6,6) + Ifges(5,4) - t171;
t252 = pkin(4) + pkin(8);
t262 = t252 * t201;
t155 = cos(qJ(4));
t116 = -t153 * mrSges(6,2) - t155 * mrSges(6,3);
t136 = t155 * qJ(5);
t202 = -pkin(4) * t153 + t136;
t247 = -t153 / 0.2e1;
t261 = t247 * t268;
t109 = t264 * t153;
t220 = qJ(5) * t153;
t98 = t155 * t252 + t220;
t38 = t154 * t109 - t152 * t98;
t229 = t154 * t38;
t39 = t109 * t152 + t154 * t98;
t234 = t152 * t39;
t176 = t229 + t234;
t260 = t137 / 0.2e1 + t140 / 0.2e1;
t151 = pkin(1) + qJ(3);
t107 = -t202 + t151;
t212 = t153 * t154;
t225 = t155 * mrSges(7,2);
t105 = mrSges(7,3) * t212 - t225;
t209 = t154 * t105;
t213 = t152 * t153;
t227 = t155 * mrSges(7,1);
t103 = -mrSges(7,3) * t213 + t227;
t216 = t152 * t103;
t110 = t264 * t155;
t85 = pkin(8) * t153 + t107;
t36 = -t110 * t154 - t152 * t85;
t37 = -t110 * t152 + t154 * t85;
t259 = -m(6) * t107 + m(7) * (t152 * t36 - t154 * t37) - t116 - t209 + t216;
t258 = 0.2e1 * m(7);
t257 = m(5) / 0.4e1;
t256 = m(6) / 0.2e1;
t255 = m(7) / 0.2e1;
t254 = -mrSges(7,1) / 0.2e1;
t253 = -Ifges(7,2) / 0.2e1;
t251 = (0.1e1 + t201) * t155 * t243;
t250 = t109 / 0.2e1;
t249 = -t152 / 0.2e1;
t248 = t152 / 0.2e1;
t246 = t153 / 0.2e1;
t245 = -t154 / 0.2e1;
t244 = t154 / 0.2e1;
t242 = Ifges(7,4) * t152;
t241 = Ifges(7,4) * t154;
t239 = Ifges(7,5) * t155;
t237 = Ifges(7,6) * t155;
t88 = t184 * t153;
t236 = qJ(5) * t88;
t235 = t152 * mrSges(7,2);
t131 = Ifges(7,4) * t212;
t65 = Ifges(7,1) * t213 + t131 + t239;
t233 = t152 * t65;
t232 = t153 * mrSges(7,1);
t231 = t153 * mrSges(7,2);
t230 = t154 * mrSges(7,1);
t180 = Ifges(7,2) * t154 + t242;
t63 = t153 * t180 + t237;
t228 = t154 * t63;
t226 = t155 * mrSges(5,2);
t224 = t155 * mrSges(7,3);
t130 = Ifges(7,5) * t212;
t4 = t36 * t105 - t37 * t103 + t155 * t130 / 0.2e1 + t109 * t88 + ((-t36 * mrSges(7,3) + t131 / 0.2e1 + t65 / 0.2e1) * t154 + (-t37 * mrSges(7,3) - t237 / 0.2e1 - t63 / 0.2e1 + (-t242 / 0.2e1 + (t253 + Ifges(7,1) / 0.2e1) * t154) * t153) * t152) * t153;
t223 = t4 * qJD(1);
t178 = -t152 * t37 - t154 * t36;
t147 = t153 ^ 2;
t149 = t155 ^ 2;
t200 = t149 + t147;
t211 = t154 * t103;
t215 = t152 * t105;
t185 = t230 - t235;
t86 = t153 * t185;
t7 = -t153 * t86 + mrSges(4,2) + mrSges(3,3) + (-t211 - t215) * t155 + (m(4) + m(3)) * qJ(2) + m(7) * (t153 * t109 + t155 * t178) + (-mrSges(5,3) - mrSges(6,1) + 0.4e1 * (m(6) / 0.4e1 + t257) * t150) * t200;
t222 = t7 * qJD(1);
t142 = t155 * mrSges(5,1);
t104 = -t152 * t224 - t232;
t106 = t154 * t224 + t231;
t161 = (t152 * t38 - t154 * t39) * t255 + t104 * t248 + t106 * t245;
t162 = m(7) * (-t155 * t262 + t220);
t117 = pkin(4) * t155 + t220;
t194 = m(6) * t117 + mrSges(6,3) * t153;
t8 = -t142 + (mrSges(5,2) - t184 / 0.2e1) * t153 + (mrSges(6,2) - t189) * t155 - t162 / 0.2e1 + t161 - t194;
t221 = t8 * qJD(1);
t169 = t216 / 0.2e1 - t209 / 0.2e1;
t174 = t153 * t189;
t159 = (t174 + t169) * t155 + t88 * t246;
t11 = t159 + t260;
t219 = t11 * qJD(1);
t186 = -t153 * mrSges(5,1) - t226;
t13 = mrSges(4,3) + 0.4e1 * (t257 + m(4) / 0.4e1) * t151 - t186 - t259;
t218 = t13 * qJD(1);
t15 = t259 * t155;
t217 = t15 * qJD(1);
t214 = t152 * t106;
t210 = t154 * t104;
t208 = t252 * t103;
t207 = t252 * t105;
t187 = t225 / 0.2e1 - t105 / 0.2e1;
t188 = t227 / 0.2e1 + t103 / 0.2e1;
t16 = t152 * t188 + t154 * t187 + t174;
t206 = t16 * qJD(1);
t19 = t152 * t187 - t154 * t188;
t205 = t19 * qJD(1);
t160 = (-t149 * t201 + t147) * t255 + (m(5) / 0.2e1 + t256) * t200;
t190 = t201 * t255;
t165 = -m(5) / 0.2e1 - m(6) / 0.2e1 + t190;
t28 = -m(4) - t160 + t165;
t204 = t28 * qJD(1);
t32 = (t195 * t258 + m(6)) * t155;
t203 = t32 * qJD(1);
t199 = qJD(4) * t153;
t191 = m(7) * t266;
t183 = Ifges(7,1) * t154 - t242;
t182 = Ifges(7,1) * t152 + t241;
t181 = -Ifges(7,2) * t152 + t241;
t64 = -Ifges(7,6) * t153 + t155 * t180;
t66 = -Ifges(7,5) * t153 + t155 * t182;
t87 = t185 * t155;
t1 = m(7) * (t109 * t110 + t36 * t38 + t37 * t39) + t38 * t103 + t36 * t104 + t39 * t105 + t37 * t106 - t109 * t87 - t110 * t86 + t117 * t116 + t151 * t142 + t194 * t107 + (t228 / 0.2e1 + t233 / 0.2e1 - t107 * mrSges(6,2) - t263 * t155) * t155 + (-t151 * mrSges(5,2) + t64 * t244 + t66 * t248 + (-Ifges(5,1) + Ifges(5,2) + Ifges(6,3) - Ifges(6,2) - Ifges(7,3)) * t155 + t263 * t153) * t153;
t167 = -t214 / 0.2e1 - t210 / 0.2e1;
t168 = t215 / 0.2e1 + t211 / 0.2e1;
t6 = ((t109 - t176) * t255 - t86 / 0.2e1 + t167) * t155 + ((t110 - t178) * t255 - t87 / 0.2e1 + t168) * t153;
t179 = t1 * qJD(1) + t6 * qJD(3);
t173 = t38 * mrSges(7,1) / 0.2e1 - t39 * mrSges(7,2) / 0.2e1;
t2 = -t236 / 0.2e1 + (-Ifges(7,3) / 0.2e1 - t252 * t189) * t153 + (0.3e1 / 0.4e1 * t237 + t109 * t254 + t63 / 0.4e1 + t207 / 0.2e1 + (-Ifges(7,1) / 0.2e1 + Ifges(7,2) / 0.4e1) * t212) * t154 + (0.3e1 / 0.4e1 * t239 + mrSges(7,2) * t250 + t65 / 0.4e1 + t131 / 0.4e1 - t208 / 0.2e1 + (0.5e1 / 0.4e1 * t241 + (t253 + Ifges(7,1) / 0.4e1) * t152) * t153) * t152 + t173;
t27 = t148 * Ifges(7,4) - qJ(5) * t185 + (-t242 + (Ifges(7,1) - Ifges(7,2)) * t154) * t152;
t175 = t2 * qJD(1) + t27 * qJD(4);
t172 = t6 * qJD(1) + qJD(3) * t251;
t170 = t230 / 0.2e1 - t235 / 0.2e1;
t18 = (t104 / 0.2e1 + t232 / 0.2e1) * t154 + (t106 / 0.2e1 - t231 / 0.2e1) * t152 + (t234 / 0.4e1 + t229 / 0.4e1 - t109 / 0.4e1) * t258;
t48 = (-0.1e1 / 0.2e1 + t195) * t243;
t83 = mrSges(6,3) + (m(6) + m(7)) * qJ(5) + t184;
t163 = qJD(1) * t18 + qJD(3) * t48 - qJD(4) * t83;
t43 = t243 / 0.2e1 + (m(6) + t191) * t153;
t41 = t153 * t170 + t185 * t246;
t35 = (t190 + t191) * t155;
t29 = t160 + t165;
t20 = -t155 * t170 + t168;
t17 = t155 * t260 - t169 - t261;
t14 = m(7) * t250 + mrSges(7,2) * t213 / 0.2e1 + t212 * t254 + t176 * t255 + (m(6) * t150 - mrSges(6,1)) * t153 - t167;
t12 = t159 - t260;
t9 = -t184 * t247 + t162 / 0.2e1 + t161 + t224 * t266;
t5 = t6 * qJD(4);
t3 = t155 * (-t238 - t240) / 0.4e1 + t185 * t250 - t233 / 0.4e1 - t228 / 0.4e1 - t152 * (-Ifges(7,2) * t213 + t131) / 0.4e1 + t236 / 0.2e1 - t207 * t244 - t208 * t249 + Ifges(7,3) * t247 + t171 * t155 + t173 + (t183 / 0.2e1 - t180 / 0.4e1) * t212 + t261 * t252 - (t181 + t182) * t213 / 0.4e1;
t10 = [qJD(2) * t7 + qJD(3) * t13 + qJD(4) * t1 + qJD(5) * t15 + qJD(6) * t4, qJD(3) * t29 + qJD(4) * t9 + qJD(5) * t35 + qJD(6) * t20 + t222, qJD(2) * t29 + qJD(6) * t12 + t218 + t5, t9 * qJD(2) + t14 * qJD(5) + t3 * qJD(6) + (pkin(4) * mrSges(6,1) + Ifges(7,5) * t245 + Ifges(7,6) * t248 + Ifges(6,4) - Ifges(5,5)) * t199 + t179 + (t110 * t184 + t66 * t244 + t64 * t249 + (t181 * t244 + t183 * t248 + Ifges(6,5) - Ifges(5,6)) * t155 + (m(6) * t202 - t116 + t186) * t150 + (m(7) * t110 - t155 * mrSges(6,1) - t87) * qJ(5) - (m(7) * t176 + t210 + t214) * t252 - t176 * mrSges(7,3)) * qJD(4), qJD(2) * t35 + qJD(4) * t14 + qJD(6) * t17 + t217, t223 + t20 * qJD(2) + t12 * qJD(3) + t3 * qJD(4) + t17 * qJD(5) + (-mrSges(7,1) * t37 - mrSges(7,2) * t36 - Ifges(7,6) * t213 + t130) * qJD(6); qJD(3) * t28 + qJD(4) * t8 + qJD(5) * t32 - qJD(6) * t19 - t222, 0, t204, t221, t203, qJD(6) * t185 - t205; -qJD(2) * t28 + qJD(6) * t11 - t218 + t5, -t204, qJD(4) * t251 (t155 * t184 - t116 - t226) * qJD(4) + t43 * qJD(5) + t41 * qJD(6) + (-mrSges(5,1) + t268) * t199 + 0.2e1 * ((t153 * t262 + t136) * t255 + t202 * t256) * qJD(4) + t172, t43 * qJD(4), t41 * qJD(4) + t155 * t267 + t219; -qJD(2) * t8 - qJD(5) * t18 - qJD(6) * t2 - t179, -t221, -qJD(5) * t48 - t172, qJD(5) * t83 - qJD(6) * t27, -t163 ((mrSges(7,2) * t252 - Ifges(7,6)) * t154 + (mrSges(7,1) * t252 - Ifges(7,5)) * t152) * qJD(6) - t175; -qJD(2) * t32 + qJD(4) * t18 - qJD(6) * t16 - t217, -t203, t48 * qJD(4), t163, 0, -t206 - t267; qJD(2) * t19 - qJD(3) * t11 + qJD(4) * t2 + qJD(5) * t16 - t223, t205, -t219, t175, t206, 0;];
Cq  = t10;
