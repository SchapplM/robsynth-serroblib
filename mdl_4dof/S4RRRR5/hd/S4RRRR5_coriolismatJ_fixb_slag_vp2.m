% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR5
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:24
% EndTime: 2019-12-31 17:27:30
% DurationCPUTime: 2.81s
% Computational Cost: add. (5557->317), mult. (12216->461), div. (0->0), fcn. (11468->6), ass. (0->179)
t183 = sin(qJ(3));
t305 = -pkin(7) - pkin(6);
t168 = t305 * t183;
t186 = cos(qJ(3));
t169 = t305 * t186;
t182 = sin(qJ(4));
t185 = cos(qJ(4));
t102 = t168 * t182 - t169 * t185;
t187 = cos(qJ(2));
t154 = -t182 * t186 - t185 * t183;
t184 = sin(qJ(2));
t131 = t154 * t184;
t259 = t131 * mrSges(5,3);
t103 = mrSges(5,2) * t187 + t259;
t231 = t186 * t184;
t235 = t183 * t184;
t130 = t182 * t235 - t185 * t231;
t153 = -t182 * t183 + t185 * t186;
t289 = pkin(3) * t183;
t223 = pkin(5) + t289;
t162 = t223 * t184;
t288 = pkin(3) * t186;
t178 = -pkin(2) - t288;
t272 = Ifges(5,4) * t154;
t211 = Ifges(5,1) * t153 + t272;
t215 = -t154 * mrSges(5,1) + t153 * mrSges(5,2);
t105 = -mrSges(5,1) * t187 + t130 * mrSges(5,3);
t303 = -t105 / 0.2e1;
t145 = Ifges(5,4) * t153;
t89 = -Ifges(5,1) * t154 + t145;
t325 = Ifges(5,2) * t154 + t145 + t89;
t221 = t185 * t168 + t169 * t182;
t326 = -t221 / 0.2e1;
t273 = Ifges(5,4) * t130;
t61 = Ifges(5,2) * t131 - Ifges(5,6) * t187 - t273;
t121 = Ifges(5,4) * t131;
t63 = -Ifges(5,1) * t130 - Ifges(5,5) * t187 + t121;
t70 = -mrSges(5,1) * t130 + mrSges(5,2) * t131;
t73 = Ifges(5,2) * t130 + t121;
t74 = Ifges(5,1) * t131 + t273;
t87 = Ifges(5,5) * t153 + Ifges(5,6) * t154;
t88 = Ifges(5,2) * t153 - t272;
t327 = -t102 * t303 + t103 * t326 - t162 * t215 / 0.2e1 + t187 * t87 / 0.4e1 - t178 * t70 / 0.2e1 - t325 * t131 / 0.4e1 - (t73 + t63) * t153 / 0.4e1 - (-t74 / 0.4e1 + t61 / 0.4e1) * t154 - (-t211 / 0.4e1 + t88 / 0.4e1) * t130;
t294 = t184 / 0.2e1;
t297 = -t154 / 0.2e1;
t11 = -t178 * t215 + t88 * t297 + t154 * t211 / 0.2e1 - t325 * t153 / 0.2e1;
t14 = -t102 * mrSges(5,1) - t221 * mrSges(5,2) + t87;
t324 = t14 * qJD(4);
t254 = t154 * mrSges(5,3);
t224 = t254 / 0.2e1;
t283 = t184 * pkin(6);
t165 = -pkin(2) * t187 - pkin(1) - t283;
t230 = t186 * t187;
t116 = pkin(5) * t230 + t183 * t165;
t95 = -pkin(7) * t235 + t116;
t247 = t185 * t95;
t151 = t186 * t165;
t219 = -pkin(7) * t231 + t151;
t81 = (-pkin(5) * t183 - pkin(3)) * t187 + t219;
t31 = t182 * t81 + t247;
t318 = t131 * t326 + t102 * t130 / 0.2e1;
t323 = t318 * mrSges(5,3) + t31 * t224 - t327;
t257 = t153 * mrSges(5,3);
t132 = t154 * t187;
t133 = t153 * t187;
t202 = -Ifges(5,5) * t133 / 0.2e1 - Ifges(5,6) * t132 / 0.2e1;
t285 = t184 * pkin(2);
t170 = -pkin(6) * t187 + t285;
t117 = pkin(5) * t235 + t186 * t170;
t82 = t184 * pkin(3) - pkin(7) * t230 + t117;
t118 = -pkin(5) * t231 + t183 * t170;
t234 = t183 * t187;
t96 = -pkin(7) * t234 + t118;
t34 = -t182 * t96 + t185 * t82;
t35 = t182 * t82 + t185 * t96;
t218 = Ifges(5,3) * t294 - t35 * mrSges(5,2) / 0.2e1 + t34 * mrSges(5,1) / 0.2e1 - t202;
t268 = Ifges(4,6) * t183;
t271 = Ifges(4,5) * t186;
t203 = t271 / 0.2e1 - t268 / 0.2e1;
t316 = Ifges(3,4) - t203;
t314 = -t117 * t183 + t118 * t186;
t217 = t186 * mrSges(4,1) - t183 * mrSges(4,2);
t250 = t182 * t95;
t30 = t185 * t81 - t250;
t228 = pkin(5) * t234;
t94 = t219 - t228;
t40 = t185 * t94 - t250;
t225 = -t40 / 0.2e1 + t30 / 0.2e1;
t311 = t183 ^ 2;
t310 = -m(5) / 0.2e1;
t309 = m(5) / 0.2e1;
t308 = mrSges(4,1) / 0.2e1;
t307 = -mrSges(4,2) / 0.2e1;
t306 = -t31 / 0.2e1;
t304 = pkin(2) * mrSges(4,2);
t302 = -t130 / 0.2e1;
t301 = t131 / 0.2e1;
t300 = t132 / 0.2e1;
t299 = t133 / 0.2e1;
t296 = -t183 / 0.2e1;
t295 = t183 / 0.2e1;
t293 = t185 / 0.2e1;
t292 = -t186 / 0.2e1;
t291 = t186 / 0.2e1;
t290 = m(5) * t162;
t158 = mrSges(4,2) * t187 - mrSges(4,3) * t235;
t287 = pkin(6) * t158;
t160 = -mrSges(4,1) * t187 - mrSges(4,3) * t231;
t286 = pkin(6) * t160;
t284 = t184 * pkin(5);
t282 = t30 * mrSges(5,2);
t281 = t31 * mrSges(5,1);
t39 = -t182 * t94 - t247;
t278 = t39 * mrSges(5,1);
t277 = t40 * mrSges(5,2);
t275 = Ifges(4,4) * t183;
t274 = Ifges(4,4) * t186;
t270 = Ifges(4,5) * t187;
t267 = Ifges(4,6) * t187;
t264 = pkin(3) * qJD(3);
t249 = t183 * mrSges(4,1);
t104 = -mrSges(5,2) * t184 + mrSges(5,3) * t132;
t106 = mrSges(5,1) * t184 - mrSges(5,3) * t133;
t115 = t151 - t228;
t210 = -Ifges(4,2) * t183 + t274;
t127 = Ifges(4,6) * t184 + t210 * t187;
t213 = Ifges(4,1) * t186 - t275;
t129 = Ifges(4,5) * t184 + t213 * t187;
t216 = t186 * mrSges(4,2) + t249;
t142 = t187 * t216;
t159 = -t184 * mrSges(4,2) - mrSges(4,3) * t234;
t161 = t184 * mrSges(4,1) - mrSges(4,3) * t230;
t163 = t223 * t187;
t197 = t213 * t184;
t128 = t197 - t270;
t232 = t186 * t128;
t196 = t210 * t184;
t126 = t196 - t267;
t236 = t183 * t126;
t62 = Ifges(5,4) * t133 + Ifges(5,2) * t132 + Ifges(5,6) * t184;
t64 = Ifges(5,1) * t133 + Ifges(5,4) * t132 + Ifges(5,5) * t184;
t71 = -mrSges(5,1) * t131 - mrSges(5,2) * t130;
t72 = -mrSges(5,1) * t132 + mrSges(5,2) * t133;
t3 = t118 * t158 + t116 * t159 + t117 * t160 + t115 * t161 + t162 * t72 + t163 * t71 + t64 * t302 + t62 * t301 + t61 * t300 + t63 * t299 + t35 * t103 + t31 * t104 + t34 * t105 + t30 * t106 + m(5) * (t162 * t163 + t30 * t34 + t31 * t35) + m(4) * (t115 * t117 + t116 * t118) + (-pkin(1) * mrSges(3,2) - t236 / 0.2e1 + t232 / 0.2e1 + t202 + t316 * t187) * t187 + (-pkin(1) * mrSges(3,1) + Ifges(5,5) * t302 + Ifges(5,6) * t301 + pkin(5) * t142 + t127 * t296 + t129 * t291 - t316 * t184 + (-Ifges(4,3) + Ifges(3,1) - Ifges(3,2) - Ifges(5,3) + (m(4) * pkin(5) + t216) * pkin(5)) * t187) * t184;
t244 = t3 * qJD(1);
t229 = Ifges(5,5) * t131 + Ifges(5,6) * t130;
t191 = t162 * t70 + (t63 / 0.2e1 + t73 / 0.2e1) * t131 + (t31 * mrSges(5,3) - t74 / 0.2e1 + t61 / 0.2e1) * t130 - t30 * t259 - t187 * t229 / 0.2e1;
t207 = Ifges(4,5) * t183 + Ifges(4,6) * t186;
t209 = Ifges(4,2) * t186 + t275;
t212 = Ifges(4,1) * t183 + t274;
t6 = t40 * t103 + t39 * t105 + m(5) * (t30 * t39 + t31 * t40) + t115 * t158 - t116 * t160 + (t128 * t296 + t126 * t292 + t187 * t207 / 0.2e1 + (pkin(5) * t217 + t209 * t295 + t212 * t292) * t184 + (t71 + t290) * t288 + (t115 * t183 - t116 * t186) * mrSges(4,3)) * t184 + t191;
t243 = t6 * qJD(1);
t7 = t30 * t103 - t31 * t105 + t191;
t242 = t7 * qJD(1);
t237 = t182 * t130;
t233 = t185 * t131;
t227 = pkin(6) * mrSges(4,3) / 0.2e1;
t226 = t306 - t39 / 0.2e1;
t220 = -mrSges(4,3) * t283 / 0.2e1;
t214 = -t153 * mrSges(5,1) - t154 * mrSges(5,2);
t208 = -t268 + t271;
t188 = (t182 * t104 / 0.2e1 + t106 * t293 + (t182 * t35 + t185 * t34) * t309) * pkin(3) + Ifges(4,3) * t294 + t117 * t308 + t118 * t307 + t218;
t194 = (t31 + t39) * t221 + (t40 - t30) * t102;
t201 = pkin(3) * t214;
t1 = (-t128 / 0.4e1 + t286 / 0.2e1 + 0.3e1 / 0.4e1 * t270 + (pkin(3) * t178 * t310 - t201 / 0.2e1 + pkin(2) * t308 + pkin(5) * t307 + (Ifges(4,2) / 0.2e1 + t227 - Ifges(4,1) / 0.4e1) * t186) * t184) * t186 + (-0.3e1 / 0.4e1 * t267 + t126 / 0.4e1 + t287 / 0.2e1 + (-t71 / 0.2e1 - t290 / 0.2e1) * pkin(3) + (-t304 / 0.2e1 + 0.3e1 / 0.2e1 * t274 - pkin(5) * mrSges(4,1) / 0.2e1 + (Ifges(4,1) / 0.2e1 + t227 - Ifges(4,2) / 0.4e1) * t183) * t184) * t183 + (t225 * t153 + t226 * t154 - t318) * mrSges(5,3) + t194 * t310 + t188 + t327;
t8 = -t183 * t201 - m(5) * t178 * t289 + pkin(2) * t249 + Ifges(4,4) * t311 + (t304 - t274 + (-Ifges(4,1) + Ifges(4,2)) * t183) * t186 + t11;
t206 = -t1 * qJD(1) - t8 * qJD(2);
t189 = t254 * t306 + t323;
t5 = t189 - t218;
t204 = t5 * qJD(1) - t11 * qJD(2);
t192 = (t103 * t293 + t182 * t303 + (t237 / 0.2e1 - t233 / 0.2e1) * mrSges(5,3)) * pkin(3);
t10 = t226 * mrSges(5,1) - t225 * mrSges(5,2) + t192;
t164 = (mrSges(5,1) * t182 + mrSges(5,2) * t185) * pkin(3);
t195 = -qJD(1) * t10 + qJD(3) * t164;
t152 = t164 * qJD(4);
t9 = -t282 / 0.2e1 - t281 / 0.2e1 - t277 / 0.2e1 + t278 / 0.2e1 + t192 + t229;
t4 = t189 + t218;
t2 = -t183 * t196 / 0.4e1 + t232 / 0.4e1 - t236 / 0.4e1 + t311 * t220 + t188 - t209 * t231 / 0.2e1 + t201 * t231 / 0.2e1 + ((t162 * t183 + t178 * t231) * pkin(3) + t194) * t309 - t212 * t235 / 0.2e1 + t39 * t224 + t287 * t296 + t71 * t289 / 0.2e1 + t216 * t284 / 0.2e1 + t286 * t292 - t217 * t285 / 0.2e1 - t225 * t257 + (-t208 / 0.4e1 + t203) * t187 + (t197 / 0.4e1 + t220 * t186) * t186 + t323;
t12 = [qJD(2) * t3 + qJD(3) * t6 + qJD(4) * t7, t2 * qJD(3) + t4 * qJD(4) + t244 + (t127 * t291 - Ifges(3,6) * t184 + t129 * t295 + t178 * t72 + t64 * t297 + t163 * t214 + t153 * t62 / 0.2e1 - pkin(2) * t142 + t88 * t300 + t89 * t299 + t102 * t104 + t221 * t106 + t35 * t257 + t34 * t254 + m(5) * (t102 * t35 + t163 * t178 + t221 * t34) + mrSges(3,2) * t284 + (m(4) * t314 + t186 * t159 - t183 * t161) * pkin(6) + (Ifges(3,5) + t212 * t291 + t209 * t296 + (-m(4) * pkin(2) - mrSges(3,1) - t217) * pkin(5)) * t187 + (-Ifges(5,5) * t154 + Ifges(5,6) * t153 + t207) * t294 + t314 * mrSges(4,3)) * qJD(2), t243 + t2 * qJD(2) + (-t116 * mrSges(4,1) - t115 * mrSges(4,2) - Ifges(4,5) * t235 - Ifges(4,6) * t231 + t229 - t277 + t278) * qJD(3) + t9 * qJD(4) + (m(5) * (t182 * t40 + t185 * t39) + (-t233 + t237) * mrSges(5,3)) * t264, t242 + t4 * qJD(2) + t9 * qJD(3) + (t229 - t281 - t282) * qJD(4); -qJD(3) * t1 + qJD(4) * t5 - t244, -qJD(3) * t8 - qJD(4) * t11, (-t217 * pkin(6) + t14 + t208) * qJD(3) + t324 + (m(5) * (-t102 * t185 + t182 * t221) + (-t153 * t185 + t154 * t182) * mrSges(5,3)) * t264 + t206, t14 * qJD(3) + t204 + t324; qJD(2) * t1 + qJD(4) * t10 - t243, -t206, -t152, -t152 - t195; -qJD(2) * t5 - qJD(3) * t10 - t242, -t204, t195, 0;];
Cq = t12;
