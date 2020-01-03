% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:02
% EndTime: 2019-12-31 21:11:07
% DurationCPUTime: 1.91s
% Computational Cost: add. (4728->212), mult. (9136->256), div. (0->0), fcn. (7861->6), ass. (0->131)
t175 = sin(qJ(5));
t176 = sin(qJ(3));
t178 = cos(qJ(5));
t179 = cos(qJ(3));
t137 = -t179 * t175 + t176 * t178;
t209 = t176 * t175 + t179 * t178;
t213 = Ifges(6,5) * t209 + Ifges(6,6) * t137;
t268 = sin(qJ(2)) * pkin(1);
t161 = pkin(7) + t268;
t170 = t176 * pkin(8);
t132 = t161 * t176 - t170;
t133 = (-pkin(8) + t161) * t179;
t211 = t132 * t175 + t133 * t178;
t73 = t132 * t178 - t133 * t175;
t300 = t211 * mrSges(6,1) + t73 * mrSges(6,2);
t8 = t213 + t300;
t305 = t8 * qJD(5);
t153 = pkin(7) * t176 - t170;
t154 = (pkin(7) - pkin(8)) * t179;
t210 = t153 * t175 + t154 * t178;
t94 = t153 * t178 - t154 * t175;
t299 = t210 * mrSges(6,1) + t94 * mrSges(6,2);
t12 = t213 + t299;
t304 = t12 * qJD(5);
t275 = m(6) / 0.2e1;
t302 = t176 ^ 2 + t179 ^ 2;
t267 = cos(qJ(2)) * pkin(1);
t104 = t137 * t267;
t105 = t209 * t267;
t181 = -pkin(3) - pkin(4);
t141 = -t175 * qJ(4) + t178 * t181;
t142 = t178 * qJ(4) + t175 * t181;
t249 = t105 * mrSges(6,2) / 0.2e1 - t104 * mrSges(6,1) / 0.2e1;
t298 = t249 + (t104 * t141 + t105 * t142) * t275;
t295 = Ifges(4,4) - Ifges(5,5);
t162 = -pkin(2) - t267;
t163 = t176 * qJ(4);
t281 = t179 * pkin(3) + t163;
t117 = t162 - t281;
t171 = t179 * pkin(4);
t106 = -t117 + t171;
t79 = mrSges(6,1) * t137 - mrSges(6,2) * t209;
t294 = t106 * t79;
t234 = pkin(2) + t281;
t118 = t171 + t234;
t293 = t118 * t79;
t292 = (t106 / 0.2e1 + t118 / 0.2e1) * t79;
t291 = t302 * t267;
t217 = t175 * mrSges(6,1) + t178 * mrSges(6,2);
t290 = qJD(5) * t217;
t289 = -t175 * t73 + t178 * t211;
t288 = -t175 * t94 + t178 * t210;
t186 = -(Ifges(6,1) - Ifges(6,2)) * t209 * t137 + (-t137 ^ 2 + t209 ^ 2) * Ifges(6,4);
t252 = t209 * mrSges(6,3);
t144 = -t179 * mrSges(5,1) - t176 * mrSges(5,3);
t282 = -t179 * mrSges(4,1) + t176 * mrSges(4,2) + t144;
t191 = (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t179 - t295 * t176;
t222 = t295 * t179;
t280 = t191 * t176 + t222 * t179;
t80 = mrSges(6,1) * t209 + mrSges(6,2) * t137;
t278 = (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t302) * t267 + (-t104 * t137 - t105 * t209) * mrSges(6,3) + (-mrSges(3,1) - t80 + t282) * t268;
t247 = qJ(4) * t179;
t134 = t181 * t176 + t247;
t277 = t134 * t80 - t186;
t276 = m(5) / 0.2e1;
t271 = m(5) * qJ(4);
t146 = pkin(3) * t176 - t247;
t270 = m(5) * t146;
t269 = m(5) * t179;
t258 = mrSges(5,3) * t179;
t251 = t137 * mrSges(6,3);
t164 = t179 * mrSges(5,2);
t148 = mrSges(4,1) * t176 + mrSges(4,2) * t179;
t243 = t146 * t144;
t183 = t243 + t277 + t280;
t147 = mrSges(5,1) * t176 - t258;
t223 = t147 + t270;
t77 = t106 * t134;
t3 = m(6) * t77 + t223 * t117 + t162 * t148 + t183 - t294;
t250 = t3 * qJD(1);
t78 = t176 * t80;
t28 = t78 + (-m(5) * t117 + m(6) * t106 - t144) * t176;
t246 = qJD(1) * t28;
t33 = t142 * mrSges(6,1) + mrSges(6,2) * t141;
t245 = qJD(5) * t33;
t240 = t291 * t161;
t11 = m(6) * (t104 * t73 + t105 * t211 - t106 * t268) + m(5) * (t117 * t268 + t240) + m(4) * (t162 * t268 + t240) + t278;
t244 = t11 * qJD(1);
t242 = t117 - t234;
t239 = t291 * pkin(7);
t233 = -t267 / 0.2e1;
t232 = t267 / 0.2e1;
t225 = m(5) * t242;
t224 = m(6) * (t106 + t118);
t221 = m(5) * t232;
t219 = t179 * t233;
t218 = t242 * t276;
t85 = t118 * t134;
t182 = (-pkin(2) / 0.2e1 + t162 / 0.2e1) * t148 - t292 + (-t234 / 0.2e1 + t117 / 0.2e1) * t147 + (t77 + t85) * t275 + t277;
t1 = t219 * t271 + ((mrSges(4,2) / 0.2e1 - mrSges(5,3) / 0.2e1) * t267 + t222) * t179 + t146 * t218 + ((mrSges(4,1) / 0.2e1 + mrSges(5,1) / 0.2e1 + pkin(3) * t276) * t267 + t191) * t176 + t182 + t243 - t298;
t6 = m(6) * t85 - pkin(2) * t148 - t223 * t234 + t183 - t293;
t216 = t1 * qJD(1) + t6 * qJD(2);
t7 = t186 + t294;
t214 = t7 * qJD(1);
t197 = (t104 * t178 + t105 * t175) * t275;
t20 = t197 - t78 + (t221 + t144 - t224 / 0.2e1 + t225 / 0.2e1) * t176;
t32 = t78 + (m(5) * t234 + m(6) * t118 - t144) * t176;
t212 = -qJD(1) * t20 + qJD(2) * t32;
t205 = m(6) * t289;
t204 = m(6) * t288;
t200 = t175 * t251 - t178 * t252 + t164;
t199 = t289 * t275;
t198 = t288 * t275;
t195 = t33 * qJD(3);
t10 = t186 + t293;
t184 = t186 + t292;
t4 = t184 + t249;
t194 = t4 * qJD(1) + t10 * qJD(2);
t17 = t199 - t205 / 0.2e1;
t22 = t198 - t204 / 0.2e1;
t48 = mrSges(5,3) + t271 + m(6) * (-t141 * t175 + t142 * t178) + t217;
t193 = t17 * qJD(1) + t22 * qJD(2) + t48 * qJD(3);
t192 = t217 * qJD(3);
t187 = qJD(3) * (-m(5) * t281 + t282);
t185 = -mrSges(5,2) * t163 - pkin(3) * t164 - t141 * t252 + t142 * t251 - t213 + (Ifges(5,4) + Ifges(4,5)) * t179 + (-Ifges(4,6) + Ifges(5,6)) * t176;
t21 = t78 + t197 + (-t144 + t221) * t176 + (t224 - t225) * t176 / 0.2e1;
t19 = t198 + pkin(7) * t269 + t204 / 0.2e1 + t200;
t16 = t199 + t161 * t269 + t205 / 0.2e1 + t200;
t5 = t184 - t249;
t2 = (t218 + t144) * t146 + t182 + mrSges(4,2) * t219 + t232 * t258 + t280 + (t270 + (mrSges(4,1) + mrSges(5,1)) * t176) * t233 + t298;
t9 = [qJD(2) * t11 + qJD(3) * t3 + qJD(4) * t28 + qJD(5) * t7, t2 * qJD(3) + t21 * qJD(4) + t5 * qJD(5) + t244 + (0.2e1 * (t104 * t94 + t105 * t210 - t118 * t268) * t275 + 0.2e1 * (-t234 * t268 + t239) * t276 + m(4) * (-pkin(2) * t268 + t239) + t278) * qJD(2), t250 + t2 * qJD(2) + (m(6) * (t141 * t211 - t142 * t73) + t185 - t300) * qJD(3) + t16 * qJD(4) + t305 + t161 * t187, qJD(2) * t21 + qJD(3) * t16 + t246, t5 * qJD(2) + t8 * qJD(3) + t214 - t305; qJD(3) * t1 - qJD(4) * t20 + qJD(5) * t4 - t244, qJD(3) * t6 + qJD(4) * t32 + qJD(5) * t10, (m(6) * (t141 * t210 - t142 * t94) + t185 - t299) * qJD(3) + t19 * qJD(4) + t304 + pkin(7) * t187 + t216, qJD(3) * t19 + t212, t12 * qJD(3) + t194 - t304; -qJD(2) * t1 + qJD(4) * t17 - t250, qJD(4) * t22 - t216, qJD(4) * t48 + t245, t193, t195 - t245; qJD(2) * t20 - qJD(3) * t17 - t246, -qJD(3) * t22 - t212, -t193 + t290, 0, t192 - t290; -qJD(2) * t4 - t214, -t194, -qJD(4) * t217 - t195, -t192, 0;];
Cq = t9;
