% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:27
% EndTime: 2019-12-31 20:49:31
% DurationCPUTime: 2.13s
% Computational Cost: add. (4527->205), mult. (9247->237), div. (0->0), fcn. (8580->6), ass. (0->129)
t297 = Ifges(5,4) - Ifges(6,5);
t164 = cos(qJ(3));
t296 = t164 ^ 2;
t162 = sin(qJ(3));
t226 = cos(pkin(8));
t194 = t226 * t162;
t225 = sin(pkin(8));
t175 = t164 * t225 + t194;
t192 = t297 * t175;
t193 = t225 * t162;
t131 = -t164 * t226 + t193;
t295 = t297 * t131;
t280 = -Ifges(5,1) - Ifges(6,1);
t273 = Ifges(6,3) + t280;
t289 = -Ifges(4,4) * t162 + (Ifges(4,1) - Ifges(4,2)) * t164;
t293 = t296 * Ifges(4,4);
t93 = mrSges(5,1) * t131 + mrSges(5,2) * t175;
t294 = (pkin(3) * t93 + t289) * t162 + t293 + (-(-Ifges(5,2) - t273) * t175 + t295) * t131;
t281 = mrSges(6,2) + mrSges(5,3);
t251 = pkin(3) * t162;
t87 = pkin(4) * t175 + qJ(5) * t131 + t251;
t92 = mrSges(6,1) * t131 - mrSges(6,3) * t175;
t247 = t87 * t92;
t290 = -t192 * t175 + t247 + t294;
t283 = m(5) + m(6);
t217 = t162 ^ 2 + t296;
t237 = mrSges(4,2) * t162;
t141 = -mrSges(4,1) * t164 + t237;
t285 = -mrSges(3,1) + t141 + t92 + t93;
t282 = mrSges(5,1) + mrSges(6,1);
t163 = sin(qJ(2));
t254 = pkin(1) * t163;
t203 = pkin(7) + t254;
t201 = t225 * pkin(3);
t146 = t201 + qJ(5);
t202 = t226 * pkin(3);
t153 = -t202 - pkin(4);
t262 = m(5) * pkin(3);
t214 = t262 / 0.2e1;
t263 = m(6) / 0.2e1;
t169 = (-t131 * t146 + t153 * t175) * t263 + (-t131 * t225 - t175 * t226) * t214;
t180 = t162 * t214 + t263 * t87;
t90 = mrSges(6,1) * t175 + t131 * mrSges(6,3);
t91 = mrSges(5,1) * t175 - t131 * mrSges(5,2);
t14 = -t169 + t180 + t90 + t91;
t216 = qJD(1) + qJD(2);
t279 = t216 * t14;
t86 = m(6) * t175;
t278 = t216 * t86;
t142 = mrSges(4,1) * t162 + mrSges(4,2) * t164;
t264 = m(5) / 0.2e1;
t275 = t264 + t263;
t137 = m(6) * t146 + mrSges(6,3);
t186 = t162 * (-qJ(4) - t203);
t145 = t164 * t203;
t156 = t164 * qJ(4);
t219 = t145 + t156;
t81 = -t226 * t186 + t219 * t225;
t218 = t164 * pkin(7) + t156;
t246 = -qJ(4) - pkin(7);
t95 = -t246 * t194 + t218 * t225;
t272 = t225 * t186 + t219 * t226;
t271 = t246 * t193 + t218 * t226;
t270 = m(6) * t153 - t226 * t262 - t282;
t269 = t225 * t262 - mrSges(5,2) + t137;
t261 = t272 / 0.2e1;
t260 = t271 / 0.2e1;
t165 = cos(qJ(2));
t253 = pkin(1) * t165;
t252 = pkin(2) * t142;
t155 = -t164 * pkin(3) - pkin(2);
t84 = t131 * pkin(4) - qJ(5) * t175 + t155;
t73 = t84 - t253;
t249 = t73 * t90;
t248 = t84 * t90;
t245 = t14 * qJD(3) - t86 * qJD(5);
t21 = t169 + t180;
t244 = t21 * qJD(3);
t178 = t281 * (t131 ^ 2 + t175 ^ 2);
t48 = t272 * t131;
t9 = t178 + t283 * (t175 * t81 - t48);
t232 = qJD(1) * t9;
t231 = t131 * mrSges(6,2);
t139 = t155 - t253;
t230 = t139 * t91;
t229 = t155 * t91;
t128 = t139 * t251;
t154 = -pkin(2) - t253;
t221 = t154 * t142;
t29 = t73 * t87;
t3 = m(5) * t128 + m(6) * t29 + t221 + t230 + t249 + t290;
t228 = t3 * qJD(1);
t104 = t175 * t253;
t105 = t131 * t253;
t172 = t281 * (t104 * t175 + t105 * t131);
t5 = t172 + t283 * (t104 * t81 - t105 * t272) + (m(4) * t154 + m(5) * t139 + m(6) * t73 + t285) * t254 + (-mrSges(3,2) + (m(4) * t203 + mrSges(4,3)) * t217) * t253;
t227 = t5 * qJD(1);
t59 = t175 * t92;
t25 = -t73 * t86 - t59;
t224 = qJD(1) * t25;
t213 = t104 * t263;
t190 = t275 * t254;
t166 = -t142 * t253 / 0.2e1 - (t146 * t263 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1 + t225 * t214) * t105 + (-t282 / 0.2e1 + t153 * t263 - t214 * t226) * t104;
t140 = t155 * t251;
t32 = t84 * t87;
t167 = t247 + (t128 + t140) * t264 + (t29 + t32) * t263;
t2 = t221 / 0.2e1 - t252 / 0.2e1 + t249 / 0.2e1 + t248 / 0.2e1 + t229 / 0.2e1 + t230 / 0.2e1 + t167 - t166 + t93 * t251 + t293 + t289 * t162 + t295 * t131 + ((Ifges(6,3) / 0.2e1 + t273 / 0.2e1 + Ifges(5,2) + t280 / 0.2e1) * t131 - t192) * t175;
t4 = m(5) * t140 + m(6) * t32 + t229 + t248 - t252 + t290;
t189 = t2 * qJD(1) + t4 * qJD(2);
t72 = t271 * t131;
t11 = t178 + t283 * (t175 * t95 - t72);
t168 = -t178 + (-m(5) / 0.2e1 - m(6) / 0.2e1) * (-t48 - t72 - (-t81 - t95) * t175);
t6 = t190 + t168;
t188 = qJD(1) * t6 - qJD(2) * t11;
t181 = t59 + (t73 + t84) * t86 / 0.2e1;
t17 = t213 + t181;
t26 = -t84 * t86 - t59;
t187 = -qJD(1) * t17 + qJD(2) * t26;
t177 = qJD(3) * t137;
t170 = Ifges(4,5) * t164 - Ifges(4,6) * t162 - t153 * t231 + (mrSges(5,3) * t202 - Ifges(6,4) - Ifges(5,5)) * t131 + (-mrSges(6,2) * t146 - mrSges(5,3) * t201 - Ifges(5,6) + Ifges(6,6)) * t175;
t79 = t86 * qJD(4);
t27 = m(6) * t260 + t263 * t271 - t231;
t22 = m(6) * t261 + t263 * t272 - t231;
t19 = t21 * qJD(4);
t18 = t213 - t181;
t13 = t14 * qJD(4);
t7 = t190 - t168;
t1 = (t84 / 0.2e1 + t73 / 0.2e1) * t90 + (-pkin(2) / 0.2e1 + t154 / 0.2e1) * t142 + t166 - (t192 + t281 * (t260 - t271 / 0.2e1 + t261 - t272 / 0.2e1)) * t175 + (t139 / 0.2e1 + t155 / 0.2e1) * t91 + t167 + t294;
t8 = [qJD(2) * t5 + qJD(3) * t3 + qJD(4) * t9 + qJD(5) * t25, t1 * qJD(3) + t7 * qJD(4) + t18 * qJD(5) + t227 + (t172 + 0.2e1 * t275 * (t104 * t95 - t105 * t271) + ((-m(4) * pkin(2) + m(5) * t155 + m(6) * t84 + t285) * t163 + (-mrSges(3,2) + (m(4) * pkin(7) + mrSges(4,3)) * t217) * t165) * pkin(1)) * qJD(2), t228 + t1 * qJD(2) + (-mrSges(4,1) * t145 + t203 * t237 - t269 * t81 + t270 * t272 + t170) * qJD(3) + t19 + t22 * qJD(5), qJD(2) * t7 + t232 + t244, qJD(2) * t18 + qJD(3) * t22 + t224; qJD(3) * t2 - qJD(4) * t6 - qJD(5) * t17 - t227, qJD(3) * t4 + qJD(4) * t11 + qJD(5) * t26, (t141 * pkin(7) - t269 * t95 + t270 * t271 + t170) * qJD(3) + t19 + t27 * qJD(5) + t189, -t188 + t244, qJD(3) * t27 + t187; -qJD(2) * t2 - t13 - t228, -t13 - t189, t137 * qJD(5), -t279, t177; qJD(2) * t6 - t232 + t245, t188 + t245, t279, 0, -t278; qJD(2) * t17 - t224 + t79, -t187 + t79, -t177, t278, 0;];
Cq = t8;
