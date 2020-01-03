% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP4
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:56
% EndTime: 2019-12-31 20:55:01
% DurationCPUTime: 2.46s
% Computational Cost: add. (6599->203), mult. (13038->259), div. (0->0), fcn. (13989->6), ass. (0->132)
t286 = qJD(2) + qJD(3);
t245 = m(6) / 0.2e1;
t285 = 0.2e1 * t245;
t157 = sin(qJ(3));
t158 = sin(qJ(2));
t159 = cos(qJ(3));
t160 = cos(qJ(2));
t138 = -t157 * t158 + t159 * t160;
t139 = -t157 * t160 - t159 * t158;
t156 = sin(pkin(8));
t209 = cos(pkin(8));
t103 = t209 * t138 + t139 * t156;
t174 = t156 * t138 - t209 * t139;
t151 = -pkin(2) * t160 - pkin(1);
t117 = -t138 * pkin(3) + t151;
t47 = -pkin(4) * t103 - qJ(5) * t174 + t117;
t282 = -m(6) * t47 + mrSges(6,1) * t103 + mrSges(6,3) * t174;
t242 = t103 / 0.2e1;
t237 = t160 * pkin(6);
t144 = t160 * pkin(7) + t237;
t196 = (-pkin(7) - pkin(6)) * t158;
t114 = t159 * t144 + t157 * t196;
t169 = t138 * qJ(4) + t114;
t252 = -t157 * t144 + t159 * t196;
t256 = t139 * qJ(4) + t252;
t262 = t156 * t256 + t209 * t169;
t221 = t103 * mrSges(6,2);
t94 = t221 / 0.2e1;
t281 = mrSges(6,2) * t242 + t262 * t285 + t94;
t284 = t281 * qJD(5);
t283 = m(5) * t117;
t260 = -t103 / 0.2e1;
t267 = (t242 + t260) * mrSges(6,2);
t280 = qJD(1) * t267;
t279 = qJD(5) * t267;
t177 = Ifges(4,5) * t138 + Ifges(4,6) * t139 + (-Ifges(5,6) + Ifges(6,6)) * t174 + (Ifges(6,4) + Ifges(5,5)) * t103;
t258 = t252 * mrSges(4,2);
t259 = t114 * mrSges(4,1);
t272 = t262 * mrSges(6,1);
t273 = t262 * mrSges(5,1);
t261 = -t156 * t169 + t209 * t256;
t274 = t261 * mrSges(6,3);
t275 = t261 * mrSges(5,2);
t278 = t177 - t258 - t259 + t274 - t272 - t273 - t275;
t277 = -t258 / 0.2e1 - t259 / 0.2e1 - t272 / 0.2e1 - t273 / 0.2e1 + t274 / 0.2e1 - t275 / 0.2e1;
t145 = t209 * t157 * pkin(2);
t240 = pkin(2) * t159;
t150 = pkin(3) + t240;
t125 = t156 * t150 + t145;
t122 = qJ(5) + t125;
t205 = t156 * t157;
t124 = -pkin(2) * t205 + t209 * t150;
t123 = -pkin(4) - t124;
t271 = t122 * t261 + t123 * t262;
t270 = -t124 * t262 + t125 * t261;
t238 = t156 * pkin(3);
t146 = qJ(5) + t238;
t195 = t209 * pkin(3);
t147 = -t195 - pkin(4);
t269 = t146 * t261 + t147 * t262;
t268 = t156 * t261 - t209 * t262;
t206 = t103 ^ 2;
t207 = t174 ^ 2;
t63 = mrSges(6,1) * t174 - t103 * mrSges(6,3);
t64 = mrSges(5,1) * t174 + t103 * mrSges(5,2);
t266 = t151 * (-mrSges(4,1) * t139 + mrSges(4,2) * t138) + t117 * t64 + t47 * t63 + (-Ifges(5,4) + Ifges(6,5)) * (-t206 + t207) + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t174 * t103;
t228 = mrSges(6,2) + mrSges(5,3);
t131 = t156 * t240 + t145;
t132 = (t209 * t159 - t205) * pkin(2);
t255 = (t157 * mrSges(4,1) + t159 * mrSges(4,2)) * pkin(2) + (mrSges(5,2) - mrSges(6,3)) * t132 + (mrSges(5,1) + mrSges(6,1)) * t131;
t248 = -m(5) / 0.2e1;
t247 = m(5) / 0.2e1;
t246 = -m(6) / 0.2e1;
t244 = m(5) * pkin(3);
t241 = -t174 / 0.2e1;
t239 = t139 * pkin(3);
t153 = t158 * pkin(2);
t227 = Ifges(4,4) * t139;
t5 = (m(5) + m(6)) * (t103 * t262 - t174 * t261) + t228 * (t206 + t207);
t226 = qJD(1) * t5;
t121 = t153 - t239;
t185 = Ifges(4,4) * t138 + (-Ifges(4,1) + Ifges(4,2)) * t139;
t57 = pkin(4) * t174 - qJ(5) * t103 - t239;
t49 = t153 + t57;
t66 = -mrSges(5,1) * t103 + mrSges(5,2) * t174;
t1 = m(4) * t151 * t153 + (-mrSges(4,2) * t153 - t227) * t139 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t158) * t158 + (-mrSges(4,1) * t153 + t185) * t138 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t160 + (Ifges(3,1) - Ifges(3,2)) * t158) * t160 + t266 - t282 * t49 + (t283 + t66) * t121;
t225 = t1 * qJD(1);
t224 = t174 * mrSges(6,2);
t223 = t174 * mrSges(5,3);
t220 = t103 * mrSges(5,3);
t2 = t185 * t138 + (-pkin(3) * t66 - t227) * t139 - t239 * t283 + t266 - t282 * t57;
t212 = t2 * qJD(1);
t165 = (t103 * t125 - t124 * t174) * t247 + (t103 * t122 + t123 * t174) * t245;
t172 = t63 + t64;
t178 = t121 * t247 + t49 * t245;
t6 = -t165 + t172 + t178;
t211 = t6 * qJD(1);
t203 = t244 / 0.2e1;
t164 = (t103 * t146 + t147 * t174) * t245 + (t103 * t156 - t174 * t209) * t203;
t175 = t239 * t248 + t57 * t245;
t9 = -t164 + t172 + t175;
t210 = t9 * qJD(1);
t16 = t282 * t174;
t208 = qJD(1) * t16;
t59 = t174 * m(6);
t204 = t59 * qJD(1);
t202 = t122 * t224;
t201 = t123 * t221;
t200 = t124 * t220;
t199 = t125 * t223;
t198 = t146 * t224;
t197 = t131 * t245;
t187 = t223 * t238;
t184 = -t131 * t261 + t132 * t262;
t183 = t195 * t220;
t17 = -m(5) * (-t124 * t131 + t125 * t132) - m(6) * (t122 * t132 + t123 * t131) + t255;
t161 = (t184 + t270) * t248 + (t184 + t271) * t246 + t202 / 0.2e1 - t201 / 0.2e1 + t200 / 0.2e1 + t199 / 0.2e1 + t228 * (t131 * t241 + t132 * t260) - t277;
t162 = t269 * t245 + t268 * t203 - t198 / 0.2e1 + t147 * t94 - t187 / 0.2e1 - t183 / 0.2e1 + t277;
t4 = t161 + t162;
t182 = -t4 * qJD(1) - t17 * qJD(2);
t115 = m(6) * t122 + mrSges(6,3);
t179 = -qJD(2) * t115 - t280;
t140 = m(6) * t146 + mrSges(6,3);
t168 = -mrSges(6,3) + (0.2e1 * qJ(5) + t125 + t238) * t246;
t76 = t197 + t168;
t170 = -qJD(2) * t76 + qJD(3) * t140 + t280;
t77 = t197 - t168;
t60 = m(6) * t241 + t174 * t245;
t11 = t164 + t175;
t10 = t165 + t178;
t3 = -t161 + t162 + t177;
t7 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t5 + qJD(5) * t16, t3 * qJD(3) + t10 * qJD(4) + t284 + t225 + (-mrSges(3,1) * t237 + Ifges(3,5) * t160 - t199 - t200 + t201 - t202 + t271 * t285 + 0.2e1 * t270 * t247 + (m(4) * (-t114 * t159 + t157 * t252) + (-t138 * t159 + t139 * t157) * mrSges(4,3)) * pkin(2) + (mrSges(3,2) * pkin(6) - Ifges(3,6)) * t158 + t278) * qJD(2), t212 + t3 * qJD(2) + (m(6) * t269 + t147 * t221 + t268 * t244 - t183 - t187 - t198 + t278) * qJD(3) + t11 * qJD(4) + t284, qJD(2) * t10 + qJD(3) * t11 + qJD(5) * t60 + t226, qJD(4) * t60 + t286 * t281 + t208; -qJD(3) * t4 - qJD(4) * t6 - t225 + t279, -qJD(3) * t17 + qJD(5) * t115, ((-t209 * t131 + t132 * t156) * t244 + m(6) * (t131 * t147 + t132 * t146) - t255) * qJD(3) + t77 * qJD(5) + t182, -t211, qJD(3) * t77 - t179; qJD(2) * t4 - qJD(4) * t9 - t212 + t279, -qJD(5) * t76 - t182, t140 * qJD(5), -t210, t170; qJD(2) * t6 + qJD(3) * t9 - qJD(5) * t59 - t226, t211, t210, 0, -t204; qJD(4) * t59 - t286 * t267 - t208, qJD(3) * t76 + t179, -t170, t204, 0;];
Cq = t7;
