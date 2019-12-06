% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:40
% EndTime: 2019-12-05 18:33:44
% DurationCPUTime: 2.53s
% Computational Cost: add. (9988->202), mult. (19245->263), div. (0->0), fcn. (20855->8), ass. (0->128)
t173 = sin(qJ(5));
t176 = cos(qJ(5));
t175 = sin(qJ(2));
t242 = t175 * pkin(1);
t164 = qJ(3) + t242;
t171 = sin(pkin(9));
t143 = (-pkin(7) - t164) * t171;
t172 = cos(pkin(9));
t168 = t172 * pkin(7);
t144 = t164 * t172 + t168;
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t109 = t177 * t143 - t144 * t174;
t154 = t171 * t177 + t172 * t174;
t150 = t154 * pkin(8);
t263 = t109 - t150;
t110 = t143 * t174 + t144 * t177;
t153 = -t171 * t174 + t172 * t177;
t243 = pkin(8) * t153;
t90 = t110 + t243;
t272 = -t173 * t90 + t176 * t263;
t53 = t173 * t263 + t176 * t90;
t295 = -t53 * mrSges(6,1) - t272 * mrSges(6,2);
t129 = t153 * t173 + t154 * t176;
t198 = t176 * t153 - t154 * t173;
t298 = Ifges(6,5) * t198 - Ifges(6,6) * t129;
t13 = t298 + t295;
t300 = t13 * qJD(5);
t158 = (-pkin(7) - qJ(3)) * t171;
t160 = qJ(3) * t172 + t168;
t133 = t158 * t174 + t160 * t177;
t104 = t133 + t243;
t132 = t177 * t158 - t160 * t174;
t262 = t132 - t150;
t273 = -t104 * t173 + t176 * t262;
t77 = t104 * t176 + t173 * t262;
t296 = -t77 * mrSges(6,1) - t273 * mrSges(6,2);
t16 = t298 + t296;
t299 = t16 * qJD(5);
t270 = t198 / 0.2e1;
t297 = 0.2e1 * t270;
t279 = t129 * mrSges(6,1);
t204 = t279 / 0.2e1;
t38 = mrSges(6,2) * t297 + 0.2e1 * t204;
t294 = t38 * qJD(3);
t130 = t154 * mrSges(5,1) + t153 * mrSges(5,2);
t211 = t198 * t173;
t212 = t129 * t176;
t256 = m(6) * pkin(4);
t283 = t198 * mrSges(6,2) + t279;
t25 = (t154 / 0.2e1 + t212 / 0.2e1 - t211 / 0.2e1) * t256 + t283 + t130;
t293 = t25 * qJD(3);
t244 = pkin(4) * t154;
t292 = m(6) * t244;
t165 = -t172 * pkin(3) - pkin(2);
t137 = -t153 * pkin(4) + t165;
t178 = cos(qJ(2));
t246 = pkin(1) * t178;
t134 = t137 - t246;
t289 = t134 * t283;
t288 = t137 * t283;
t287 = (t134 / 0.2e1 + t137 / 0.2e1) * t283;
t222 = t129 * mrSges(6,3);
t138 = t154 * t246;
t139 = t153 * t246;
t93 = -t138 * t176 - t139 * t173;
t94 = -t138 * t173 + t139 * t176;
t229 = t93 * mrSges(6,1) / 0.2e1 - t94 * mrSges(6,2) / 0.2e1;
t284 = t138 * mrSges(5,1) / 0.2e1 + t139 * mrSges(5,2) / 0.2e1 - t229;
t183 = Ifges(6,4) * t198 * t297 + (-Ifges(6,4) * t129 + (Ifges(6,1) - Ifges(6,2)) * t270 + (-Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t198) * t129;
t209 = t171 ^ 2 + t172 ^ 2;
t264 = t209 * mrSges(4,3);
t278 = -mrSges(3,2) + t264;
t224 = t198 * mrSges(6,3);
t208 = qJD(1) + qJD(2);
t266 = t208 * t25;
t265 = t208 * t38;
t200 = t209 * qJ(3);
t82 = -mrSges(6,1) * t198 + mrSges(6,2) * t129;
t260 = -mrSges(4,1) * t172 - mrSges(5,1) * t153 + mrSges(4,2) * t171 + mrSges(5,2) * t154 - mrSges(3,1) + t82;
t259 = t154 ^ 2;
t258 = m(5) / 0.2e1;
t257 = m(6) / 0.2e1;
t245 = pkin(2) * t175;
t231 = t25 * qJD(4) + t38 * qJD(5);
t39 = t204 - t279 / 0.2e1;
t67 = (t154 + t211 - t212) * t256 / 0.2e1;
t230 = t67 * qJD(4) + t39 * qJD(5);
t227 = pkin(4) * qJD(4);
t157 = t165 - t246;
t181 = t82 * t244 + (Ifges(5,4) * t153 + (Ifges(5,1) - Ifges(5,2)) * t154) * t153 - t259 * Ifges(5,4) + t183;
t3 = t157 * t130 + t134 * t292 + t181 + t289;
t219 = t3 * qJD(1);
t7 = t183 + t289;
t218 = t7 * qJD(1);
t185 = t129 * t222 + t198 * t224 + (t153 ^ 2 + t259) * mrSges(5,3) + t264;
t199 = t209 * t164;
t18 = m(6) * (-t129 * t272 + t198 * t53) + m(5) * (-t109 * t154 + t110 * t153) + m(4) * t199 + t185;
t216 = qJD(1) * t18;
t184 = t94 * t224 - t93 * t222 + (t138 * t154 + t139 * t153) * mrSges(5,3);
t11 = m(6) * (t272 * t93 + t53 * t94) + m(5) * (-t109 * t138 + t110 * t139) + m(4) * (t199 * t178 - t245) * pkin(1) + t184 + (m(5) * t157 + m(6) * t134 + t260) * t242 + (-m(4) * t242 + t278) * t246;
t213 = t11 * qJD(1);
t205 = t176 * t224;
t201 = (t134 + t137) * t154;
t179 = t287 + (t157 / 0.2e1 + t165 / 0.2e1) * t130 + t181;
t190 = m(6) * (t173 * t94 + t176 * t93);
t1 = t179 + 0.2e1 * (-t190 / 0.4e1 + m(6) * t201 / 0.4e1) * pkin(4) + t284;
t6 = t165 * t130 + t137 * t292 + t181 + t288;
t197 = t1 * qJD(1) + t6 * qJD(2);
t182 = t183 + t287;
t4 = t182 - t229;
t8 = t183 + t288;
t196 = t4 * qJD(1) + t8 * qJD(2);
t195 = -pkin(4) * t173 * t222 + Ifges(5,5) * t153 - Ifges(5,6) * t154 + t298;
t19 = m(6) * (-t129 * t273 + t198 * t77) + m(5) * (-t132 * t154 + t133 * t153) + m(4) * t200 + t185;
t180 = -m(4) * (t199 + t200) / 0.2e1 - m(5) * ((-t109 - t132) * t154 + (t110 + t133) * t153) / 0.2e1 - m(6) * ((t53 + t77) * t198 - (t272 + t273) * t129) / 0.2e1 - t185;
t189 = (t257 + t258 + m(4) / 0.2e1) * t242;
t9 = t189 + t180;
t194 = qJD(1) * t9 - qJD(2) * t19;
t156 = (mrSges(6,1) * t173 + mrSges(6,2) * t176) * pkin(4);
t187 = qJD(4) * t156;
t155 = t156 * qJD(5);
t64 = t67 * qJD(3);
t36 = t39 * qJD(3);
t10 = t189 - t180;
t5 = t182 + t229;
t2 = t179 - t284 + (t201 * t257 + t190 / 0.2e1) * pkin(4);
t12 = [qJD(2) * t11 + qJD(3) * t18 + qJD(4) * t3 + qJD(5) * t7, t10 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t213 + (t184 + 0.2e1 * (t273 * t93 + t77 * t94) * t257 + 0.2e1 * (-t132 * t138 + t133 * t139) * t258 + (t278 * t178 + m(4) * (t178 * t200 - t245) + (m(5) * t165 + m(6) * t137 + t260) * t175) * pkin(1)) * qJD(2), qJD(2) * t10 + t216 + t230, t219 + t2 * qJD(2) + t64 + (-t110 * mrSges(5,1) - t109 * mrSges(5,2) + t195 + t295) * qJD(4) + t300 + (-t205 + m(6) * (t173 * t272 - t176 * t53)) * t227, t5 * qJD(2) + t13 * qJD(4) + t218 + t300 + t36; -qJD(3) * t9 + qJD(4) * t1 + qJD(5) * t4 - t213, qJD(3) * t19 + qJD(4) * t6 + qJD(5) * t8, -t194 + t230, t64 + (-t133 * mrSges(5,1) - t132 * mrSges(5,2) + t195 + t296) * qJD(4) + t299 + (-t205 + m(6) * (t173 * t273 - t176 * t77)) * t227 + t197, t16 * qJD(4) + t196 + t299 + t36; qJD(2) * t9 - t216 + t231, t194 + t231, 0, t266, t265; -qJD(2) * t1 - t219 - t293, -t197 - t293, -t266, -t155, -t155 - t187; -qJD(2) * t4 - t218 - t294, -t196 - t294, -t265, t187, 0;];
Cq = t12;
