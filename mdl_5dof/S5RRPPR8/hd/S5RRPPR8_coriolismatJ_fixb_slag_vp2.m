% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:10
% EndTime: 2019-12-31 19:38:14
% DurationCPUTime: 1.97s
% Computational Cost: add. (4539->172), mult. (8315->225), div. (0->0), fcn. (8809->6), ass. (0->91)
t138 = sin(pkin(8));
t140 = sin(qJ(2));
t142 = cos(qJ(2));
t171 = cos(pkin(8));
t108 = -t140 * t138 - t142 * t171;
t109 = -t142 * t138 + t140 * t171;
t139 = sin(qJ(5));
t141 = cos(qJ(5));
t151 = t108 * t139 + t109 * t141;
t156 = t141 * t108 - t109 * t139;
t206 = pkin(6) - qJ(4);
t119 = t206 * t140;
t122 = t206 * t142;
t202 = t138 * t119 + t171 * t122;
t209 = pkin(7) * t108 + t202;
t91 = -t171 * t119 + t138 * t122;
t60 = -pkin(7) * t109 - t91;
t215 = t139 * t60 + t141 * t209;
t173 = t141 * t60;
t174 = t139 * t209;
t226 = -t173 + t174;
t235 = -t215 * mrSges(6,1) + t226 * mrSges(6,2) + Ifges(6,5) * t156 - Ifges(6,6) * t151;
t237 = t235 * qJD(5);
t201 = -t142 * pkin(2) - t140 * qJ(3);
t118 = -pkin(1) + t201;
t104 = t142 * pkin(3) - t118;
t87 = -pkin(4) * t108 + t104;
t236 = m(6) * t87 - mrSges(6,1) * t156 + mrSges(6,2) * t151;
t232 = m(5) * t104 - mrSges(5,1) * t108 + mrSges(5,2) * t109;
t204 = t151 * mrSges(6,1);
t219 = t156 * mrSges(6,2);
t230 = -Ifges(6,4) * t156 ^ 2 + t87 * (-t219 - t204);
t106 = -t139 * t138 + t141 * t171;
t107 = -t141 * t138 - t139 * t171;
t229 = t106 * t215 - t107 * t226;
t184 = t156 * mrSges(6,3);
t217 = t219 / 0.2e1 + t204 / 0.2e1;
t223 = Ifges(6,4) * t151 + (-Ifges(6,1) + Ifges(6,2)) * t156;
t143 = -pkin(2) - pkin(3);
t116 = -qJ(3) * t138 + t171 * t143;
t117 = t171 * qJ(3) + t138 * t143;
t197 = m(6) / 0.2e1;
t199 = m(5) / 0.2e1;
t115 = -pkin(4) + t116;
t84 = t115 * t141 - t117 * t139;
t85 = t115 * t139 + t117 * t141;
t218 = (t117 * t108 - t116 * t109) * t199 + (-t151 * t84 + t156 * t85) * t197 + t217;
t216 = t138 * t91 + t171 * t202;
t214 = t106 / 0.2e1;
t157 = -t107 * mrSges(6,1) + t106 * mrSges(6,2);
t208 = t157 * qJD(5);
t207 = Ifges(3,4) - Ifges(4,5);
t200 = -m(5) / 0.2e1;
t198 = -m(6) / 0.2e1;
t196 = -t107 / 0.2e1;
t180 = t151 * mrSges(6,3);
t8 = -t151 * t180 - t156 * t184 - m(6) * (t151 * t226 + t156 * t215) - m(5) * (t108 * t202 + t109 * t91) + (-t108 ^ 2 - t109 ^ 2) * mrSges(5,3);
t176 = qJD(1) * t8;
t133 = t142 * qJ(3);
t105 = t143 * t140 + t133;
t120 = -t142 * mrSges(4,1) - t140 * mrSges(4,3);
t160 = m(4) * t118 + t120;
t88 = -pkin(4) * t109 + t105;
t1 = t160 * (pkin(2) * t140 - t133) + (-pkin(1) * mrSges(3,2) - t118 * mrSges(4,3) + t207 * t142) * t142 + (-t104 * mrSges(5,1) + Ifges(5,4) * t109) * t109 + (-t104 * mrSges(5,2) - Ifges(5,4) * t108 + (-Ifges(5,1) + Ifges(5,2)) * t109) * t108 + (-pkin(1) * mrSges(3,1) + t118 * mrSges(4,1) + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t142 - t207 * t140) * t140 + t223 * t151 + t232 * t105 + t230 + t236 * t88;
t175 = t1 * qJD(1);
t2 = t215 * t180 + (-t215 * mrSges(6,3) + t223) * t151 + t230;
t172 = t2 * qJD(1);
t20 = (-t160 + t232 + t236) * t140;
t170 = qJD(1) * t20;
t19 = t85 * mrSges(6,1) + mrSges(6,2) * t84;
t169 = qJD(5) * t19;
t148 = t105 * t199 + t88 * t197 - t217;
t11 = -t109 * mrSges(5,1) - t108 * mrSges(5,2) + t148 - t218;
t167 = t11 * qJD(1);
t166 = t157 * qJD(2);
t145 = (-t106 * t151 - t107 * t156) * t197 + (t138 * t108 - t171 * t109) * t199;
t25 = (t198 + t200) * t140 + t145;
t165 = t25 * qJD(1);
t26 = 0.2e1 * t217;
t164 = t26 * qJD(1);
t3 = (-t226 / 0.2e1 + t174 / 0.2e1 - t173 / 0.2e1) * mrSges(6,2);
t153 = -t3 * qJD(1) + t19 * qJD(2);
t23 = m(4) * qJ(3) + t138 * mrSges(5,1) + t171 * mrSges(5,2) + mrSges(4,3) + m(6) * (t85 * t106 + t84 * t107) + m(5) * (-t116 * t138 + t117 * t171) + t157;
t144 = t229 * t198 + t216 * t200;
t146 = t229 * t197 + t216 * t199;
t6 = t144 + t146;
t152 = -t6 * qJD(1) + t23 * qJD(2);
t24 = t145 + (m(5) + m(6)) * t140 / 0.2e1;
t12 = t148 + t218;
t5 = t184 * t214 + t180 * t196 + (m(4) * pkin(6) + mrSges(4,2)) * t142 + (t151 * t196 + t156 * t214) * mrSges(6,3) - t144 + t146 + (t171 * t108 + t138 * t109) * mrSges(5,3);
t4 = [qJD(2) * t1 + qJD(3) * t20 - qJD(4) * t8 - qJD(5) * t2, t5 * qJD(3) + t12 * qJD(4) - t237 + t175 + (-t202 * mrSges(5,1) + t91 * mrSges(5,2) + Ifges(5,5) * t108 - Ifges(5,6) * t109 + 0.2e1 * (t215 * t84 + t226 * t85) * t197 + 0.2e1 * (t116 * t202 + t117 * t91) * t199 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t142 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t140 + (m(4) * t201 - t142 * mrSges(3,1) + t140 * mrSges(3,2) + t120) * pkin(6) + (t151 * t85 + t156 * t84) * mrSges(6,3) + (t108 * t116 + t109 * t117) * mrSges(5,3) + t235) * qJD(2), qJD(2) * t5 + qJD(4) * t24 + t170, qJD(2) * t12 + qJD(3) * t24 - t176, -qJD(2) * t235 - t172 + t237; -qJD(3) * t6 - qJD(4) * t11 - qJD(5) * t3 - t175, qJD(3) * t23 + t169, t152, -t167, t153 - t169; qJD(2) * t6 + qJD(4) * t25 - t170, -t152 + t208, 0, t165, t166 - t208; qJD(2) * t11 - qJD(3) * t25 + qJD(5) * t26 + t176, t167, -t165, 0, t164; qJD(2) * t3 - qJD(4) * t26 + t172, -qJD(3) * t157 - t153, -t166, -t164, 0;];
Cq = t4;
