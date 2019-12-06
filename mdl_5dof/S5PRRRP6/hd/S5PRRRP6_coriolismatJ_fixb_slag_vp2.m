% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:51
% EndTime: 2019-12-05 16:50:55
% DurationCPUTime: 1.78s
% Computational Cost: add. (2971->183), mult. (7173->249), div. (0->0), fcn. (6681->6), ass. (0->116)
t139 = cos(qJ(2));
t213 = -t139 / 0.2e1;
t221 = m(5) * pkin(3);
t136 = sin(qJ(3));
t258 = t136 * t213;
t138 = cos(qJ(3));
t128 = -pkin(3) * t138 - pkin(2);
t211 = sin(qJ(4));
t176 = t211 * t136;
t212 = cos(qJ(4));
t113 = -t138 * t212 + t176;
t178 = t212 * t136;
t114 = -t138 * t211 - t178;
t190 = t114 * qJ(5);
t163 = t113 * pkin(4) + t190;
t51 = t128 + t163;
t66 = mrSges(6,1) * t113 + mrSges(6,3) * t114;
t238 = m(6) * t51 + t66;
t165 = (Ifges(5,6) - Ifges(6,6)) * t114 + (-Ifges(6,4) - Ifges(5,5)) * t113;
t246 = pkin(7) + pkin(6);
t185 = t246 * t138;
t225 = -t176 * t246 + t185 * t212;
t239 = t225 * mrSges(6,1);
t240 = t225 * mrSges(5,1);
t78 = t178 * t246 + t185 * t211;
t253 = t78 * mrSges(6,3);
t254 = t78 * mrSges(5,2);
t257 = t165 - t239 - t240 + t254 - t253;
t256 = t239 / 0.2e1 + t240 / 0.2e1 + t253 / 0.2e1 - t254 / 0.2e1;
t193 = t138 * mrSges(4,2);
t223 = -m(6) / 0.2e1;
t245 = t139 / 0.2e1;
t64 = -mrSges(6,1) * t114 + mrSges(6,3) * t113;
t65 = -mrSges(5,1) * t114 - mrSges(5,2) * t113;
t226 = (t64 + t65) * t245;
t210 = pkin(3) * t136;
t63 = -pkin(4) * t114 + qJ(5) * t113;
t54 = t63 + t210;
t252 = t258 * t221 + t54 * t139 * t223 - (t136 * mrSges(4,1) + t193) * t245 - t226;
t249 = -pkin(4) * t225 - qJ(5) * t78;
t241 = -Ifges(6,5) + Ifges(5,4);
t248 = ((Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t114 + t241 * t113) * t113 - t114 ^ 2 * t241 + t128 * t65 + t51 * t64;
t222 = m(6) / 0.2e1;
t247 = 0.2e1 * t222;
t243 = mrSges(5,1) + mrSges(6,1);
t242 = -mrSges(5,2) + mrSges(6,3);
t130 = m(6) * qJ(5) + mrSges(6,3);
t237 = qJD(4) * t130;
t236 = t130 * qJD(5);
t99 = t113 * t139;
t235 = (mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1) * t99;
t198 = t113 * mrSges(6,2);
t105 = -t198 / 0.2e1;
t180 = t211 * pkin(3);
t181 = t212 * pkin(3);
t229 = t113 * t181 + t114 * t180;
t228 = -t138 * mrSges(4,1) + t136 * mrSges(4,2);
t220 = m(6) * pkin(3);
t219 = -mrSges(6,2) / 0.2e1;
t98 = t114 * t139;
t218 = t98 / 0.2e1;
t137 = sin(qJ(2));
t96 = t113 * t137;
t97 = t114 * t137;
t216 = m(6) * (pkin(4) * t96 + qJ(5) * t97);
t215 = m(6) * t225;
t214 = m(6) * t96;
t201 = qJD(5) * t214;
t184 = t136 ^ 2 + t138 ^ 2;
t187 = t137 * t139;
t12 = m(4) * (-0.1e1 + t184) * t187 + (m(5) + m(6)) * (t96 * t99 + t97 * t98 - t187);
t189 = t12 * qJD(1);
t182 = t216 / 0.2e1;
t123 = t180 + qJ(5);
t179 = t123 * t114 * mrSges(6,2);
t175 = t114 * t213;
t172 = -t225 * t99 - t78 * t98;
t170 = t218 * t243 + t235;
t166 = 0.2e1 * t105;
t67 = mrSges(5,1) * t113 - mrSges(5,2) * t114;
t1 = (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t136 + pkin(3) * t67) * t136 + m(5) * t128 * t210 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t138 + (Ifges(4,1) - Ifges(4,2)) * t136) * t138 + t248 + t238 * t54;
t127 = -t181 - pkin(4);
t142 = (-t123 * t99 - t127 * t98) * t222 + (-t211 * t99 + t212 * t98) * t221 / 0.2e1 + mrSges(4,1) * t258 + t193 * t213 + t170;
t3 = t142 - t252;
t162 = -t3 * qJD(1) + t1 * qJD(2);
t2 = t238 * t63 + t248;
t144 = -t63 * t139 * t222 - t226;
t153 = m(6) * (pkin(4) * t98 - t99 * qJ(5));
t5 = -t235 - (mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1) * t98 - t153 / 0.2e1 + t144;
t161 = t5 * qJD(1) + t2 * qJD(2);
t160 = t123 * t97 - t127 * t96;
t21 = t238 * t114;
t47 = (t175 + t218) * m(6);
t157 = -qJD(1) * t47 - qJD(2) * t21;
t117 = m(6) * t123 + mrSges(6,3);
t156 = qJD(3) * t117;
t155 = -qJD(3) * t130 - t237;
t152 = t242 * t97 + t243 * t96;
t149 = t211 * t97 + t212 * t96;
t145 = m(6) * (-pkin(3) * t149 + t160);
t13 = t182 - t145 / 0.2e1;
t146 = -t180 * t243 + t181 * t242;
t32 = -(t123 * t212 + t127 * t211) * t220 - t146;
t140 = ((t180 - t123) * t78 + (t127 + t181) * t225) * t222 + t179 / 0.2e1 + t127 * t105 + t229 * t219 - t256;
t143 = pkin(4) * t105 + t190 * t219 + t249 * t223 + t256;
t8 = t140 + t143;
t148 = -t13 * qJD(1) + t8 * qJD(2) - t32 * qJD(3);
t103 = t180 * t247 + t130;
t46 = m(6) * t175 - t222 * t98;
t25 = t166 + t215;
t22 = t215 / 0.2e1 + t225 * t222 + t166;
t9 = t145 / 0.2e1 + t182 + t152;
t7 = t140 - t143 + t165;
t6 = t153 / 0.2e1 + t144 + t170;
t4 = t142 + t252;
t10 = [qJD(2) * t12, t4 * qJD(3) + t6 * qJD(4) + t46 * qJD(5) + t189 + ((mrSges(4,3) * t184 - mrSges(3,2)) * t139 + (-mrSges(3,1) + t66 + t67 + t228) * t137 + m(4) * (pkin(6) * t139 * t184 - t137 * pkin(2)) + m(5) * (t128 * t137 + t172) + (t137 * t51 + t172) * t247 + (mrSges(6,2) + mrSges(5,3)) * (t113 * t99 + t114 * t98)) * qJD(2), t4 * qJD(2) + (m(6) * t160 + t137 * t228 + t149 * t221 + t152) * qJD(3) + t9 * qJD(4) - t201, t6 * qJD(2) + t9 * qJD(3) + (t152 + t216) * qJD(4) - t201, t46 * qJD(2) + 0.2e1 * (-qJD(3) / 0.2e1 - qJD(4) / 0.2e1) * t214; -qJD(3) * t3 + qJD(4) * t5 + qJD(5) * t47 - t189, qJD(3) * t1 + qJD(4) * t2 + qJD(5) * t21, (-t127 * t198 + t179 + m(6) * (-t123 * t78 + t127 * t225) + (-t211 * t78 - t212 * t225) * t221 + Ifges(4,5) * t138 - Ifges(4,6) * t136 + t228 * pkin(6) + t229 * mrSges(5,3) + t257) * qJD(3) + t7 * qJD(4) + t22 * qJD(5) + t162, t7 * qJD(3) + (m(6) * t249 + t163 * mrSges(6,2) + t257) * qJD(4) + t25 * qJD(5) + t161, qJD(3) * t22 + qJD(4) * t25 - t157; qJD(2) * t3 - qJD(4) * t13, qJD(4) * t8 - t162, -qJD(4) * t32 + qJD(5) * t117, ((-pkin(4) * t211 + qJ(5) * t212) * t220 + t146) * qJD(4) + t103 * qJD(5) + t148, qJD(4) * t103 + t156; -qJD(2) * t5 + qJD(3) * t13, -qJD(3) * t8 - t161, -t148 + t236, t236, -t155; -t47 * qJD(2), t157, -t156 - t237, t155, 0;];
Cq = t10;
