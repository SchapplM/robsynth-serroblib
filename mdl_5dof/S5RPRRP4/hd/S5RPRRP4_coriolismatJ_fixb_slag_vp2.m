% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:20
% EndTime: 2020-01-03 11:49:28
% DurationCPUTime: 2.65s
% Computational Cost: add. (5778->261), mult. (12661->347), div. (0->0), fcn. (12427->6), ass. (0->137)
t247 = mrSges(5,1) + mrSges(6,1);
t250 = mrSges(5,2) + mrSges(6,2);
t154 = sin(qJ(4));
t155 = sin(qJ(3));
t156 = cos(qJ(3));
t230 = cos(qJ(4));
t137 = -t154 * t156 - t230 * t155;
t153 = cos(pkin(8));
t124 = t137 * t153;
t183 = t230 * pkin(3);
t147 = t183 + pkin(4);
t198 = t153 * t156;
t199 = t153 * t155;
t136 = -t154 * t155 + t230 * t156;
t125 = t136 * t153;
t203 = t125 * t154;
t235 = m(5) * pkin(3);
t236 = m(6) / 0.2e1;
t239 = t247 * t124 / 0.2e1 - t250 * t125 / 0.2e1;
t249 = (pkin(3) * t203 + t147 * t124) * t236 + (t230 * t124 + t203) * t235 / 0.2e1 - mrSges(4,1) * t199 / 0.2e1 - mrSges(4,2) * t198 / 0.2e1 + t239;
t152 = sin(pkin(8));
t122 = t136 * t152;
t100 = -mrSges(6,1) * t153 - t122 * mrSges(6,3);
t101 = -mrSges(5,1) * t153 - t122 * mrSges(5,3);
t243 = t100 + t101;
t248 = t243 / 0.2e1;
t246 = mrSges(6,3) + mrSges(5,3);
t245 = Ifges(6,4) + Ifges(5,4);
t123 = t137 * t152;
t216 = t123 * mrSges(6,3);
t98 = mrSges(6,2) * t153 + t216;
t217 = t123 * mrSges(5,3);
t99 = mrSges(5,2) * t153 + t217;
t244 = t98 + t99;
t114 = t122 * qJ(5);
t138 = -pkin(2) * t153 - pkin(6) * t152 - pkin(1);
t132 = t156 * t138;
t200 = t152 * t156;
t187 = pkin(7) * t200;
t84 = -t187 + t132 + (-qJ(2) * t155 - pkin(3)) * t153;
t109 = qJ(2) * t198 + t155 * t138;
t201 = t152 * t155;
t95 = -pkin(7) * t201 + t109;
t87 = t154 * t95;
t44 = t230 * t84 - t87;
t29 = -t114 + t44;
t241 = t250 * t230;
t170 = (Ifges(5,5) + Ifges(6,5)) * t123 + (-Ifges(5,6) - Ifges(6,6)) * t122;
t228 = pkin(3) * t154;
t240 = t246 * t122 * t228;
t108 = -qJ(2) * t199 + t132;
t94 = t108 - t187;
t54 = t230 * t94 - t87;
t220 = t54 * mrSges(5,2);
t89 = t230 * t95;
t53 = -t154 * t94 - t89;
t221 = t53 * mrSges(5,1);
t36 = -t114 + t54;
t224 = t36 * mrSges(6,2);
t205 = t123 * qJ(5);
t35 = t53 - t205;
t225 = t35 * mrSges(6,1);
t238 = t225 / 0.2e1 - t224 / 0.2e1 + t221 / 0.2e1 - t220 / 0.2e1 + (t35 * t236 - t216 / 0.2e1) * pkin(4);
t166 = -t250 * t136 + t247 * t137;
t237 = m(5) / 0.2e1;
t234 = m(6) * pkin(4);
t229 = m(6) * t137;
t227 = t29 * mrSges(6,2);
t45 = t154 * t84 + t89;
t30 = t45 + t205;
t226 = t30 * mrSges(6,1);
t223 = t44 * mrSges(5,2);
t222 = t45 * mrSges(5,1);
t25 = -pkin(4) * t153 + t29;
t219 = -t25 + t29;
t133 = t153 * mrSges(4,2) - mrSges(4,3) * t201;
t134 = -t153 * mrSges(4,1) - mrSges(4,3) * t200;
t135 = pkin(3) * t201 + t152 * qJ(2);
t115 = t123 * mrSges(6,2);
t70 = t122 * mrSges(6,1) + t115;
t85 = -t123 * pkin(4) + t135;
t160 = t135 * (mrSges(5,1) * t122 + mrSges(5,2) * t123) + t85 * t70 + ((-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * t153 + t245 * t123) * t123 - t25 * t216 - t44 * t217 - t170 * t153 / 0.2e1;
t163 = t30 * mrSges(6,3) + t45 * mrSges(5,3) + t245 * t122 + (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t153 + (-Ifges(6,1) + Ifges(6,2) - Ifges(5,1) + Ifges(5,2)) * t123;
t72 = -mrSges(6,1) * t123 + mrSges(6,2) * t122;
t73 = -mrSges(5,1) * t123 + mrSges(5,2) * t122;
t96 = pkin(3) * t200 + t122 * pkin(4);
t1 = m(5) * (t44 * t53 + t45 * t54) + t160 - t163 * t122 + m(6) * (t25 * t35 + t30 * t36 + t85 * t96) + t108 * t133 - t109 * t134 + t96 * t72 + t36 * t98 + t54 * t99 + t35 * t100 + t53 * t101 + ((t108 * mrSges(4,3) + Ifges(4,5) * t153 + (-mrSges(4,2) * qJ(2) + Ifges(4,4) * t155) * t152) * t155 + (-t109 * mrSges(4,3) + Ifges(4,6) * t153 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t156 + (-Ifges(4,1) + Ifges(4,2)) * t155) * t152 + (m(5) * t135 + t73) * pkin(3)) * t156) * t152;
t218 = t1 * qJD(1);
t2 = t29 * t98 + t44 * t99 - t45 * t101 + (m(6) * t219 - t100) * t30 - ((-m(6) * t85 - t72) * pkin(4) + t163) * t122 + t160;
t211 = t2 * qJD(1);
t195 = t156 * t133;
t196 = t155 * t134;
t23 = t137 * t25;
t161 = t23 * t236 - t196 / 0.2e1 + t195 / 0.2e1 + ((t44 - t54) * t237 - t36 * t236) * t137 + ((t45 + t53) * t237 + (t30 + t35) * t236) * t136;
t6 = t161 + t244 * t136 / 0.2e1 + t137 * t248 + t246 * (t137 * t122 / 0.2e1 - t136 * t123 / 0.2e1) - t249;
t210 = t6 * qJD(1);
t180 = -mrSges(6,3) / 0.2e1 - mrSges(5,3) / 0.2e1;
t162 = (t100 / 0.2e1 + t101 / 0.2e1 - t180 * t122) * t137 + (t98 / 0.2e1 + t99 / 0.2e1 + t180 * t123) * t136;
t158 = (-t137 * t29 + t23) * t236 + t162;
t8 = (mrSges(5,2) / 0.2e1 + mrSges(6,2) / 0.2e1) * t125 + (-t234 / 0.2e1 - mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t124 + t158;
t209 = t8 * qJD(1);
t150 = t152 ^ 2;
t148 = t150 * qJ(2);
t151 = t153 ^ 2;
t167 = t155 * mrSges(4,1) + t156 * mrSges(4,2);
t9 = (t195 - t196) * t153 + t244 * t125 + t243 * t124 + (t151 + t150) * mrSges(3,3) + (t167 * t152 + t72 + t73) * t152 + m(6) * (t124 * t25 + t125 * t30 + t152 * t85) + m(5) * (t124 * t44 + t125 * t45 + t135 * t152) + m(4) * (t148 + (-t108 * t155 + t109 * t156) * t153) + m(3) * (qJ(2) * t151 + t148);
t208 = t9 * qJD(1);
t15 = m(6) * (-t122 * t25 + t123 * t30) + t123 * t98 - t122 * t100;
t207 = qJD(1) * t15;
t206 = t122 * t136;
t204 = t123 * t137;
t202 = t147 * t122;
t185 = t123 * t228;
t18 = 0.2e1 * (t185 / 0.4e1 - t202 / 0.4e1 - t96 / 0.4e1) * m(6) - t70;
t194 = t18 * qJD(1);
t41 = (-t206 / 0.2e1 - t204 / 0.2e1 - t152 / 0.2e1) * m(6);
t193 = t41 * qJD(1);
t179 = mrSges(6,1) + t234;
t56 = -t179 * t122 - t115;
t192 = t56 * qJD(1);
t188 = t234 / 0.2e1;
t182 = t147 * t216;
t169 = t183 - t147;
t168 = t183 * t217;
t26 = (pkin(4) / 0.2e1 + t183 / 0.2e1 - t147 / 0.2e1) * t229;
t157 = -m(6) * (-t147 * t30 + (t219 * t154 + t230 * t30) * pkin(3)) / 0.2e1 + t227 / 0.2e1 + t226 / 0.2e1 + t223 / 0.2e1 + t222 / 0.2e1 + t182 / 0.2e1 + t168 / 0.2e1 + t228 * t248 + t240 / 0.2e1 - t244 * t183 / 0.2e1;
t4 = t157 + t238;
t68 = t241 * pkin(3) + (-m(6) * t169 + t247) * t228;
t165 = t4 * qJD(1) + t26 * qJD(2) + t68 * qJD(3);
t128 = t136 * t228;
t40 = (-t204 - t206 + t152) * t236;
t37 = (t185 - t202 + t96) * t236;
t17 = -t169 * t229 / 0.2e1 + t137 * t188 + t166;
t7 = t124 * t188 + t158 + t239;
t5 = t161 + t162 + t249;
t3 = -t157 + t170 + t238;
t10 = [qJD(2) * t9 + qJD(3) * t1 + qJD(4) * t2 + qJD(5) * t15, t208 + 0.2e1 * (t236 + t237) * (t124 * t136 - t125 * t137) * qJD(2) + t5 * qJD(3) + t7 * qJD(4) + t40 * qJD(5), t218 + t5 * qJD(2) + (-Ifges(4,5) * t201 - Ifges(4,6) * t200 + m(6) * (t147 * t35 + t36 * t228) + t225 - t220 + t221 - t224 - t182 + (t154 * t54 + t230 * t53) * t235 - t109 * mrSges(4,1) - t108 * mrSges(4,2) - t168 + t170 - t240) * qJD(3) + t3 * qJD(4) + t37 * qJD(5), t211 + t7 * qJD(2) + t3 * qJD(3) + (-t222 - t226 - t223 - t227 + (-m(6) * t30 - t216) * pkin(4) + t170) * qJD(4), qJD(2) * t40 + qJD(3) * t37 + t207; qJD(3) * t6 + qJD(4) * t8 + qJD(5) * t41 - t208, 0, t210 + (-t167 + t166) * qJD(3) + t17 * qJD(4) + 0.2e1 * ((t137 * t183 + t128) * t237 + (t147 * t137 + t128) * t236) * qJD(3), t209 + t17 * qJD(3) + (pkin(4) * t229 + t166) * qJD(4), t193; -qJD(2) * t6 - qJD(4) * t4 + qJD(5) * t18 - t218, -qJD(4) * t26 - t210, -t68 * qJD(4), ((-mrSges(5,1) - t179) * t154 - t241) * qJD(4) * pkin(3) - t165, t194; -qJD(2) * t8 + qJD(3) * t4 + qJD(5) * t56 - t211, qJD(3) * t26 - t209, t165, 0, t192; -qJD(2) * t41 - qJD(3) * t18 - qJD(4) * t56 - t207, -t193, -t194, -t192, 0;];
Cq = t10;
