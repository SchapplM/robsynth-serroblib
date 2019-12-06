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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:05:35
% EndTime: 2019-12-05 18:05:42
% DurationCPUTime: 1.93s
% Computational Cost: add. (5778->257), mult. (12661->338), div. (0->0), fcn. (12427->6), ass. (0->137)
t244 = mrSges(5,1) + mrSges(6,1);
t248 = mrSges(5,2) + mrSges(6,2);
t247 = mrSges(5,3) + mrSges(6,3);
t150 = sin(qJ(4));
t151 = sin(qJ(3));
t152 = cos(qJ(3));
t228 = cos(qJ(4));
t133 = -t150 * t152 - t228 * t151;
t230 = t133 / 0.2e1;
t243 = Ifges(6,4) + Ifges(5,4);
t149 = cos(pkin(8));
t148 = sin(pkin(8));
t119 = t133 * t148;
t214 = t119 * mrSges(6,3);
t94 = t149 * mrSges(6,2) + t214;
t215 = t119 * mrSges(5,3);
t95 = t149 * mrSges(5,2) + t215;
t240 = t95 + t94;
t132 = -t150 * t151 + t228 * t152;
t118 = t132 * t148;
t96 = -t149 * mrSges(6,1) - t118 * mrSges(6,3);
t97 = -t149 * mrSges(5,1) - t118 * mrSges(5,3);
t239 = t97 + t96;
t120 = t133 * t149;
t181 = t228 * pkin(3);
t143 = t181 + pkin(4);
t121 = t132 * t149;
t167 = t244 * t120 / 0.2e1 - t248 * t121 / 0.2e1;
t197 = t149 * t152;
t198 = t149 * t151;
t204 = t121 * t150;
t233 = m(5) * pkin(3);
t234 = m(6) / 0.2e1;
t242 = -(pkin(3) * t204 + t143 * t120) * t234 - (t228 * t120 + t204) * t233 / 0.2e1 + mrSges(4,1) * t198 / 0.2e1 + mrSges(4,2) * t197 / 0.2e1 - t167;
t110 = t118 * qJ(5);
t134 = -t149 * pkin(2) - t148 * pkin(6) - pkin(1);
t128 = t152 * t134;
t199 = t148 * t152;
t183 = pkin(7) * t199;
t80 = -t183 + t128 + (-qJ(2) * t151 - pkin(3)) * t149;
t105 = qJ(2) * t197 + t151 * t134;
t200 = t148 * t151;
t91 = -pkin(7) * t200 + t105;
t83 = t150 * t91;
t40 = t228 * t80 - t83;
t27 = -t110 + t40;
t23 = -t149 * pkin(4) + t27;
t217 = t23 - t27;
t176 = m(6) * t217;
t238 = t248 * t228;
t166 = (Ifges(5,5) + Ifges(6,5)) * t119 + (-Ifges(5,6) - Ifges(6,6)) * t118;
t226 = pkin(3) * t150;
t237 = t247 * t118 * t226;
t104 = -qJ(2) * t198 + t128;
t90 = t104 - t183;
t50 = t228 * t90 - t83;
t218 = t50 * mrSges(5,2);
t85 = t228 * t91;
t49 = -t150 * t90 - t85;
t219 = t49 * mrSges(5,1);
t34 = -t110 + t50;
t222 = t34 * mrSges(6,2);
t205 = t119 * qJ(5);
t33 = t49 - t205;
t223 = t33 * mrSges(6,1);
t236 = t223 / 0.2e1 - t222 / 0.2e1 + t219 / 0.2e1 - t218 / 0.2e1 + (t33 * t234 - t214 / 0.2e1) * pkin(4);
t162 = -t248 * t132 + t244 * t133;
t235 = m(5) / 0.2e1;
t232 = m(6) * pkin(4);
t227 = m(6) * t133;
t225 = t27 * mrSges(6,2);
t41 = t150 * t80 + t85;
t28 = t41 + t205;
t224 = t28 * mrSges(6,1);
t221 = t40 * mrSges(5,2);
t220 = t41 * mrSges(5,1);
t129 = t149 * mrSges(4,2) - mrSges(4,3) * t200;
t130 = -t149 * mrSges(4,1) - mrSges(4,3) * t199;
t131 = pkin(3) * t200 + t148 * qJ(2);
t111 = t119 * mrSges(6,2);
t66 = t118 * mrSges(6,1) + t111;
t81 = -t119 * pkin(4) + t131;
t156 = t131 * (t118 * mrSges(5,1) + t119 * mrSges(5,2)) + t81 * t66 + ((-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * t149 + t243 * t119) * t119 - t23 * t214 - t40 * t215 - t166 * t149 / 0.2e1;
t157 = t41 * mrSges(5,3) + t28 * mrSges(6,3) + t243 * t118 + (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t149 + (-Ifges(6,1) - Ifges(5,1) + Ifges(6,2) + Ifges(5,2)) * t119;
t68 = -t119 * mrSges(6,1) + t118 * mrSges(6,2);
t69 = -t119 * mrSges(5,1) + t118 * mrSges(5,2);
t92 = pkin(3) * t199 + t118 * pkin(4);
t1 = m(6) * (t23 * t33 + t28 * t34 + t81 * t92) + ((t104 * mrSges(4,3) + Ifges(4,5) * t149 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t151) * t148) * t151 + (-t105 * mrSges(4,3) + Ifges(4,6) * t149 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t152 + (-Ifges(4,1) + Ifges(4,2)) * t151) * t148 + (m(5) * t131 + t69) * pkin(3)) * t152) * t148 + m(5) * (t40 * t49 + t41 * t50) - t157 * t118 - t105 * t130 + t104 * t129 + t50 * t95 + t33 * t96 + t49 * t97 + t92 * t68 + t34 * t94 + t156;
t216 = t1 * qJD(1);
t2 = t40 * t95 - t41 * t97 + t27 * t94 + (-t96 - t176) * t28 - ((-m(6) * t81 - t68) * pkin(4) + t157) * t118 + t156;
t209 = t2 * qJD(1);
t161 = t240 * t132 / 0.2e1 + t239 * t230 + t247 * (t118 * t230 - t132 * t119 / 0.2e1);
t193 = t152 * t129;
t194 = t151 * t130;
t154 = ((t40 - t50) * t133 + (t41 + t49) * t132) * t235 + ((t23 - t34) * t133 + (t28 + t33) * t132) * t234 - t194 / 0.2e1 + t193 / 0.2e1 + t161;
t6 = t154 + t242;
t208 = t6 * qJD(1);
t186 = t232 / 0.2e1;
t160 = t120 * t186 + t167;
t178 = mrSges(5,3) / 0.2e1 + mrSges(6,3) / 0.2e1;
t7 = (-t95 / 0.2e1 - t94 / 0.2e1 + t178 * t119) * t132 + (-t97 / 0.2e1 - t176 / 0.2e1 - t96 / 0.2e1 - t178 * t118) * t133 + t160;
t207 = t7 * qJD(1);
t146 = t148 ^ 2;
t144 = t146 * qJ(2);
t147 = t149 ^ 2;
t163 = t151 * mrSges(4,1) + t152 * mrSges(4,2);
t9 = (t193 - t194) * t149 + t240 * t121 + t239 * t120 + (t147 + t146) * mrSges(3,3) + (t163 * t148 + t68 + t69) * t148 + m(5) * (t40 * t120 + t41 * t121 + t131 * t148) + m(6) * (t23 * t120 + t28 * t121 + t81 * t148) + m(4) * (t144 + (-t104 * t151 + t105 * t152) * t149) + m(3) * (t147 * qJ(2) + t144);
t206 = t9 * qJD(1);
t203 = t132 * t118;
t202 = t133 * t119;
t201 = t143 * t118;
t15 = m(6) * (-t23 * t118 + t28 * t119) - t118 * t96 + t119 * t94;
t196 = t15 * qJD(1);
t184 = t119 * t226;
t18 = 0.2e1 * (t184 / 0.4e1 - t201 / 0.4e1 - t92 / 0.4e1) * m(6) - t66;
t192 = t18 * qJD(1);
t37 = (-t203 / 0.2e1 - t202 / 0.2e1 - t148 / 0.2e1) * m(6);
t191 = t37 * qJD(1);
t177 = mrSges(6,1) + t232;
t52 = -t177 * t118 - t111;
t190 = t52 * qJD(1);
t180 = t143 * t214;
t165 = t181 - t143;
t164 = t181 * t215;
t24 = (pkin(4) / 0.2e1 + t181 / 0.2e1 - t143 / 0.2e1) * t227;
t153 = -m(6) * (-t143 * t28 + (-t217 * t150 + t228 * t28) * pkin(3)) / 0.2e1 + t225 / 0.2e1 + t224 / 0.2e1 + t221 / 0.2e1 + t220 / 0.2e1 + t180 / 0.2e1 + t164 / 0.2e1 + t239 * t226 / 0.2e1 + t237 / 0.2e1 - t240 * t181 / 0.2e1;
t4 = t153 + t236;
t64 = t238 * pkin(3) + (-m(6) * t165 + t244) * t226;
t159 = t4 * qJD(1) + t24 * qJD(2) + t64 * qJD(3);
t124 = t132 * t226;
t36 = (-t202 - t203 + t148) * t234;
t35 = (t184 - t201 + t92) * t234;
t17 = -t165 * t227 / 0.2e1 + t133 * t186 + t162;
t8 = t176 * t230 + t160 + t161;
t5 = t154 - t242;
t3 = -t153 + t166 + t236;
t10 = [qJD(2) * t9 + qJD(3) * t1 + qJD(4) * t2 + qJD(5) * t15, t206 + 0.2e1 * (t235 + t234) * (t132 * t120 - t133 * t121) * qJD(2) + t5 * qJD(3) + t8 * qJD(4) + t36 * qJD(5), t216 + t5 * qJD(2) + (-Ifges(4,5) * t200 - Ifges(4,6) * t199 + (t150 * t50 + t228 * t49) * t233 + m(6) * (t143 * t33 + t34 * t226) - t180 - t222 + t223 - t218 + t219 - t104 * mrSges(4,2) - t105 * mrSges(4,1) - t164 + t166 - t237) * qJD(3) + t3 * qJD(4) + t35 * qJD(5), t209 + t8 * qJD(2) + t3 * qJD(3) + (-t220 - t224 - t221 - t225 + (-m(6) * t28 - t214) * pkin(4) + t166) * qJD(4), t36 * qJD(2) + t35 * qJD(3) + t196; t6 * qJD(3) - t7 * qJD(4) + t37 * qJD(5) - t206, 0, t208 + (-t163 + t162) * qJD(3) + t17 * qJD(4) + 0.2e1 * ((t133 * t181 + t124) * t235 + (t143 * t133 + t124) * t234) * qJD(3), -t207 + t17 * qJD(3) + (pkin(4) * t227 + t162) * qJD(4), t191; -t6 * qJD(2) - t4 * qJD(4) + t18 * qJD(5) - t216, -t24 * qJD(4) - t208, -t64 * qJD(4), ((-mrSges(5,1) - t177) * t150 - t238) * qJD(4) * pkin(3) - t159, t192; t7 * qJD(2) + t4 * qJD(3) + t52 * qJD(5) - t209, t24 * qJD(3) + t207, t159, 0, t190; -t37 * qJD(2) - t18 * qJD(3) - t52 * qJD(4) - t196, -t191, -t192, -t190, 0;];
Cq = t10;
