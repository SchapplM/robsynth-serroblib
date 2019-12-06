% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:44
% EndTime: 2019-12-05 17:07:47
% DurationCPUTime: 1.55s
% Computational Cost: add. (3636->174), mult. (7518->233), div. (0->0), fcn. (6798->6), ass. (0->113)
t119 = sin(qJ(5));
t120 = sin(qJ(4));
t122 = cos(qJ(5));
t123 = cos(qJ(4));
t101 = -t119 * t120 + t122 * t123;
t102 = -t119 * t123 - t122 * t120;
t165 = Ifges(6,4) * t102;
t231 = Ifges(6,2) - Ifges(6,1);
t236 = t101 * t231 - t165;
t124 = cos(qJ(3));
t192 = t124 * pkin(2);
t82 = t102 * t192;
t83 = t101 * t192;
t29 = t101 * t82 - t102 * t83;
t235 = m(6) * t29;
t121 = sin(qJ(3));
t193 = t121 * pkin(2);
t114 = pkin(7) + t193;
t177 = pkin(8) + t114;
t98 = t177 * t120;
t99 = t177 * t123;
t148 = -t119 * t99 - t122 * t98;
t68 = -t119 * t98 + t122 * t99;
t234 = -t68 * mrSges(6,1) - t148 * mrSges(6,2);
t200 = -pkin(8) - pkin(7);
t110 = t200 * t120;
t111 = t200 * t123;
t146 = t122 * t110 + t111 * t119;
t79 = t110 * t119 - t111 * t122;
t233 = -t79 * mrSges(6,1) - t146 * mrSges(6,2);
t172 = Ifges(6,5) * t101 + Ifges(6,6) * t102;
t11 = t172 + t234;
t230 = t11 * qJD(5);
t15 = t172 + t233;
t229 = t15 * qJD(5);
t226 = t29 / 0.4e1;
t115 = -t123 * pkin(4) - pkin(3);
t106 = t115 - t192;
t225 = m(6) * (t106 + t115) * t120;
t195 = pkin(4) * t120;
t224 = m(6) * t195;
t222 = t68 * mrSges(6,3);
t220 = t79 * mrSges(6,3);
t219 = Ifges(5,1) - Ifges(5,2);
t117 = t120 ^ 2;
t156 = t123 ^ 2 + t117;
t218 = t156 * mrSges(5,3) - mrSges(4,2);
t71 = -t102 * mrSges(6,1) + t101 * mrSges(6,2);
t217 = t71 * qJD(5);
t107 = -mrSges(5,1) * t123 + mrSges(5,2) * t120;
t173 = t82 * mrSges(6,1) / 0.2e1 - t83 * mrSges(6,2) / 0.2e1;
t209 = (t101 * t83 + t102 * t82) * mrSges(6,3);
t161 = t102 * mrSges(6,3);
t41 = t68 * t161;
t54 = t106 * t71;
t56 = t79 * t161;
t64 = t115 * t71;
t208 = t41 / 0.2e1 + t54 / 0.2e1 + t56 / 0.2e1 + t64 / 0.2e1;
t207 = 0.2e1 * m(6);
t204 = t68 / 0.2e1;
t203 = -t68 / 0.2e1;
t202 = t79 / 0.2e1;
t201 = -t79 / 0.2e1;
t197 = pkin(3) * t121;
t196 = pkin(4) * t119;
t175 = -t41 - t54;
t174 = -t56 - t64;
t169 = Ifges(5,4) * t120;
t168 = Ifges(5,4) * t123;
t167 = Ifges(6,4) * t101;
t166 = Ifges(6,4) * t101 ^ 2;
t164 = pkin(4) * qJD(4);
t163 = t101 * mrSges(6,3);
t159 = t120 * mrSges(5,1);
t158 = t123 * mrSges(5,2);
t143 = -t101 * mrSges(6,1) - t102 * mrSges(6,2);
t10 = -mrSges(4,1) * t193 + t143 * t193 + t107 * t193 + m(6) * (t106 * t193 + t148 * t82 + t68 * t83) + m(5) * (-t197 + (t156 * t114 - t193) * t124) * pkin(2) + t209 + t218 * t192;
t157 = t10 * qJD(2);
t154 = qJD(2) * t235;
t153 = t122 * t163;
t70 = t143 * t195;
t149 = t117 * Ifges(5,4) - t70;
t145 = -t192 / 0.2e1 + pkin(3) / 0.2e1;
t144 = t158 + t159;
t40 = t148 * t163;
t55 = t146 * t163;
t90 = (-pkin(3) - t192) * t144;
t142 = -t40 / 0.2e1 - t55 / 0.2e1 - t70 - t90 / 0.2e1;
t141 = Ifges(5,5) * t123 - Ifges(5,6) * t120 + t161 * t196 + t172;
t140 = m(6) * (t119 * t83 + t122 * t82);
t138 = -t219 * t120 - t168;
t137 = -t102 * t231 - t167;
t5 = -t40 - t90 - t106 * t224 + (t165 + t222) * t102 + t138 * t123 + (mrSges(6,3) * t148 + t137) * t101 + t149 + t175;
t133 = t5 * qJD(2);
t8 = -t166 + (-t236 + t222) * t102 + t175;
t132 = t8 * qJD(2);
t131 = t167 * t101 + t236 * t102 + t208;
t1 = (t145 * mrSges(5,1) + t169) * t120 + 0.2e1 * (t140 / 0.4e1 - t225 / 0.4e1) * pkin(4) + (t165 + (t202 + t204) * mrSges(6,3)) * t102 + (t145 * mrSges(5,2) + t138) * t123 + ((t146 / 0.2e1 + t148 / 0.2e1) * mrSges(6,3) + t137) * t101 + t142 + t173 - t208;
t6 = -t55 + pkin(3) * t159 - t115 * t224 + (t165 + t220) * t102 + (pkin(3) * mrSges(5,2) + t138) * t123 + (mrSges(6,3) * t146 + t137) * t101 + t149 + t174;
t129 = t1 * qJD(2) + t6 * qJD(3);
t125 = (t203 + t201) * t161 + t131;
t4 = t125 - t173;
t9 = -t166 + (-t236 + t220) * t102 + t174;
t128 = -t4 * qJD(2) + t9 * qJD(3);
t105 = (mrSges(6,1) * t119 + mrSges(6,2) * t122) * pkin(4);
t13 = (t203 + t204) * mrSges(6,1);
t16 = (t201 + t202) * mrSges(6,1);
t126 = -t13 * qJD(2) - t16 * qJD(3) + t105 * qJD(4);
t100 = t105 * qJD(5);
t7 = qJD(3) * t226 * t207;
t3 = t125 + t173;
t2 = (-t158 / 0.2e1 - t159 / 0.2e1) * t192 - pkin(3) * t144 / 0.2e1 + t168 * t123 + t131 - t142 + t173 - (t146 + t148) * t163 / 0.2e1 + (-t79 - t68) * t161 / 0.2e1 + (t123 * t219 - t169) * t120 + (t225 + t140) * pkin(4) / 0.2e1;
t12 = [0, t7, t154 / 0.2e1, -t217 + (-t144 - t71 + (t122 * pkin(4) * t102 + t101 * t196) * t207 / 0.2e1) * qJD(4), -qJD(4) * t71 - t217; t7, qJD(3) * t10 - qJD(4) * t5 - qJD(5) * t8, t157 + t2 * qJD(4) + t3 * qJD(5) + qJD(3) * t209 + (qJD(1) * t226 + (t146 * t82 + t79 * t83) * qJD(3) / 0.2e1) * t207 + (t218 * t124 + m(5) * (t156 * pkin(7) * t124 - t197) + (m(6) * t115 - mrSges(4,1) + t107 + t143) * t121) * qJD(3) * pkin(2), t2 * qJD(3) + (t107 * t114 + t141 + t234) * qJD(4) + t230 + (-t153 + m(6) * (t119 * t148 - t122 * t68)) * t164 - t133, t3 * qJD(3) + t11 * qJD(4) - t132 + t230; -t154 / 0.2e1, -qJD(1) * t235 / 0.2e1 - t1 * qJD(4) + t4 * qJD(5) - t157, -qJD(4) * t6 - qJD(5) * t9, (t107 * pkin(7) + t141 + t233) * qJD(4) + t229 + (-t153 + m(6) * (t119 * t146 - t122 * t79)) * t164 - t129, t15 * qJD(4) - t128 + t229; 0, t1 * qJD(3) + t13 * qJD(5) + t133, t16 * qJD(5) + t129, -t100, -t100 - t126; 0, -t4 * qJD(3) - t13 * qJD(4) + t132, -t16 * qJD(4) + t128, t126, 0;];
Cq = t12;
