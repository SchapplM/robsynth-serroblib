% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:26
% EndTime: 2019-12-05 15:16:29
% DurationCPUTime: 1.20s
% Computational Cost: add. (2664->143), mult. (7285->220), div. (0->0), fcn. (7523->8), ass. (0->92)
t119 = sin(qJ(5));
t122 = cos(qJ(5));
t118 = sin(pkin(9));
t123 = cos(qJ(4));
t120 = sin(qJ(4));
t124 = cos(qJ(3));
t160 = t120 * t124;
t163 = cos(pkin(9));
t90 = -t118 * t160 - t163 * t123;
t157 = t123 * t124;
t91 = t118 * t157 - t163 * t120;
t143 = -t119 * t91 + t122 * t90;
t57 = t119 * t90 + t122 * t91;
t17 = -t57 * mrSges(6,1) - t143 * mrSges(6,2);
t209 = t17 * qJD(5);
t158 = t122 * t123;
t135 = t119 * t120 - t158;
t136 = t119 * t123 + t122 * t120;
t190 = -pkin(7) - pkin(6);
t105 = t190 * t120;
t106 = t190 * t123;
t142 = t122 * t105 + t119 * t106;
t73 = t119 * t105 - t122 * t106;
t18 = -t73 * mrSges(6,1) - t142 * mrSges(6,2) - Ifges(6,5) * t135 - Ifges(6,6) * t136;
t208 = t18 * qJD(5);
t121 = sin(qJ(3));
t161 = t120 * t121;
t85 = t119 * t161 - t121 * t158;
t86 = t136 * t121;
t144 = t85 * mrSges(6,1) + t86 * mrSges(6,2);
t207 = t144 * qJD(5);
t191 = m(6) * pkin(4);
t206 = -mrSges(6,1) / 0.2e1;
t205 = mrSges(6,2) / 0.2e1;
t114 = -t123 * pkin(4) - pkin(3);
t204 = m(6) * t114;
t116 = t120 ^ 2;
t154 = t123 ^ 2 + t116;
t203 = pkin(6) * t154;
t62 = mrSges(6,1) * t136 - mrSges(6,2) * t135;
t9 = -Ifges(6,4) * t136 ^ 2 - (-Ifges(6,4) * t135 - (-Ifges(6,1) + Ifges(6,2)) * t136) * t135 + t114 * t62;
t199 = m(5) * (-0.1e1 + t154);
t164 = t123 * mrSges(5,2);
t167 = t120 * mrSges(5,1);
t104 = t164 + t167;
t103 = -t123 * mrSges(5,1) + t120 * mrSges(5,2);
t87 = t136 * t124;
t88 = t135 * t124;
t170 = t88 * t205 + t87 * t206;
t162 = t118 * t121;
t74 = t136 * t162;
t75 = t135 * t162;
t134 = t75 * t205 + t74 * t206;
t195 = t154 * mrSges(5,3) - mrSges(4,2);
t193 = m(5) / 0.2e1;
t192 = m(6) / 0.2e1;
t188 = -t124 / 0.2e1;
t159 = t121 * t124;
t156 = t124 ^ 2 * t118;
t151 = t62 * t188;
t63 = mrSges(6,1) * t135 + mrSges(6,2) * t136;
t148 = -mrSges(4,1) + t103 + t63;
t145 = t162 / 0.2e1;
t141 = t62 * t145;
t12 = (-t143 * t87 - t88 * t57 - t86 * t74 - t85 * t75 - t156) * t192 + (t91 * t157 - t90 * t160 - t156) * t193 + 0.2e1 * (m(6) / 0.4e1 - t199 / 0.4e1) * t118 * t121 ^ 2;
t107 = t118 ^ 2 * t159;
t24 = m(6) * (t143 * t74 + t57 * t75 + t107) + m(5) * (t107 + (t120 * t90 - t123 * t91) * t162);
t139 = t24 * qJD(1) + t12 * qJD(2);
t28 = t159 * t199 + (t85 * t88 + t86 * t87 - t159) * m(6);
t138 = t12 * qJD(1) + t28 * qJD(2);
t6 = t141 - t134;
t133 = -t164 / 0.2e1 - t167 / 0.2e1;
t132 = (t119 * t75 + t122 * t74) * t191;
t125 = t118 * pkin(4) * t161 * t192;
t129 = t104 / 0.2e1 + t62 / 0.2e1 + t133;
t3 = t129 * t162 - t132 / 0.2e1 + t125 + t134;
t4 = -pkin(3) * t104 - t116 * Ifges(5,4) + (Ifges(5,4) * t123 + (Ifges(5,1) - Ifges(5,2)) * t120) * t123 + t9 + (t204 + t63) * pkin(4) * t120;
t126 = t160 * t191;
t128 = (-t119 * t88 - t122 * t87) * t191 / 0.2e1 + t170;
t7 = t129 * t124 + t126 / 0.2e1 + t128;
t131 = t3 * qJD(1) - t7 * qJD(2) + t4 * qJD(3);
t10 = t151 - t170;
t5 = t141 + t134;
t130 = t5 * qJD(1) + t10 * qJD(2) + t9 * qJD(3);
t102 = (mrSges(6,1) * t119 + mrSges(6,2) * t122) * pkin(4);
t127 = t102 * qJD(4);
t96 = t102 * qJD(5);
t11 = t151 + t170;
t8 = t104 * t188 - t126 / 0.2e1 + t133 * t124 + t128 + t151;
t2 = t132 / 0.2e1 + t125 + t6 + 0.2e1 * t104 * t145;
t1 = t12 * qJD(3);
t13 = [t24 * qJD(3), t1, t2 * qJD(4) + t6 * qJD(5) + t139 + (m(6) * (t142 * t74 + t73 * t75) + ((-m(5) * pkin(3) + t148 + t204) * t124 + (-m(5) * t203 - t195) * t121) * t118 + (-t75 * t135 - t74 * t136) * mrSges(6,3)) * qJD(3), t2 * qJD(3) + (-t91 * mrSges(5,1) - t90 * mrSges(5,2) + (t119 * t143 - t122 * t57) * t191 + t17) * qJD(4) + t209, t6 * qJD(3) + t17 * qJD(4) + t209; t1, t28 * qJD(3), t8 * qJD(4) + t11 * qJD(5) + t138 + ((t135 * t88 + t136 * t87) * mrSges(6,3) + t195 * t124 + t148 * t121 + 0.2e1 * (-pkin(3) * t121 + t124 * t203) * t193 + 0.2e1 * (t114 * t121 - t142 * t87 - t73 * t88) * t192) * qJD(3), t8 * qJD(3) + (t121 * t103 + t144 + (-t119 * t86 + t122 * t85) * t191) * qJD(4) + t207, t11 * qJD(3) + qJD(4) * t144 + t207; t3 * qJD(4) + t5 * qJD(5) - t139, -t7 * qJD(4) + t10 * qJD(5) - t138, t4 * qJD(4) + t9 * qJD(5), t208 + t131 + (Ifges(5,5) * t123 - Ifges(5,6) * t120 + (m(6) * (t119 * t142 - t122 * t73) + (-t119 * t136 + t122 * t135) * mrSges(6,3)) * pkin(4) + t103 * pkin(6) + t18) * qJD(4), t18 * qJD(4) + t130 + t208; -t3 * qJD(3), t7 * qJD(3), -t131, -t96, -t96 - t127; -t5 * qJD(3), -t10 * qJD(3), -t130, t127, 0;];
Cq = t13;
