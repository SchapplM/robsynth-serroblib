% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:39
% EndTime: 2019-12-31 17:47:40
% DurationCPUTime: 0.81s
% Computational Cost: add. (1787->149), mult. (4153->236), div. (0->0), fcn. (3744->6), ass. (0->103)
t90 = cos(pkin(8));
t166 = m(6) * t90;
t165 = pkin(3) + qJ(2);
t93 = cos(qJ(5));
t152 = -t93 / 0.2e1;
t92 = sin(qJ(5));
t153 = t92 / 0.2e1;
t91 = cos(pkin(7));
t142 = t90 * t91;
t138 = t92 * t91;
t89 = sin(pkin(7));
t144 = t89 * t93;
t88 = sin(pkin(8));
t65 = t88 * t138 + t144;
t149 = t65 * mrSges(6,3);
t44 = -mrSges(6,2) * t142 + t149;
t154 = -t44 / 0.2e1;
t141 = t91 * t93;
t145 = t89 * t92;
t66 = t88 * t141 - t145;
t148 = t66 * mrSges(6,3);
t45 = mrSges(6,1) * t142 + t148;
t159 = t45 * t152 + t92 * t154 + (t93 * t66 / 0.2e1 + t65 * t153) * mrSges(6,3);
t164 = t159 * t88;
t161 = -t92 ^ 2 - t93 ^ 2;
t163 = t88 * (0.1e1 + t161) * t166;
t162 = mrSges(4,1) + mrSges(3,3);
t114 = -qJ(3) * t89 - pkin(1);
t68 = (-pkin(2) - qJ(4)) * t91 + t114;
t76 = t165 * t89;
t111 = -t68 * t88 + t76 * t90;
t135 = t90 * t68 + t88 * t76;
t136 = t93 * t44;
t139 = t92 * t45;
t30 = pkin(6) * t89 + t135;
t133 = t165 * t91;
t46 = (pkin(4) * t90 + pkin(6) * t88) * t91 + t133;
t22 = -t92 * t30 + t46 * t93;
t23 = t30 * t93 + t92 * t46;
t29 = -pkin(4) * t89 - t111;
t158 = (-m(5) * t111 + m(6) * t29) * t88 + (t136 - t139 + m(5) * t135 - m(6) * (t22 * t92 - t23 * t93)) * t90;
t85 = t91 ^ 2;
t157 = 0.2e1 * t89;
t156 = m(5) / 0.2e1;
t31 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t155 = t31 / 0.2e1;
t151 = m(4) * t89;
t150 = t91 * t163;
t147 = t88 * t89;
t146 = t88 * t91;
t143 = t90 * t89;
t140 = t92 * mrSges(6,1);
t137 = t93 * mrSges(6,2);
t134 = Ifges(6,5) * t65 + Ifges(6,6) * t66;
t82 = t88 ^ 2;
t84 = t90 ^ 2;
t132 = t82 + t84;
t120 = t142 / 0.2e1;
t1 = t134 * t120 + t29 * t31 + t22 * t44 - t23 * t45 + (-t22 * mrSges(6,3) + Ifges(6,4) * t65 + Ifges(6,5) * t120) * t65 + (Ifges(6,6) * t120 + t23 * mrSges(6,3) - Ifges(6,4) * t66 + (Ifges(6,2) - Ifges(6,1)) * t65) * t66;
t131 = t1 * qJD(1);
t110 = mrSges(5,1) * t89 + mrSges(5,3) * t146;
t100 = t89 * t110;
t32 = -mrSges(6,1) * t65 - mrSges(6,2) * t66;
t64 = -t88 * t145 + t141;
t67 = t88 * t144 + t138;
t70 = -mrSges(5,2) * t89 - mrSges(5,3) * t142;
t83 = t89 ^ 2;
t2 = -t32 * t143 + m(6) * (-t29 * t143 + t22 * t64 + t23 * t67) + t67 * t44 + t64 * t45 + t90 * t100 + m(5) * (t133 * t91 + (t111 * t90 + t135 * t88) * t89) + t70 * t147 + t162 * t83 + (m(4) + m(3)) * (t83 + t85) * qJ(2) + (t90 * mrSges(5,1) - mrSges(5,2) * t88 + t162) * t85;
t130 = t2 * qJD(1);
t3 = t32 * t147 - t88 * t100 + t70 * t143 + (-pkin(2) * t91 + t114) * t151 + (t91 * mrSges(4,2) - t89 * mrSges(4,3) + t158) * t89;
t129 = t3 * qJD(1);
t109 = t64 * mrSges(6,1) / 0.2e1 - t67 * mrSges(6,2) / 0.2e1;
t94 = t88 * t155 + t159 * t90;
t5 = t94 - t109;
t128 = t5 * qJD(1);
t107 = t137 / 0.2e1 + t140 / 0.2e1;
t7 = (t107 * t89 + t155) * t90 - t164;
t127 = t7 * qJD(1);
t9 = (mrSges(6,2) * t120 + t154 + t149 / 0.2e1) * t93 + (mrSges(6,1) * t120 - t148 / 0.2e1 + t45 / 0.2e1) * t92;
t126 = t9 * qJD(1);
t112 = -t64 * t92 + t67 * t93;
t104 = t112 * t88;
t101 = m(6) * (t161 * t84 - t82);
t98 = t101 / 0.4e1;
t16 = -t151 - m(6) * t104 / 0.2e1 + (t98 - m(6) * t84 / 0.4e1 + (-t82 / 0.2e1 - t84 / 0.2e1) * m(5)) * t157;
t125 = t16 * qJD(1);
t105 = m(6) * (t64 * t93 + t92 * t67);
t115 = m(5) * t132;
t19 = -t105 / 0.2e1 + (-m(5) / 0.2e1 + t101 / 0.2e1 - t115 / 0.2e1) * t91;
t124 = t19 * qJD(1);
t121 = t150 / 0.2e1;
t103 = qJD(5) * (-mrSges(6,1) * t93 + mrSges(6,2) * t92);
t4 = t70 * t142 + (-t110 + t32) * t146 + t158 * t91;
t102 = -t4 * qJD(1) + qJD(3) * t121;
t99 = t107 * t90;
t95 = t98 - t115 / 0.4e1;
t28 = qJD(4) * t121;
t18 = t105 / 0.2e1 + (t156 + 0.2e1 * t95) * t91;
t17 = m(6) * (t84 * t89 + t104) / 0.2e1 + t132 * t89 * t156 + t95 * t157;
t10 = t136 / 0.2e1 + t148 * t153 - t139 / 0.2e1 + t149 * t152 + t91 * t99;
t8 = -t90 * t31 / 0.2e1 + t89 * t99 + t164;
t6 = t94 + t109;
t11 = [qJD(2) * t2 - qJD(3) * t3 - qJD(4) * t4 + qJD(5) * t1, t130 + t17 * qJD(3) + t18 * qJD(4) + t6 * qJD(5) + (t112 - t147) * qJD(2) * t166, t89 * qJD(3) * t163 + t17 * qJD(2) + t8 * qJD(5) - t129 + t28, t18 * qJD(2) + t10 * qJD(5) + t102, t131 + t6 * qJD(2) + t8 * qJD(3) + t10 * qJD(4) + (-mrSges(6,1) * t23 - mrSges(6,2) * t22 + t134) * qJD(5); qJD(3) * t16 + qJD(4) * t19 + qJD(5) * t5 - t130, 0, t125, t124, t90 * t103 + t128; -qJD(2) * t16 - qJD(5) * t7 + t129 + t28, -t125, 0, qJD(1) * t121, t88 * t103 - t127; -t19 * qJD(2) - t9 * qJD(5) - t102, -t124, -qJD(1) * t150 / 0.2e1, 0, -t126 + (-t137 - t140) * qJD(5); -qJD(2) * t5 + qJD(3) * t7 + qJD(4) * t9 - t131, -t128, t127, t126, 0;];
Cq = t11;
