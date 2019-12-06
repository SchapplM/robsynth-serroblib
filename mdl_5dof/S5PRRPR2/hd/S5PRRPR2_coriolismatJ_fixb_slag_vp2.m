% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:23
% EndTime: 2019-12-05 16:17:24
% DurationCPUTime: 0.80s
% Computational Cost: add. (1936->157), mult. (4556->228), div. (0->0), fcn. (3450->6), ass. (0->112)
t166 = qJD(2) + qJD(3);
t96 = cos(qJ(5));
t135 = t96 * mrSges(6,2);
t94 = sin(qJ(5));
t138 = t94 * mrSges(6,1);
t119 = t135 + t138;
t92 = sin(pkin(9));
t90 = t92 ^ 2;
t93 = cos(pkin(9));
t91 = t93 ^ 2;
t165 = t90 * t119 + (t90 + t91) * mrSges(5,3);
t97 = cos(qJ(3));
t147 = pkin(2) * t97;
t95 = sin(qJ(3));
t146 = t95 * pkin(2);
t164 = -mrSges(5,1) * t93 + mrSges(5,2) * t92 - mrSges(4,1);
t120 = mrSges(6,1) * t96 - mrSges(6,2) * t94;
t150 = t96 / 0.2e1;
t154 = t96 ^ 2;
t155 = t94 ^ 2;
t136 = t94 * t92;
t126 = mrSges(6,3) * t136;
t75 = t93 * mrSges(6,2) - t126;
t133 = t96 * t92;
t125 = mrSges(6,3) * t133;
t76 = -t93 * mrSges(6,1) - t125;
t99 = (t154 / 0.2e1 + t155 / 0.2e1) * t90 * mrSges(6,3) + (t93 * t120 / 0.2e1 + t94 * t75 / 0.2e1 + t76 * t150) * t92;
t163 = qJD(1) * t99;
t86 = qJ(4) + t146;
t79 = t90 * t86;
t87 = t90 * qJ(4);
t131 = t79 + t87;
t140 = t86 * t93;
t77 = -t93 * pkin(4) - t92 * pkin(7) - pkin(3);
t73 = t77 - t147;
t44 = -t94 * t140 + t73 * t96;
t45 = t96 * t140 + t73 * t94;
t129 = qJ(4) * t93;
t56 = -t94 * t129 + t77 * t96;
t57 = t96 * t129 + t77 * t94;
t162 = ((t45 + t57) * t96 + (-t44 - t56) * t94) * t93 + t131;
t107 = (-t135 / 0.2e1 - t138 / 0.2e1) * t93;
t134 = t96 * t75;
t137 = t94 * t76;
t158 = -t134 / 0.2e1 + t137 / 0.2e1;
t101 = t107 + t158;
t161 = qJD(4) * t101;
t104 = t107 - t158;
t160 = t104 * qJD(4);
t159 = -mrSges(4,2) + t165;
t52 = (t154 - 0.1e1 + t155) * t93 * t92;
t148 = m(6) * t52;
t123 = t148 / 0.2e1;
t132 = qJD(4) * t123 - qJD(5) * t99;
t152 = -m(6) / 0.2e1;
t122 = qJD(1) * t152;
t157 = -qJD(5) * t101 + t52 * t122;
t156 = qJD(1) * t123 + qJD(5) * t104;
t153 = m(5) / 0.2e1;
t151 = m(6) / 0.2e1;
t139 = t93 * t97;
t70 = (-t94 * t139 + t95 * t96) * pkin(2);
t71 = (t96 * t139 + t94 * t95) * pkin(2);
t26 = (-pkin(2) * t139 - t70 * t94 + t71 * t96) * t92;
t149 = m(6) * t26;
t145 = Ifges(6,4) * t94;
t144 = Ifges(6,4) * t96;
t143 = Ifges(6,5) * t93;
t142 = Ifges(6,6) * t93;
t141 = t86 * t91;
t130 = qJ(4) * t91;
t109 = (t134 - t137) * t93 + t165;
t10 = m(6) * (t79 + (-t44 * t94 + t45 * t96) * t93) + m(5) * (t79 + t141) + t109;
t128 = qJD(2) * t10;
t127 = t90 * t147;
t124 = t149 / 0.2e1;
t118 = t70 * t76 + t71 * t75;
t103 = (t144 + (Ifges(6,1) - Ifges(6,2)) * t94) * t92 - t142;
t108 = (-Ifges(6,4) * t136 - t143) * t94;
t31 = t44 * t75;
t32 = t45 * t76;
t35 = t44 * t126;
t112 = t90 * t120;
t53 = t86 * t112;
t3 = -t53 - t31 + t32 - t35 + (t108 + (t45 * mrSges(6,3) + t103) * t96) * t92;
t117 = -t3 * qJD(2) - t163;
t13 = m(6) * (t87 + (-t56 * t94 + t57 * t96) * t93) + m(5) * (t87 + t130) + t109;
t100 = ((qJ(4) + t86) * t91 + t131) * t153 + t109;
t102 = (t96 * t70 + t94 * t71) * t151 + t146 * t153;
t5 = t162 * t152 - t100 + t102;
t116 = qJD(2) * t5 - qJD(3) * t13;
t115 = t166 * t99;
t114 = t166 * t101;
t113 = t70 * mrSges(6,1) / 0.2e1 - t71 * mrSges(6,2) / 0.2e1;
t111 = (-Ifges(6,5) * t94 - Ifges(6,6) * t96) * t92;
t74 = t86 * t127;
t7 = (t159 * t97 + t164 * t95) * pkin(2) + m(6) * (t44 * t70 + t45 * t71 + t74) + m(5) * (t141 * t147 + t74 + (-pkin(3) - t147) * t146) + t118;
t110 = -t7 * qJD(2) + t26 * t122;
t38 = t56 * t75;
t39 = t57 * t76;
t48 = t56 * t126;
t65 = qJ(4) * t112;
t98 = -(-t143 + (Ifges(6,1) * t96 - t145) * t92) * t136 / 0.2e1 - (-t142 + (-Ifges(6,2) * t94 + t144) * t92) * t133 / 0.2e1 - t93 * t111 / 0.2e1 + (-t45 / 0.2e1 - t57 / 0.2e1) * t125 + t31 / 0.2e1 - t32 / 0.2e1 + t35 / 0.2e1 + t38 / 0.2e1 - t39 / 0.2e1 + t48 / 0.2e1 + t53 / 0.2e1 + t65 / 0.2e1 + (-t94 * (-Ifges(6,2) * t96 - t145) / 0.2e1 + (-Ifges(6,1) * t94 - t144) * t150) * t90;
t1 = t98 - t113;
t4 = -t65 - t38 + t39 - t48 + (t108 + (t57 * mrSges(6,3) + t103) * t96) * t92;
t106 = t1 * qJD(2) - t4 * qJD(3) - t163;
t84 = qJ(4) * t127;
t24 = 0.2e1 * (qJD(2) / 0.4e1 + qJD(3) / 0.4e1) * t148;
t8 = qJD(3) * t124 + t132;
t6 = t162 * t151 + t100 + t102;
t2 = t98 + t113;
t9 = [0, t8, qJD(2) * t124 + t132, t24, -t120 * qJD(5) * t92 - t115; t8, qJD(3) * t7 + qJD(4) * t10 - qJD(5) * t3, t6 * qJD(4) + t2 * qJD(5) - t110 + (t118 + 0.2e1 * (t56 * t70 + t57 * t71 + t84) * t151 + 0.2e1 * t84 * t153 + (-m(5) * pkin(3) + t164) * t146 + (m(5) * t130 + t159) * t147) * qJD(3), qJD(3) * t6 + t128 + t156, t2 * qJD(3) + t160 + (-mrSges(6,1) * t45 - mrSges(6,2) * t44 + t111) * qJD(5) + t117; -qJD(2) * t149 / 0.2e1 + t132, -t5 * qJD(4) + t1 * qJD(5) + t110, qJD(4) * t13 - qJD(5) * t4, -t116 + t156, t160 + (-mrSges(6,1) * t57 - mrSges(6,2) * t56 + t111) * qJD(5) + t106; -t24, qJD(3) * t5 - t128 + t157, t116 + t157, 0, -t119 * qJD(5) - t114; t115, -qJD(3) * t1 - t117 + t161, -t106 + t161, t114, 0;];
Cq = t9;
