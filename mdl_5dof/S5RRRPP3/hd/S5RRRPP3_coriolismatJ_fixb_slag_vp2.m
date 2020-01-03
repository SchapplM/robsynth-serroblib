% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:10
% EndTime: 2019-12-31 20:53:12
% DurationCPUTime: 0.90s
% Computational Cost: add. (1660->175), mult. (3175->202), div. (0->0), fcn. (2069->4), ass. (0->101)
t96 = sin(qJ(3));
t98 = cos(qJ(3));
t165 = t96 ^ 2 + t98 ^ 2;
t156 = m(5) / 0.2e1;
t99 = cos(qJ(2));
t152 = t99 * pkin(1);
t135 = t96 * qJ(4);
t115 = -t98 * pkin(3) - t135;
t57 = -pkin(2) + t115;
t38 = t57 - t152;
t59 = t98 * mrSges(5,2) - t96 * mrSges(5,3);
t164 = (t38 + t57) * t156 + t59;
t163 = t165 * t152;
t162 = Ifges(6,6) - Ifges(5,6) - Ifges(4,4);
t139 = qJ(4) * t98;
t61 = t96 * pkin(3) - t139;
t161 = m(5) * t61;
t160 = mrSges(5,3) + mrSges(6,2);
t60 = -t98 * mrSges(4,1) + t96 * mrSges(4,2);
t159 = t60 - mrSges(3,1);
t104 = (Ifges(4,1) - Ifges(4,2) + Ifges(5,2) - Ifges(6,2) - Ifges(5,3) + Ifges(6,3)) * t98 + t162 * t96;
t116 = t162 * t98;
t158 = t104 * t96 - t116 * t98;
t62 = -mrSges(6,2) * t96 - mrSges(6,3) * t98;
t157 = t164 + t62;
t155 = m(6) / 0.2e1;
t95 = pkin(3) + qJ(5);
t154 = m(6) * t95;
t153 = pkin(7) * t98;
t151 = mrSges(4,2) * t98;
t150 = mrSges(5,2) * t96;
t51 = qJ(5) * t96 + t61;
t149 = t51 * t62;
t97 = sin(qJ(2));
t81 = pkin(1) * t97 + pkin(7);
t148 = t81 * t98;
t147 = t96 * mrSges(6,1);
t84 = t98 * mrSges(5,1);
t83 = t98 * mrSges(6,1);
t146 = -mrSges(5,1) - mrSges(6,1);
t141 = t83 + t84;
t140 = m(6) * qJD(3);
t35 = -t95 * t98 - pkin(2) - t135;
t32 = t35 - t152;
t122 = m(6) * t32 + t62;
t109 = m(5) * t38 + t122 + t59;
t9 = t109 * t96;
t138 = qJD(1) * t9;
t101 = t61 * t59 + t149 + t158;
t63 = -mrSges(5,3) * t98 - t150;
t123 = t63 + t161;
t19 = t32 * t51;
t58 = -mrSges(6,2) * t98 + mrSges(6,3) * t96;
t64 = mrSges(4,1) * t96 + t151;
t82 = -pkin(2) - t152;
t3 = m(6) * t19 + t123 * t38 + t32 * t58 + t82 * t64 + t101;
t137 = t3 * qJD(1);
t105 = -mrSges(3,2) + t165 * (mrSges(4,3) - t146);
t52 = (pkin(4) + t81) * t96;
t92 = t98 * pkin(4);
t53 = t92 + t148;
t5 = ((m(4) * t82 + t109 + t159) * t97 + (m(6) * (t52 * t96 + t53 * t98) + t105) * t99) * pkin(1) + (m(5) + m(4)) * t163 * t81;
t136 = t5 * qJD(1);
t18 = t122 * t98;
t134 = qJD(1) * t18;
t76 = mrSges(6,3) + t154;
t133 = t76 * qJD(3);
t77 = (m(5) + m(6)) * qJ(4) + t160;
t132 = t77 * qJD(3);
t127 = t152 / 0.2e1;
t126 = t35 / 0.2e1 + t32 / 0.2e1;
t124 = m(6) * (-t32 - t35);
t121 = m(6) * t35 + t62;
t120 = m(6) * t127;
t117 = t124 / 0.2e1;
t114 = (t155 + t156) * t152;
t22 = t35 * t51;
t100 = t149 + t126 * t58 + (t57 / 0.2e1 + t38 / 0.2e1) * t63 + t164 * t61 + (t82 / 0.2e1 - pkin(2) / 0.2e1) * t64 + (t19 + t22) * t155;
t2 = ((mrSges(4,2) / 0.2e1 - mrSges(6,2) / 0.2e1 - mrSges(5,3) / 0.2e1 + 0.2e1 * (-m(5) / 0.4e1 - m(6) / 0.4e1) * qJ(4)) * t152 - t116) * t98 + ((mrSges(4,1) / 0.2e1 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1 + pkin(3) * t156 + t154 / 0.2e1) * t152 + t104) * t96 + t100;
t4 = m(6) * t22 - pkin(2) * t64 + t123 * t57 + t35 * t58 + t101;
t113 = -t2 * qJD(1) - t4 * qJD(2);
t108 = m(5) * t57 + t121 + t59;
t12 = t108 * t96;
t7 = (t114 - t124 / 0.2e1 + t157) * t96;
t112 = qJD(1) * t7 + qJD(2) * t12;
t10 = (t62 + (t127 + t126) * m(6)) * t98;
t20 = t121 * t98;
t111 = qJD(1) * t10 + qJD(2) * t20;
t107 = (t146 * qJ(4) - Ifges(4,6)) * qJD(3) * t96;
t106 = -pkin(3) * t84 - t95 * t83 + (Ifges(6,4) + Ifges(5,5)) * t96 + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t98;
t102 = qJD(3) * (m(5) * t115 + t59 + t60);
t72 = t92 + t153;
t71 = (pkin(4) + pkin(7)) * t96;
t33 = -m(6) * t71 - t147;
t26 = -m(6) * t52 - t147;
t23 = m(5) * t153 + m(6) * t72 + t141;
t16 = m(5) * t148 + m(6) * t53 + t141;
t11 = (t117 - t62 + t120) * t98;
t8 = (t117 + t114 - t157) * t96;
t1 = (-t95 * t96 + t139) * t120 + t100 + t158 - (t151 + (mrSges(6,3) + mrSges(4,1)) * t96) * t152 / 0.2e1 + (t160 * t98 + t150 - t161) * t127;
t6 = [qJD(2) * t5 + qJD(3) * t3 - qJD(4) * t9 - qJD(5) * t18, t1 * qJD(3) + t8 * qJD(4) + t11 * qJD(5) + t136 + (((-m(4) * pkin(2) + t108 + t159) * t97 + (m(6) * (t71 * t96 + t72 * t98) + t105) * t99) * pkin(1) + 0.2e1 * (t156 + m(4) / 0.2e1) * t163 * pkin(7)) * qJD(2), t137 + t1 * qJD(2) + (m(6) * (-qJ(4) * t52 - t53 * t95) - t53 * mrSges(6,3) - t52 * mrSges(6,2) + t106) * qJD(3) + t16 * qJD(4) + t26 * qJD(5) + t107 + t81 * t102, qJD(2) * t8 + qJD(3) * t16 - t138, qJD(2) * t11 + qJD(3) * t26 - t134; qJD(3) * t2 - qJD(4) * t7 - qJD(5) * t10 - t136, qJD(3) * t4 - qJD(4) * t12 - qJD(5) * t20, (m(6) * (-qJ(4) * t71 - t72 * t95) - t71 * mrSges(6,2) - t72 * mrSges(6,3) + t106) * qJD(3) + t23 * qJD(4) + t33 * qJD(5) + t107 + pkin(7) * t102 - t113, qJD(3) * t23 - t112, qJD(3) * t33 - t111; -qJD(2) * t2 - t137, t113, qJD(4) * t77 + qJD(5) * t76, t132, t133; qJD(2) * t7 + t138, t112, -m(6) * qJD(5) - t132, 0, -t140; qJD(2) * t10 + t134, t111, m(6) * qJD(4) - t133, t140, 0;];
Cq = t6;
