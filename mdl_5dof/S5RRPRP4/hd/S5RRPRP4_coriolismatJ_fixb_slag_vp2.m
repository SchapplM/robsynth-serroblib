% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:38
% EndTime: 2019-12-31 19:52:40
% DurationCPUTime: 0.81s
% Computational Cost: add. (918->147), mult. (1899->175), div. (0->0), fcn. (1079->4), ass. (0->90)
t119 = Ifges(6,5) - Ifges(5,4);
t75 = sin(qJ(4));
t158 = t119 * t75;
t77 = cos(qJ(4));
t47 = t75 * mrSges(6,1) - t77 * mrSges(6,3);
t92 = pkin(4) * t77 + qJ(5) * t75;
t24 = t92 * t47;
t72 = t77 ^ 2;
t157 = t119 * t72 + t24;
t151 = t75 * mrSges(5,1) + t77 * mrSges(5,2) + t47;
t156 = mrSges(4,3) + t151;
t149 = t75 * pkin(4) - t77 * qJ(5);
t80 = qJD(4) * (-m(6) * t149 - t151);
t115 = t75 ^ 2 + t72;
t154 = (mrSges(6,2) + mrSges(5,3)) * t115 + mrSges(3,1);
t152 = -Ifges(5,1) - Ifges(6,1);
t150 = -Ifges(5,2) - t152;
t76 = sin(qJ(2));
t70 = t76 * pkin(1);
t102 = t70 / 0.2e1;
t95 = -qJ(3) - t149;
t31 = t70 - t95;
t96 = -(t31 - t95) * m(6) / 0.2e1;
t144 = m(5) / 0.4e1;
t148 = t144 + m(4) / 0.4e1;
t122 = t77 * mrSges(6,1);
t124 = t75 * mrSges(6,3);
t93 = t122 + t124;
t23 = t95 * t93;
t123 = t77 * mrSges(5,1);
t125 = t75 * mrSges(5,2);
t94 = t123 - t125;
t45 = qJ(3) * t94;
t147 = t23 / 0.2e1 - t45 / 0.2e1;
t142 = m(6) / 0.2e1;
t141 = m(6) / 0.4e1;
t140 = -t31 / 0.2e1;
t139 = t31 / 0.2e1;
t100 = t115 * t76;
t32 = pkin(1) * t100;
t138 = t32 / 0.2e1;
t103 = -t70 / 0.2e1;
t135 = m(6) * t31;
t134 = m(6) * t95;
t133 = m(6) * t75;
t78 = cos(qJ(2));
t132 = t78 * pkin(1);
t131 = m(6) * qJ(5);
t59 = t70 + qJ(3);
t126 = t59 * t78;
t120 = t78 * mrSges(3,2);
t79 = -pkin(2) - pkin(7);
t117 = t115 * t79 * t70;
t101 = -pkin(2) - t132;
t57 = -pkin(7) + t101;
t97 = t57 * t100;
t98 = mrSges(4,2) * t70 + t156 * t132;
t3 = -t98 + (t120 - m(6) * (t31 * t78 + t97) - m(5) * (t97 + t126) - m(4) * (t101 * t76 + t126)) * pkin(1) + t154 * t70;
t114 = t3 * qJD(1);
t83 = ((Ifges(6,3) + Ifges(5,2) / 0.2e1 - t150 / 0.2e1 + t152 / 0.2e1) * t77 - t158) * t75 + t157;
t84 = t59 * t94;
t86 = m(6) * t92;
t4 = t84 + (t86 + t93) * t31 + t83;
t113 = t4 * qJD(1);
t8 = 0.4e1 * t148 * t59 + t135 + t156;
t111 = t8 * qJD(1);
t16 = (t47 + t135) * t77;
t109 = t16 * qJD(1);
t64 = mrSges(6,3) + t131;
t108 = t64 * qJD(4);
t104 = t102 + qJ(3);
t82 = t158 + (-Ifges(6,3) + t150) * t77;
t1 = -t24 + (mrSges(6,1) * t140 - t59 * mrSges(5,1) / 0.2e1 - t119 * t77 + pkin(4) * t96 + (mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1 + pkin(4) * t142) * t70) * t77 + (mrSges(6,3) * t140 + t59 * mrSges(5,2) / 0.2e1 + qJ(5) * t96 + (mrSges(6,3) / 0.2e1 - mrSges(5,2) / 0.2e1 + t131 / 0.2e1) * t70 + t82) * t75 + t147;
t85 = t95 * t92;
t5 = m(6) * t85 + t82 * t75 - t157 + t23 - t45;
t91 = t1 * qJD(1) + t5 * qJD(2);
t14 = t134 - t156 + (-m(5) - m(4)) * qJ(3);
t6 = (t138 - t104) * m(5) + (t102 - t104) * m(4) + (t138 + t103 + t95) * m(6) - t156;
t89 = t6 * qJD(1) + t14 * qJD(2);
t17 = (t47 - t134) * t77;
t9 = (t47 + (t103 + t139 - t95 / 0.2e1) * m(6)) * t77;
t88 = t9 * qJD(1) + t17 * qJD(2);
t81 = ((-Ifges(6,4) - Ifges(5,5)) * t75 + (-Ifges(5,6) + Ifges(6,6)) * t77 + t149 * mrSges(6,2)) * qJD(4);
t60 = qJ(3) * t132;
t44 = (m(6) * t79 - mrSges(6,2)) * t75;
t29 = (m(6) * t57 - mrSges(6,2)) * t75;
t10 = (m(6) * t103 - t47 + t96) * t77;
t7 = 0.4e1 * t149 * t141 + m(4) * t102 + 0.2e1 * (t141 + t144) * t32 + t156 + (t141 + t148) * (0.2e1 * t70 + 0.4e1 * qJ(3));
t2 = (t123 / 0.2e1 + t122 / 0.2e1 + t124 / 0.2e1 - t125 / 0.2e1 + t86 / 0.2e1) * t70 + (t31 * t92 - t85) * t142 + t84 / 0.2e1 + t93 * t139 + t83 - t147;
t11 = [-t3 * qJD(2) + t8 * qJD(3) + t4 * qJD(4) - t16 * qJD(5), t7 * qJD(3) + t2 * qJD(4) + t10 * qJD(5) - t114 + (t98 + (-t154 * t76 - t120) * pkin(1) + 0.2e1 * (-t132 * t95 + t117) * t142 + m(5) * (t60 + t117) + m(4) * (-pkin(2) * t70 + t60)) * qJD(2), t7 * qJD(2) + t111, t2 * qJD(2) + t29 * qJD(5) + t57 * t80 + t113 + t81, t10 * qJD(2) + t29 * qJD(4) - t109; -t6 * qJD(3) - t1 * qJD(4) - t9 * qJD(5) + t114, -t14 * qJD(3) - t5 * qJD(4) - t17 * qJD(5), -t89, t44 * qJD(5) + t79 * t80 + t81 - t91, t44 * qJD(4) - t88; t6 * qJD(2) - t111, t89, 0, qJD(5) * t133 + t80, qJD(4) * t133; t1 * qJD(2) - t113, t91, 0, t64 * qJD(5), t108; t9 * qJD(2) + t109, t88, 0, -t108, 0;];
Cq = t11;
