% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:56
% EndTime: 2019-12-31 18:25:57
% DurationCPUTime: 0.65s
% Computational Cost: add. (1758->112), mult. (2874->150), div. (0->0), fcn. (2471->6), ass. (0->72)
t81 = cos(qJ(5));
t78 = t81 ^ 2;
t80 = sin(qJ(5));
t110 = t80 ^ 2 + t78;
t92 = -Ifges(6,4) * t80 + (Ifges(6,1) - Ifges(6,2)) * t81;
t135 = t92 * t80;
t100 = t81 * mrSges(6,1) - t80 * mrSges(6,2);
t106 = sin(pkin(8));
t107 = cos(pkin(8));
t123 = sin(qJ(3));
t124 = cos(qJ(3));
t64 = -t106 * t124 - t107 * t123;
t40 = t64 * t100;
t133 = -t78 * Ifges(6,4) - t135;
t132 = t110 * mrSges(6,3);
t63 = t106 * t123 - t107 * t124;
t85 = -t123 * mrSges(4,1) + t64 * mrSges(5,1) - t124 * mrSges(4,2) + (mrSges(5,2) - t132) * t63;
t128 = t85 + t40;
t82 = -pkin(1) - pkin(2);
t68 = t124 * qJ(2) + t123 * t82;
t101 = t106 * t68;
t67 = -t123 * qJ(2) + t124 * t82;
t35 = t107 * t67 - t101;
t102 = t110 * t35;
t10 = m(6) * (-0.1e1 + t110) * t64 * t63;
t127 = t10 * qJD(2);
t126 = m(5) * pkin(3);
t125 = -t35 / 0.2e1;
t121 = Ifges(6,4) * t81;
t66 = -pkin(3) + t67;
t89 = t107 * t66 - t101;
t31 = pkin(4) - t89;
t114 = t81 * mrSges(6,2);
t115 = t80 * mrSges(6,1);
t70 = t114 + t115;
t120 = t31 * t70;
t119 = t63 * t70;
t75 = -t107 * pkin(3) - pkin(4);
t118 = t75 * t70;
t49 = t107 * t68;
t111 = t106 * t66 + t49;
t32 = -pkin(7) + t111;
t34 = t106 * t67 + t49;
t87 = -t68 * mrSges(4,1) - t35 * mrSges(5,2) - t67 * mrSges(4,2) + (-mrSges(5,1) - t100) * t34;
t3 = mrSges(6,3) * t102 - m(6) * (t32 * t102 + t31 * t34) - m(5) * (t111 * t35 - t89 * t34) + t87;
t109 = t3 * qJD(1);
t103 = t110 * t63;
t4 = -mrSges(3,3) - m(3) * qJ(2) - m(6) * (-t32 * t103 - t31 * t64) - m(5) * (-t111 * t63 + t89 * t64) - m(4) * (-t67 * t123 + t68 * t124) + t128;
t108 = t4 * qJD(1);
t13 = -t120 - t133;
t105 = t13 * qJD(1);
t93 = t114 / 0.2e1 + t115 / 0.2e1;
t19 = (t70 / 0.2e1 + t93) * t63;
t104 = t19 * qJD(1);
t95 = -Ifges(6,5) * t81 + Ifges(6,6) * t80;
t83 = -t40 / 0.2e1 + m(5) * ((-t35 + t89) * t64 + (t34 - t111) * t63) / 0.2e1 + m(6) * ((-t31 - t102) * t64 + (-t110 * t32 + t34) * t63) / 0.2e1;
t86 = (-t106 * t63 + t107 * t64) * t126;
t74 = t106 * pkin(3) + pkin(7);
t88 = m(6) * (-t74 * t103 - t75 * t64);
t84 = t88 / 0.2e1 + t40 / 0.2e1 + t86 / 0.2e1;
t1 = t83 - t84 - t85;
t94 = t1 * qJD(1) + t127;
t91 = t93 * t63;
t17 = -t118 + t133;
t18 = (-t70 / 0.2e1 + t93) * t63;
t6 = (mrSges(6,2) * t125 + t121) * t81 + (t75 / 0.2e1 - t31 / 0.2e1) * t70 + (mrSges(6,1) * t125 + t92) * t80;
t90 = t6 * qJD(1) + t18 * qJD(2) + t17 * qJD(3);
t21 = -t119 / 0.2e1 + t91;
t20 = t119 / 0.2e1 + t91;
t7 = -t118 / 0.2e1 + t120 / 0.2e1 - t93 * t35 - t121 * t81 - t135;
t2 = t83 + t84;
t5 = [-t4 * qJD(2) - t3 * qJD(3) + t13 * qJD(5), t2 * qJD(3) + t21 * qJD(5) - t108 + t127, -t109 + t2 * qJD(2) + (m(6) * (t74 * t102 + t75 * t34) - t107 * t34 * t126 + t87 + (t106 * t126 + t132) * t35) * qJD(3) + t7 * qJD(5), 0, t105 + t21 * qJD(2) + t7 * qJD(3) + (-t100 * t32 + t95) * qJD(5); t1 * qJD(3) - t19 * qJD(5) + t108, t10 * qJD(3), (t86 + t88 + t128) * qJD(3) + t20 * qJD(5) + t94, 0, t20 * qJD(3) + qJD(5) * t40 - t104; -qJD(2) * t1 - qJD(5) * t6 + t109, -qJD(5) * t18 - t94, -t17 * qJD(5), 0, (-t100 * t74 - t95) * qJD(5) - t90; 0, 0, 0, 0, -t70 * qJD(5); qJD(2) * t19 + qJD(3) * t6 - t105, qJD(3) * t18 + t104, t90, 0, 0;];
Cq = t5;
