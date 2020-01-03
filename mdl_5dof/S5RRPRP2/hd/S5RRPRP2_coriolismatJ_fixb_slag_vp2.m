% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:28
% EndTime: 2019-12-31 19:49:30
% DurationCPUTime: 0.53s
% Computational Cost: add. (999->126), mult. (2077->163), div. (0->0), fcn. (1500->6), ass. (0->85)
t73 = sin(qJ(4));
t75 = cos(qJ(4));
t125 = t73 ^ 2 + t75 ^ 2;
t76 = cos(qJ(2));
t113 = t76 * pkin(1);
t71 = sin(pkin(8));
t74 = sin(qJ(2));
t57 = t71 * t74 * pkin(1);
t72 = cos(pkin(8));
t33 = t72 * t113 - t57;
t124 = t125 * t33;
t96 = t73 * qJ(5);
t86 = -t75 * pkin(4) - t96;
t83 = -pkin(3) + t86;
t60 = pkin(2) + t113;
t90 = t72 * t60 - t57;
t18 = t83 - t90;
t114 = t72 * pkin(2);
t34 = t83 - t114;
t123 = t18 + t34;
t43 = -t75 * mrSges(6,1) - t73 * mrSges(6,3);
t122 = -t75 * mrSges(5,1) + t73 * mrSges(5,2) + t43;
t106 = t72 * t74;
t32 = (t71 * t76 + t106) * pkin(1);
t121 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t125) * t33 + (-t74 * mrSges(3,1) - t76 * mrSges(3,2)) * pkin(1) + (-mrSges(4,1) + t122) * t32;
t120 = t33 / 0.2e1;
t45 = -t73 * pkin(4) + t75 * qJ(5);
t116 = m(6) * t45;
t115 = m(6) * t75;
t112 = Ifges(5,4) * t73;
t111 = Ifges(6,5) * t75;
t28 = -pkin(3) - t90;
t103 = t75 * mrSges(5,2);
t105 = t73 * mrSges(5,1);
t47 = t103 + t105;
t110 = t28 * t47;
t59 = -pkin(3) - t114;
t109 = t59 * t47;
t104 = t73 * mrSges(6,1);
t62 = t75 * mrSges(6,2);
t102 = t75 * mrSges(6,3);
t89 = pkin(1) * t106 + t71 * t60;
t29 = pkin(7) + t89;
t101 = t124 * t29;
t58 = t71 * pkin(2) + pkin(7);
t100 = t124 * t58;
t1 = m(6) * (t18 * t32 + t101) + m(5) * (t28 * t32 + t101) + m(4) * (-t90 * t32 + t89 * t33) + t121;
t98 = t1 * qJD(1);
t27 = t45 * t43;
t64 = Ifges(6,5) * t73;
t48 = -Ifges(6,3) * t75 + t64;
t49 = Ifges(6,3) * t73 + t111;
t50 = Ifges(5,2) * t75 + t112;
t67 = Ifges(5,4) * t75;
t51 = -Ifges(5,2) * t73 + t67;
t52 = Ifges(6,1) * t73 - t111;
t53 = Ifges(6,1) * t75 + t64;
t54 = Ifges(5,1) * t73 + t67;
t55 = Ifges(5,1) * t75 - t112;
t81 = -t73 * t50 / 0.2e1 - t27 - t75 * t49 / 0.2e1 + (t48 + t53 + t55) * t73 / 0.2e1 + (t51 + t52 + t54) * t75 / 0.2e1;
t46 = -t102 + t104;
t91 = t46 - t116;
t4 = t91 * t18 + t110 + t81;
t97 = t4 * qJD(1);
t10 = (m(6) * t18 + t43) * t73;
t95 = t10 * qJD(1);
t61 = m(6) * qJ(5) + mrSges(6,3);
t94 = t61 * qJD(4);
t93 = -t18 / 0.2e1 - t34 / 0.2e1;
t92 = m(6) * t123;
t87 = -t92 / 0.2e1;
t77 = (-t103 / 0.2e1 - t105 / 0.2e1 - t104 / 0.2e1 + t102 / 0.2e1 + t116 / 0.2e1) * t33;
t2 = t27 + (-t28 / 0.2e1 - t59 / 0.2e1) * t47 + t93 * t46 + t45 * t92 / 0.2e1 + (t49 / 0.2e1 - t54 / 0.2e1 - t51 / 0.2e1 - t52 / 0.2e1) * t75 + (-t48 / 0.2e1 - t55 / 0.2e1 + t50 / 0.2e1 - t53 / 0.2e1) * t73 + t77;
t5 = t91 * t34 + t109 + t81;
t85 = t2 * qJD(1) - t5 * qJD(2);
t15 = (m(6) * t34 + t43) * t73;
t6 = (t43 + (t120 - t93) * m(6)) * t73;
t84 = t6 * qJD(1) + t15 * qJD(2);
t79 = (-mrSges(6,2) * t96 - pkin(4) * t62 + (Ifges(6,4) + Ifges(5,5)) * t75 + (-Ifges(5,6) + Ifges(6,6)) * t73) * qJD(4);
t78 = qJD(4) * (m(6) * t86 + t122);
t31 = t58 * t115 + t62;
t17 = t29 * t115 + t62;
t7 = (m(6) * t120 - t43 + t87) * t73;
t3 = t81 + t77 + t45 * t87 + t109 / 0.2e1 + t110 / 0.2e1 + t123 * t46 / 0.2e1;
t8 = [t1 * qJD(2) + t4 * qJD(4) - t10 * qJD(5), t3 * qJD(4) + t7 * qJD(5) + t98 + (m(6) * (t34 * t32 + t100) + m(5) * (t59 * t32 + t100) + m(4) * (-t32 * t72 + t33 * t71) * pkin(2) + t121) * qJD(2), 0, t3 * qJD(2) + t17 * qJD(5) + t29 * t78 + t79 + t97, t7 * qJD(2) + t17 * qJD(4) - t95; -t2 * qJD(4) - t6 * qJD(5) - t98, t5 * qJD(4) - t15 * qJD(5), 0, t31 * qJD(5) + t58 * t78 + t79 - t85, t31 * qJD(4) - t84; 0, 0, 0, (t116 + (-mrSges(5,2) + mrSges(6,3)) * t75) * qJD(4) + ((-mrSges(5,1) - mrSges(6,1)) * qJD(4) + m(6) * qJD(5)) * t73, m(6) * qJD(4) * t73; t2 * qJD(2) - t97, t85, 0, t61 * qJD(5), t94; t6 * qJD(2) + t95, t84, 0, -t94, 0;];
Cq = t8;
