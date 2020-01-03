% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:02
% EndTime: 2019-12-31 17:14:03
% DurationCPUTime: 0.49s
% Computational Cost: add. (685->111), mult. (1434->143), div. (0->0), fcn. (896->4), ass. (0->70)
t55 = cos(qJ(3));
t83 = Ifges(5,5) - Ifges(4,4);
t110 = t83 * t55;
t52 = t55 ^ 2;
t53 = sin(qJ(3));
t109 = t53 ^ 2 + t52;
t108 = Ifges(4,1) + Ifges(5,1);
t94 = cos(qJ(2)) * pkin(1);
t107 = t109 * t94;
t106 = -Ifges(5,3) + t108;
t34 = -t55 * mrSges(5,1) - t53 * mrSges(5,3);
t105 = -t55 * mrSges(4,1) + t53 * mrSges(4,2) + t34;
t95 = sin(qJ(2)) * pkin(1);
t104 = (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t109) * t94 + (-mrSges(3,1) + t105) * t95;
t103 = -m(5) / 0.2e1;
t102 = m(5) / 0.2e1;
t101 = -mrSges(5,1) / 0.2e1;
t76 = t53 * qJ(4);
t67 = -t55 * pkin(3) - t76;
t32 = -pkin(2) + t67;
t20 = t32 - t94;
t100 = t20 / 0.2e1;
t99 = t53 / 0.2e1;
t97 = m(5) * t55;
t85 = t55 * mrSges(4,2);
t87 = t53 * mrSges(4,1);
t36 = t85 + t87;
t96 = pkin(2) * t36;
t93 = m(5) * qJ(4);
t92 = Ifges(4,4) * t53;
t45 = -pkin(2) - t94;
t88 = t45 * t36;
t86 = t53 * mrSges(5,1);
t47 = t55 * mrSges(5,2);
t84 = t55 * mrSges(5,3);
t82 = t20 + t32;
t44 = pkin(6) + t95;
t81 = t107 * t44;
t80 = t107 * pkin(6);
t66 = t53 * pkin(3) - t55 * qJ(4);
t17 = t66 * t34;
t61 = t17 - t53 * t92 / 0.2e1 + (0.2e1 * Ifges(5,5) * t53 - t92) * t99 + (t106 * t99 + (-Ifges(4,2) - Ifges(5,3) / 0.2e1 + t108 / 0.2e1) * t53 - t110) * t55;
t68 = -t84 + t86;
t3 = t88 + (m(5) * t66 + t68) * t20 + t61;
t78 = t3 * qJD(1);
t5 = m(4) * (t45 * t95 + t81) + m(5) * (t20 * t95 + t81) + t104;
t77 = t5 * qJD(1);
t10 = (m(5) * t20 + t34) * t53;
t75 = t10 * qJD(1);
t46 = mrSges(5,3) + t93;
t74 = t46 * qJD(3);
t71 = pkin(3) * t103;
t70 = t94 / 0.2e1;
t69 = m(5) * t82;
t16 = t32 * t68;
t60 = -t83 * t53 + (Ifges(4,2) - t106) * t55;
t1 = -t17 - t16 / 0.2e1 + (-t45 / 0.2e1 + pkin(2) / 0.2e1) * t36 + (mrSges(5,3) * t100 + t110 + qJ(4) * t69 / 0.2e1 + (-mrSges(4,2) / 0.2e1 + mrSges(5,3) / 0.2e1 + t93 / 0.2e1) * t94) * t55 + (t20 * t101 + t82 * t71 + (-mrSges(4,1) / 0.2e1 + t101 + t71) * t94 + t60) * t53;
t62 = t32 * t66;
t4 = -m(5) * t62 + t83 * t52 + t60 * t53 - t16 - t17 + t96;
t64 = t1 * qJD(1) + t4 * qJD(2);
t11 = (m(5) * t32 + t34) * t53;
t6 = (t34 + (t70 + t100 + t32 / 0.2e1) * m(5)) * t53;
t63 = t6 * qJD(1) + t11 * qJD(2);
t58 = (-mrSges(5,2) * t76 - pkin(3) * t47 + (Ifges(5,4) + Ifges(4,5)) * t55 + (-Ifges(4,6) + Ifges(5,6)) * t53) * qJD(3);
t57 = qJD(3) * (m(5) * t67 + t105);
t33 = pkin(6) * t97 + t47;
t19 = t44 * t97 + t47;
t7 = -t69 * t99 + (m(5) * t70 - t34) * t53;
t2 = t61 + (-t85 / 0.2e1 - t87 / 0.2e1 - t86 / 0.2e1 + t84 / 0.2e1 + t66 * t103) * t94 + (t20 * t66 + t62) * t102 + t68 * t100 + t88 / 0.2e1 - t96 / 0.2e1 + t16 / 0.2e1;
t8 = [t5 * qJD(2) + t3 * qJD(3) - t10 * qJD(4), t2 * qJD(3) + t7 * qJD(4) + t77 + (m(4) * (-pkin(2) * t95 + t80) + 0.2e1 * (t32 * t95 + t80) * t102 + t104) * qJD(2), t2 * qJD(2) + t19 * qJD(4) + t44 * t57 + t58 + t78, t7 * qJD(2) + t19 * qJD(3) - t75; -t1 * qJD(3) - t6 * qJD(4) - t77, -t4 * qJD(3) - t11 * qJD(4), pkin(6) * t57 + t33 * qJD(4) + t58 - t64, t33 * qJD(3) - t63; t1 * qJD(2) - t78, t64, t46 * qJD(4), t74; t6 * qJD(2) + t75, t63, -t74, 0;];
Cq = t8;
