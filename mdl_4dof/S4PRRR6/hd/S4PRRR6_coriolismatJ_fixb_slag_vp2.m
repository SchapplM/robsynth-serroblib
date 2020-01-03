% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:41
% EndTime: 2019-12-31 16:34:42
% DurationCPUTime: 0.56s
% Computational Cost: add. (1186->85), mult. (2978->128), div. (0->0), fcn. (2792->6), ass. (0->54)
t117 = -pkin(6) - pkin(5);
t71 = sin(qJ(3));
t61 = t117 * t71;
t74 = cos(qJ(3));
t62 = t117 * t74;
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t39 = t61 * t70 - t62 * t73;
t58 = -t70 * t74 - t73 * t71;
t87 = t70 * t71 - t73 * t74;
t90 = t73 * t61 + t62 * t70;
t7 = -t39 * mrSges(5,1) - t90 * mrSges(5,2) - Ifges(5,5) * t87 + Ifges(5,6) * t58;
t126 = t7 * qJD(4);
t72 = sin(qJ(2));
t46 = t87 * t72;
t47 = t58 * t72;
t91 = t46 * mrSges(5,1) - t47 * mrSges(5,2);
t125 = t91 * qJD(4);
t81 = -t58 * mrSges(5,1) - t87 * mrSges(5,2);
t92 = -pkin(3) * t74 - pkin(2);
t4 = (Ifges(5,4) * t58 + (-Ifges(5,1) + Ifges(5,2)) * t87) * t58 - Ifges(5,4) * t87 ^ 2 - t92 * t81;
t115 = pkin(3) * t71;
t124 = m(5) * t115;
t75 = cos(qJ(2));
t120 = t72 * t75;
t89 = -t74 * mrSges(4,1) + t71 * mrSges(4,2);
t48 = t58 * t75;
t49 = t87 * t75;
t100 = t48 * mrSges(5,1) / 0.2e1 + t49 * mrSges(5,2) / 0.2e1;
t68 = t71 ^ 2;
t69 = t74 ^ 2;
t118 = m(5) / 0.2e1;
t106 = t71 * mrSges(4,1);
t101 = t74 * mrSges(4,2);
t98 = t68 + t69;
t12 = m(4) * (-0.1e1 + t98) * t120 + m(5) * (t46 * t49 + t47 * t48 - t120);
t97 = t12 * qJD(1);
t78 = t75 * t81;
t95 = -t78 / 0.2e1;
t88 = t101 + t106;
t82 = t87 * mrSges(5,1) - t58 * mrSges(5,2);
t1 = -t82 * t115 - t92 * t124 + pkin(2) * t88 + (-t69 + t68) * Ifges(4,4) + (-Ifges(4,1) + Ifges(4,2)) * t71 * t74 + t4;
t77 = t75 * t124;
t79 = (t48 * t73 - t49 * t70) * pkin(3) * t118 + t100;
t2 = t78 / 0.2e1 + t77 / 0.2e1 + t79;
t86 = -t2 * qJD(1) - t1 * qJD(2);
t5 = t95 - t100;
t85 = t5 * qJD(1) - t4 * qJD(2);
t60 = (mrSges(5,1) * t70 + mrSges(5,2) * t73) * pkin(3);
t80 = t60 * qJD(3);
t57 = t60 * qJD(4);
t6 = t95 + t100;
t3 = -t77 / 0.2e1 + t79 + t95 + (-t88 / 0.2e1 - t101 / 0.2e1 - t106 / 0.2e1) * t75;
t8 = [t12 * qJD(2), t3 * qJD(3) + t6 * qJD(4) + t97 + ((t48 * t58 + t49 * t87) * mrSges(5,3) + (t98 * mrSges(4,3) - mrSges(3,2)) * t75 + (-mrSges(3,1) + t82 + t89) * t72 + m(4) * (t98 * t75 * pkin(5) - t72 * pkin(2)) + 0.2e1 * (-t39 * t49 + t48 * t90 + t92 * t72) * t118) * qJD(2), t3 * qJD(2) + (t72 * t89 + t91 + m(5) * (t46 * t73 + t47 * t70) * pkin(3)) * qJD(3) + t125, t6 * qJD(2) + qJD(3) * t91 + t125; -qJD(3) * t2 + qJD(4) * t5 - t97, -qJD(3) * t1 - qJD(4) * t4, t126 + t86 + (Ifges(4,5) * t74 - Ifges(4,6) * t71 + (m(5) * (-t39 * t73 + t70 * t90) + (t70 * t58 + t73 * t87) * mrSges(5,3)) * pkin(3) + t89 * pkin(5) + t7) * qJD(3), t7 * qJD(3) + t126 + t85; t2 * qJD(2), -t86, -t57, -t57 - t80; -t5 * qJD(2), -t85, t80, 0;];
Cq = t8;
