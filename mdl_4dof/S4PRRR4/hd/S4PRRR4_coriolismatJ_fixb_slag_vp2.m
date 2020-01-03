% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR4
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:31
% EndTime: 2019-12-31 16:32:32
% DurationCPUTime: 0.34s
% Computational Cost: add. (903->49), mult. (1971->73), div. (0->0), fcn. (1801->4), ass. (0->29)
t53 = sin(qJ(4));
t54 = sin(qJ(3));
t55 = cos(qJ(4));
t56 = cos(qJ(3));
t46 = -t53 * t56 - t55 * t54;
t80 = -pkin(6) - pkin(5);
t49 = t80 * t54;
t50 = t80 * t56;
t33 = t53 * t49 - t55 * t50;
t73 = t33 * mrSges(5,3);
t88 = Ifges(5,4) * t46 + t73;
t45 = -t53 * t54 + t55 * t56;
t61 = t55 * t49 + t53 * t50;
t4 = -t33 * mrSges(5,1) - t61 * mrSges(5,2) + Ifges(5,5) * t45 + Ifges(5,6) * t46;
t62 = -t46 * mrSges(5,1) + t45 * mrSges(5,2);
t85 = t62 * qJD(4);
t79 = t54 * pkin(3);
t70 = Ifges(5,1) - Ifges(5,2);
t64 = -t56 * pkin(3) - pkin(2);
t69 = -t46 * t73 - t64 * t62;
t60 = t54 * mrSges(4,1) + t56 * mrSges(4,2);
t1 = pkin(2) * t60 - m(5) * t64 * t79 + (t54 ^ 2 - t56 ^ 2) * Ifges(4,4) + t69 + (-Ifges(4,1) + Ifges(4,2)) * t56 * t54 + (mrSges(5,2) * t79 + t88) * t46 + (mrSges(5,1) * t79 - Ifges(5,4) * t45 + t70 * t46) * t45;
t59 = t1 * qJD(2);
t2 = -Ifges(5,4) * t45 ^ 2 + (t70 * t45 + t88) * t46 + t69;
t58 = t2 * qJD(2);
t48 = (mrSges(5,1) * t53 + mrSges(5,2) * t55) * pkin(3);
t57 = t48 * qJD(3);
t44 = t48 * qJD(4);
t3 = [0, 0, -t85 + (-t60 - t62 + (t45 * t53 + t46 * t55) * pkin(3) * m(5)) * qJD(3), -t62 * qJD(3) - t85; 0, -t1 * qJD(3) - t2 * qJD(4), t4 * qJD(4) - t59 + (Ifges(4,5) * t56 - Ifges(4,6) * t54 + (m(5) * (-t33 * t55 + t53 * t61) + (-t55 * t45 + t53 * t46) * mrSges(5,3)) * pkin(3) + (-t56 * mrSges(4,1) + t54 * mrSges(4,2)) * pkin(5) + t4) * qJD(3), -t58 + (qJD(4) + qJD(3)) * t4; 0, t59, -t44, -t44 - t57; 0, t58, t57, 0;];
Cq = t3;
