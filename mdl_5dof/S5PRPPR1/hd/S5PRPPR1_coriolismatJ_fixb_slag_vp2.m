% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:50
% EndTime: 2019-12-05 15:21:52
% DurationCPUTime: 0.54s
% Computational Cost: add. (1755->105), mult. (4199->172), div. (0->0), fcn. (4169->6), ass. (0->68)
t66 = sin(pkin(9));
t68 = cos(pkin(9));
t70 = sin(qJ(5));
t71 = cos(qJ(5));
t53 = -t70 * t66 + t71 * t68;
t105 = t66 ^ 2;
t104 = t68 ^ 2;
t103 = m(6) / 0.2e1;
t67 = sin(pkin(8));
t101 = t67 / 0.2e1;
t69 = cos(pkin(8));
t100 = -t69 / 0.2e1;
t99 = t69 / 0.2e1;
t54 = t71 * t66 + t70 * t68;
t46 = t54 * t67;
t35 = t69 * mrSges(6,2) - t46 * mrSges(6,3);
t97 = t46 * t35;
t48 = t53 * t67;
t36 = -t69 * mrSges(6,1) - t48 * mrSges(6,3);
t96 = t48 * t36;
t95 = t53 * t46;
t94 = t54 * t48;
t93 = t66 * t67;
t92 = t67 * t68;
t88 = -Ifges(6,5) * t46 - Ifges(6,6) * t48;
t58 = -t69 * pkin(3) - t67 * qJ(4) - pkin(2);
t87 = qJ(3) * t69;
t38 = t66 * t58 + t68 * t87;
t52 = t68 * t58;
t30 = -pkin(6) * t92 + t52 + (-qJ(3) * t66 - pkin(4)) * t69;
t34 = -pkin(6) * t93 + t38;
t17 = t71 * t30 - t70 * t34;
t18 = t70 * t30 + t71 * t34;
t37 = -t66 * t87 + t52;
t55 = t69 * mrSges(5,2) - mrSges(5,3) * t93;
t56 = -t69 * mrSges(5,1) - mrSges(5,3) * t92;
t7 = -t97 - t96 + m(6) * (-t17 * t48 - t18 * t46) + (-t66 * t55 - t68 * t56 + m(5) * (-t37 * t68 - t38 * t66)) * t67;
t86 = t7 * qJD(2);
t47 = t54 * t69;
t49 = t53 * t69;
t79 = -t104 - t105;
t12 = (t47 * t46 + t49 * t48) * t103 + 0.2e1 * (m(5) * (-t79 - 0.1e1) / 0.4e1 - m(6) / 0.4e1) * t69 * t67;
t85 = t12 * qJD(2);
t78 = (-t54 * t46 - t53 * t48) * t103;
t14 = t78 + (-m(6) / 0.2e1 + (-t105 / 0.2e1 - t104 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t67;
t84 = t14 * qJD(2);
t24 = t48 * mrSges(6,1) - t46 * mrSges(6,2);
t83 = t24 * qJD(2);
t82 = t24 * qJD(5);
t81 = m(5) * t101;
t57 = (pkin(4) * t66 + qJ(3)) * t67;
t1 = t88 * t100 + t17 * t35 - t18 * t36 + t57 * t24 - (-t17 * mrSges(6,3) - Ifges(6,4) * t46 + Ifges(6,5) * t100) * t46 + (Ifges(6,6) * t99 - t18 * mrSges(6,3) - Ifges(6,4) * t48 - (Ifges(6,1) - Ifges(6,2)) * t46) * t48;
t3 = t24 * t99 + t97 / 0.2e1 + t96 / 0.2e1 + (t48 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * mrSges(6,3);
t77 = -t3 * qJD(1) + t1 * qJD(2);
t64 = t67 ^ 2;
t63 = t64 * qJ(3);
t65 = t69 ^ 2;
t6 = t49 * t35 - t47 * t36 + (t68 * t55 - t66 * t56) * t69 + (t65 + t64) * mrSges(4,3) + (t46 * mrSges(6,1) + t48 * mrSges(6,2) + (mrSges(5,1) * t66 + mrSges(5,2) * t68) * t67) * t67 + m(6) * (-t17 * t47 + t18 * t49 + t57 * t67) + m(5) * (t63 + (-t37 * t66 + t38 * t68) * t69) + m(4) * (t65 * qJ(3) + t63);
t76 = t12 * qJD(1) + t6 * qJD(2);
t75 = -t47 * mrSges(6,1) / 0.2e1 - t49 * mrSges(6,2) / 0.2e1;
t74 = -t53 * t35 / 0.2e1 + t54 * t36 / 0.2e1;
t73 = t3 * qJD(2);
t4 = (t94 / 0.2e1 - t95 / 0.2e1) * mrSges(6,3) + t74 + t75;
t72 = t4 * qJD(2);
t13 = m(6) * t101 + t79 * t81 + t78 + t81;
t5 = -t74 + t75 - (t94 - t95) * mrSges(6,3) / 0.2e1;
t2 = t12 * qJD(3) - t3 * qJD(5);
t8 = [0, t2, t85, 0, -t73 - t82; t2, t6 * qJD(3) + t7 * qJD(4) + t1 * qJD(5), m(6) * (-t53 * t47 + t54 * t49) * qJD(3) + t13 * qJD(4) + t5 * qJD(5) + t76, t13 * qJD(3) + t86, t5 * qJD(3) + (-t18 * mrSges(6,1) - t17 * mrSges(6,2) + t88) * qJD(5) + t77; -t85, t14 * qJD(4) - t4 * qJD(5) - t76, 0, t84, (-t54 * mrSges(6,1) - t53 * mrSges(6,2)) * qJD(5) - t72; 0, -t14 * qJD(3) + t82 - t86, -t84, 0, t83; t73, t4 * qJD(3) - t24 * qJD(4) - t77, t72, -t83, 0;];
Cq = t8;
