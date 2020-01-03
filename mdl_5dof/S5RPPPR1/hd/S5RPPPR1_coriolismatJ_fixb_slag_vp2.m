% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:09
% EndTime: 2020-01-03 11:20:11
% DurationCPUTime: 0.56s
% Computational Cost: add. (2138->106), mult. (4582->170), div. (0->0), fcn. (4552->8), ass. (0->69)
t67 = sin(pkin(9));
t69 = cos(pkin(9));
t71 = sin(qJ(5));
t72 = cos(qJ(5));
t55 = -t71 * t67 + t72 * t69;
t105 = t55 / 0.2e1;
t104 = t67 ^ 2;
t103 = t69 ^ 2;
t102 = m(6) / 0.2e1;
t56 = t72 * t67 + t71 * t69;
t101 = -t56 / 0.2e1;
t68 = sin(pkin(8));
t100 = t68 / 0.2e1;
t70 = cos(pkin(8));
t99 = -t70 / 0.2e1;
t98 = t70 / 0.2e1;
t46 = t56 * t68;
t37 = t70 * mrSges(6,2) - t46 * mrSges(6,3);
t96 = t46 * t37;
t48 = t55 * t68;
t38 = -t70 * mrSges(6,1) - t48 * mrSges(6,3);
t95 = t48 * t38;
t62 = sin(pkin(7)) * pkin(1) + qJ(3);
t94 = t62 * t70;
t93 = t67 * t68;
t92 = t68 * t69;
t88 = -Ifges(6,5) * t46 - Ifges(6,6) * t48;
t54 = -cos(pkin(7)) * pkin(1) - pkin(2) - t68 * qJ(4) - t70 * pkin(3);
t35 = t67 * t54 + t69 * t94;
t51 = t69 * t54;
t24 = -pkin(6) * t92 + t51 + (-t62 * t67 - pkin(4)) * t70;
t29 = -pkin(6) * t93 + t35;
t15 = t72 * t24 - t71 * t29;
t16 = t71 * t24 + t72 * t29;
t34 = -t67 * t94 + t51;
t57 = t70 * mrSges(5,2) - mrSges(5,3) * t93;
t58 = -t70 * mrSges(5,1) - mrSges(5,3) * t92;
t8 = m(6) * (-t15 * t48 - t16 * t46) - t96 - t95 + (-t67 * t57 - t69 * t58 + m(5) * (-t34 * t69 - t35 * t67)) * t68;
t87 = t8 * qJD(1);
t47 = t56 * t70;
t49 = t55 * t70;
t80 = -t103 - t104;
t12 = (t47 * t46 + t49 * t48) * t102 + 0.2e1 * (m(5) * (-t80 - 0.1e1) / 0.4e1 - m(6) / 0.4e1) * t70 * t68;
t86 = t12 * qJD(1);
t79 = (-t56 * t46 - t55 * t48) * t102;
t14 = t79 + (-m(6) / 0.2e1 + (-t104 / 0.2e1 - t103 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t68;
t85 = t14 * qJD(1);
t25 = t48 * mrSges(6,1) - t46 * mrSges(6,2);
t84 = t25 * qJD(1);
t83 = t25 * qJD(5);
t82 = m(5) * t100;
t53 = (pkin(4) * t67 + t62) * t68;
t1 = t15 * t37 - t16 * t38 + t88 * t99 + t53 * t25 - (-t15 * mrSges(6,3) - Ifges(6,4) * t46 + Ifges(6,5) * t99) * t46 + (-t16 * mrSges(6,3) + Ifges(6,6) * t98 - Ifges(6,4) * t48 - (Ifges(6,1) - Ifges(6,2)) * t46) * t48;
t3 = t25 * t98 + t96 / 0.2e1 + t95 / 0.2e1 + (t48 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * mrSges(6,3);
t78 = t1 * qJD(1) - t3 * qJD(2);
t65 = t68 ^ 2;
t60 = t65 * t62;
t66 = t70 ^ 2;
t6 = t49 * t37 - t47 * t38 + (t69 * t57 - t67 * t58) * t70 + (t66 + t65) * mrSges(4,3) + (t46 * mrSges(6,1) + t48 * mrSges(6,2) + (mrSges(5,1) * t67 + mrSges(5,2) * t69) * t68) * t68 + m(6) * (-t15 * t47 + t16 * t49 + t53 * t68) + m(5) * (t60 + (-t34 * t67 + t35 * t69) * t70) + m(4) * (t66 * t62 + t60);
t77 = t6 * qJD(1) + t12 * qJD(2);
t76 = -t47 * mrSges(6,1) / 0.2e1 - t49 * mrSges(6,2) / 0.2e1;
t75 = t3 * qJD(1);
t73 = (t48 * t101 + t46 * t105) * mrSges(6,3) + t37 * t105 + t38 * t101;
t4 = t73 - t76;
t74 = t4 * qJD(1);
t13 = m(6) * t100 + t80 * t82 + t79 + t82;
t5 = t73 + t76;
t2 = t12 * qJD(3) - t3 * qJD(5);
t7 = [t6 * qJD(3) + t8 * qJD(4) + t1 * qJD(5), t2, m(6) * (-t55 * t47 + t56 * t49) * qJD(3) + t13 * qJD(4) + t5 * qJD(5) + t77, t13 * qJD(3) + t87, t5 * qJD(3) + (-t16 * mrSges(6,1) - t15 * mrSges(6,2) + t88) * qJD(5) + t78; t2, 0, t86, 0, -t75 - t83; t14 * qJD(4) + t4 * qJD(5) - t77, -t86, 0, t85, (-t56 * mrSges(6,1) - t55 * mrSges(6,2)) * qJD(5) + t74; -t14 * qJD(3) + t83 - t87, 0, -t85, 0, t84; -t4 * qJD(3) - t25 * qJD(4) - t78, t75, -t74, -t84, 0;];
Cq = t7;
