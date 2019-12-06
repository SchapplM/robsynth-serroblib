% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:36
% EndTime: 2019-12-05 17:35:38
% DurationCPUTime: 0.50s
% Computational Cost: add. (1400->111), mult. (3002->162), div. (0->0), fcn. (2403->6), ass. (0->66)
t57 = sin(qJ(4));
t97 = t57 / 0.2e1;
t58 = cos(qJ(4));
t94 = t58 ^ 2;
t95 = t57 ^ 2;
t96 = -t94 - t95 + 0.1e1;
t93 = m(6) * pkin(4);
t70 = mrSges(6,1) + t93;
t92 = -mrSges(6,2) / 0.2e1;
t55 = sin(pkin(8));
t91 = m(6) * t55;
t52 = sin(pkin(7)) * pkin(1) + qJ(3);
t90 = t52 * t57;
t89 = t55 * t57;
t88 = t55 * t58;
t56 = cos(pkin(8));
t39 = -t56 * mrSges(6,1) - mrSges(6,3) * t88;
t87 = t57 * t39;
t38 = t56 * mrSges(6,2) - mrSges(6,3) * t89;
t86 = t58 * t38;
t85 = Ifges(6,4) + Ifges(5,4);
t84 = Ifges(5,5) + Ifges(6,5);
t83 = Ifges(6,6) + Ifges(5,6);
t82 = qJ(5) * t55;
t37 = -cos(pkin(7)) * pkin(1) - pkin(2) - t55 * pkin(6) - t56 * pkin(3);
t29 = t58 * t37;
t64 = -t58 * t82 + t29;
t13 = (-pkin(4) - t90) * t56 + t64;
t17 = t52 * t56 * t58 + t57 * t37;
t15 = -t57 * t82 + t17;
t9 = (m(6) * (-t13 * t58 - t15 * t57) - t57 * t38 - t58 * t39) * t55;
t81 = qJD(1) * t9;
t51 = mrSges(6,2) * t89;
t74 = t56 * t90;
t14 = t64 - t74;
t69 = m(6) * (-t13 + t14);
t59 = t39 / 0.2e1 - t69 / 0.2e1;
t65 = mrSges(6,1) / 0.2e1 + t93 / 0.2e1;
t72 = mrSges(6,3) * t55 / 0.2e1;
t3 = ((t38 / 0.2e1 + t57 * t72) * t57 + (t58 * t72 + t59) * t58) * t55 + (-t51 / 0.2e1 + t65 * t88) * t56;
t80 = t3 * qJD(1);
t4 = (-t38 / 0.2e1 + (-mrSges(5,2) + t92) * t56) * t58 + ((-mrSges(5,1) - t65) * t56 + t59) * t57;
t79 = t4 * qJD(1);
t78 = qJD(4) * t55;
t12 = -0.2e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t96 * t56 * t55;
t77 = t12 * qJD(1);
t20 = t70 * t88 - t51;
t76 = t20 * qJD(1);
t22 = (-t95 / 0.2e1 - t94 / 0.2e1 - 0.1e1 / 0.2e1) * t91;
t75 = t22 * qJD(1);
t31 = (pkin(4) * t57 + t52) * t55;
t68 = m(6) * t31 + (t57 * mrSges(6,1) + t58 * mrSges(6,2)) * t55;
t66 = -mrSges(5,1) - t70;
t16 = t29 - t74;
t63 = t17 * mrSges(5,1) + t16 * mrSges(5,2);
t62 = t57 * mrSges(5,1) + t58 * mrSges(5,2);
t1 = t14 * t38 - t31 * t51 + t63 * t56 + (t69 - t39) * t15 + ((t13 * mrSges(6,3) + t84 * t56 + (-t52 * mrSges(5,2) + t85 * t57) * t55) * t57 + (t31 * mrSges(6,1) - t15 * mrSges(6,3) + t83 * t56 + (t52 * mrSges(5,1) - t85 * t58) * t55 + t68 * pkin(4) + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t89) * t58) * t55;
t61 = t1 * qJD(1) - t3 * qJD(2);
t40 = t55 ^ 2 * t52;
t54 = t56 ^ 2;
t6 = t54 * mrSges(4,3) + m(5) * t40 + m(4) * (t52 * t54 + t40) + (t86 - t87 + t62 * t56 + m(6) * (-t13 * t57 + t15 * t58) + m(5) * (-t16 * t57 + t17 * t58)) * t56 + ((mrSges(4,3) + t62) * t55 + t68) * t55;
t60 = qJD(1) * t6 + qJD(2) * t12;
t21 = t96 * t91 / 0.2e1;
t5 = t69 * t97 + t86 / 0.2e1 - t87 / 0.2e1 + (-t65 * t57 + t92 * t58) * t56 + (-t58 * t89 / 0.2e1 + t88 * t97) * mrSges(5,3);
t2 = qJD(3) * t12 - qJD(4) * t3;
t7 = [qJD(3) * t6 + qJD(4) * t1 + qJD(5) * t9, t2, qJD(4) * t5 + qJD(5) * t21 + t60, t5 * qJD(3) + (-t14 * mrSges(6,2) - t70 * t15 - t63) * qJD(4) + (-t83 * t58 + (mrSges(6,3) * pkin(4) - t84) * t57) * t78 + t61, qJD(3) * t21 + t81; t2, 0, t77, -t80 + t51 * qJD(4) + (mrSges(5,2) * t57 + t66 * t58) * t78, 0; -qJD(4) * t4 + qJD(5) * t22 - t60, -t77, 0, -t79 + ((-mrSges(5,2) - mrSges(6,2)) * t58 + t66 * t57) * qJD(4), t75; qJD(3) * t4 - qJD(5) * t20 - t61, t80, t79, 0, -t76; -qJD(3) * t22 + qJD(4) * t20 - t81, 0, -t75, t76, 0;];
Cq = t7;
