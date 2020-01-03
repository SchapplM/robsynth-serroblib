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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:25:11
% EndTime: 2020-01-03 11:25:14
% DurationCPUTime: 0.55s
% Computational Cost: add. (1400->115), mult. (3002->170), div. (0->0), fcn. (2403->6), ass. (0->68)
t59 = cos(qJ(4));
t97 = t59 ^ 2;
t58 = sin(qJ(4));
t98 = t58 ^ 2;
t100 = -t97 - t98 + 0.1e1;
t96 = m(6) * pkin(4);
t99 = -mrSges(6,1) - t96;
t95 = -mrSges(6,2) / 0.2e1;
t94 = -t58 / 0.2e1;
t56 = sin(pkin(8));
t93 = m(6) * t56;
t53 = sin(pkin(7)) * pkin(1) + qJ(3);
t92 = t53 * t58;
t91 = t56 * t58;
t90 = t56 * t59;
t57 = cos(pkin(8));
t39 = -t57 * mrSges(6,1) - mrSges(6,3) * t90;
t89 = t58 * t39;
t38 = t57 * mrSges(6,2) - mrSges(6,3) * t91;
t88 = t59 * t38;
t87 = mrSges(5,2) + mrSges(6,2);
t86 = Ifges(6,4) + Ifges(5,4);
t85 = Ifges(5,5) + Ifges(6,5);
t84 = Ifges(6,6) + Ifges(5,6);
t52 = mrSges(6,1) * t90;
t74 = t59 * t96;
t83 = -t56 * t74 - t52;
t82 = qJ(5) * t56;
t37 = -cos(pkin(7)) * pkin(1) - pkin(2) - t56 * pkin(6) - t57 * pkin(3);
t29 = t59 * t37;
t65 = -t59 * t82 + t29;
t13 = (-pkin(4) - t92) * t57 + t65;
t73 = t57 * t92;
t14 = t65 - t73;
t69 = m(6) * (-t13 + t14);
t60 = t39 / 0.2e1 - t69 / 0.2e1;
t71 = mrSges(6,3) * t56 / 0.2e1;
t3 = ((t38 / 0.2e1 + t58 * t71) * t58 + (t59 * t71 + t60) * t59) * t56 + (t52 / 0.2e1 + (mrSges(6,2) * t94 + t74 / 0.2e1) * t56) * t57;
t81 = t3 * qJD(1);
t66 = -mrSges(6,1) / 0.2e1 - t96 / 0.2e1;
t4 = (-t38 / 0.2e1 + (t95 - mrSges(5,2)) * t57) * t59 + ((-mrSges(5,1) + t66) * t57 + t60) * t58;
t80 = t4 * qJD(1);
t17 = t59 * t57 * t53 + t58 * t37;
t15 = -t58 * t82 + t17;
t9 = (m(6) * (-t13 * t59 - t15 * t58) - t58 * t38 - t59 * t39) * t56;
t79 = t9 * qJD(1);
t78 = qJD(4) * t56;
t12 = -0.2e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t100 * t57 * t56;
t77 = t12 * qJD(1);
t20 = mrSges(6,2) * t91 + t83;
t76 = t20 * qJD(1);
t22 = (-t98 / 0.2e1 - t97 / 0.2e1 - 0.1e1 / 0.2e1) * t93;
t75 = t22 * qJD(1);
t31 = (pkin(4) * t58 + t53) * t56;
t68 = m(6) * t31 + (t58 * mrSges(6,1) + t59 * mrSges(6,2)) * t56;
t16 = t29 - t73;
t64 = t17 * mrSges(5,1) + t16 * mrSges(5,2);
t63 = t58 * mrSges(5,1) + t59 * mrSges(5,2);
t1 = t14 * t38 + t31 * t52 + t64 * t57 + (t69 - t39) * t15 + ((-t31 * mrSges(6,2) + t13 * mrSges(6,3) + t85 * t57 + (-t53 * mrSges(5,2) + t86 * t58) * t56) * t58 + (-t15 * mrSges(6,3) + t84 * t57 + (t53 * mrSges(5,1) - t86 * t59) * t56 + t68 * pkin(4) + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t91) * t59) * t56;
t62 = t1 * qJD(1) - t3 * qJD(2);
t40 = t56 ^ 2 * t53;
t55 = t57 ^ 2;
t6 = t55 * mrSges(4,3) + m(5) * t40 + m(4) * (t55 * t53 + t40) + (t88 - t89 + t63 * t57 + m(6) * (-t58 * t13 + t15 * t59) + m(5) * (-t16 * t58 + t17 * t59)) * t57 + ((mrSges(4,3) + t63) * t56 + t68) * t56;
t61 = t6 * qJD(1) + t12 * qJD(2);
t21 = t100 * t93 / 0.2e1;
t5 = t58 * t69 / 0.2e1 + t88 / 0.2e1 - t89 / 0.2e1 - t59 * mrSges(5,3) * t91 / 0.2e1 + (-t57 * mrSges(5,1) - mrSges(5,3) * t90) * t94 + (t95 * t59 + (-mrSges(5,1) / 0.2e1 + t66) * t58) * t57;
t2 = t12 * qJD(3) - t3 * qJD(4);
t7 = [t6 * qJD(3) + t1 * qJD(4) + t9 * qJD(5), t2, t5 * qJD(4) + t21 * qJD(5) + t61, t5 * qJD(3) + (-t14 * mrSges(6,2) + t99 * t15 - t64) * qJD(4) + (-t84 * t59 + (mrSges(6,3) * pkin(4) - t85) * t58) * t78 + t62, t21 * qJD(3) + t79; t2, 0, t77, -t81 + t83 * qJD(4) + (-mrSges(5,1) * t59 + t87 * t58) * t78, 0; -t4 * qJD(4) + t22 * qJD(5) - t61, -t77, 0, -t80 + (-t87 * t59 + (-mrSges(5,1) + t99) * t58) * qJD(4), t75; t4 * qJD(3) + t20 * qJD(5) - t62, t81, t80, 0, t76; -t22 * qJD(3) - t20 * qJD(4) - t79, 0, -t75, -t76, 0;];
Cq = t7;
