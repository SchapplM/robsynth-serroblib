% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:28
% EndTime: 2019-12-05 14:59:29
% DurationCPUTime: 0.34s
% Computational Cost: add. (822->57), mult. (2693->113), div. (0->0), fcn. (2689->8), ass. (0->50)
t92 = m(6) * pkin(6);
t57 = sin(qJ(4));
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t76 = t56 ^ 2 + t58 ^ 2;
t69 = -0.1e1 + t76;
t91 = t57 * t69;
t45 = -t58 * mrSges(6,1) + t56 * mrSges(6,2);
t90 = -mrSges(5,1) + t45;
t89 = t56 * mrSges(6,1) + t58 * mrSges(6,2);
t71 = cos(pkin(9));
t54 = sin(pkin(9));
t59 = cos(qJ(4));
t80 = t54 * t59;
t40 = -t56 * t80 - t58 * t71;
t41 = -t56 * t71 + t58 * t80;
t65 = t40 * t56 - t41 * t58;
t78 = t57 * t54;
t20 = (t65 + t80) * t78;
t55 = sin(pkin(8));
t67 = t55 * t71;
t72 = cos(pkin(8));
t38 = t57 * t67 + t72 * t59;
t39 = -t72 * t57 + t59 * t67;
t81 = t54 * t55;
t29 = -t39 * t56 + t58 * t81;
t30 = t39 * t58 + t56 * t81;
t64 = t29 * t56 - t30 * t58 + t39;
t84 = m(6) / 0.2e1;
t3 = (t65 * t38 + (t38 * t59 + t64 * t57) * t54) * t84;
t15 = -t65 * t59 + (-t69 * t57 ^ 2 - t59 ^ 2) * t54;
t83 = t15 / 0.2e1;
t88 = -(t20 * qJD(2) + qJD(3) * t83) * m(6) - t3 * qJD(1);
t32 = t59 * t91;
t7 = (-t38 * t91 - t64 * t59) * t84;
t87 = -t7 * qJD(1) + (-t15 * qJD(2) / 0.2e1 - t32 * qJD(3)) * m(6);
t86 = t76 * mrSges(6,3) - mrSges(5,2);
t85 = -m(6) * pkin(4) + t90;
t75 = m(6) * qJD(4);
t19 = (-pkin(4) * mrSges(6,2) + Ifges(6,4) * t58) * t58 + (-pkin(4) * mrSges(6,1) - Ifges(6,4) * t56 + (Ifges(6,1) - Ifges(6,2)) * t58) * t56;
t70 = t19 * qJD(4);
t6 = m(6) * t64 * t38;
t60 = t6 * qJD(1) + t3 * qJD(2) + t7 * qJD(3);
t28 = t89 * t59;
t13 = t89 * t78;
t11 = t75 * t83;
t5 = t89 * t38;
t2 = t7 * qJD(4);
t1 = t3 * qJD(4);
t4 = [t6 * qJD(4), t1, t2, (t85 * t39 + (mrSges(5,2) - (mrSges(6,3) + t92) * t76) * t38) * qJD(4) + t5 * qJD(5) + t60, t5 * qJD(4) + (-t30 * mrSges(6,1) - t29 * mrSges(6,2)) * qJD(5); t1, t20 * t75, t11, t13 * qJD(5) + (t85 * t59 + (-t76 * t92 - t86) * t57) * qJD(4) * t54 - t88, t13 * qJD(4) + (-t41 * mrSges(6,1) - t40 * mrSges(6,2)) * qJD(5); t2, t11, t32 * t75, -t28 * qJD(5) + (t90 * t57 + t86 * t59 + (t76 * t59 * pkin(6) - pkin(4) * t57) * m(6)) * qJD(4) - t87, t45 * qJD(5) * t57 - t28 * qJD(4); -t60, t88, t87, t19 * qJD(5), t70 + (Ifges(6,5) * t58 - Ifges(6,6) * t56 + t45 * pkin(6)) * qJD(5); 0, 0, 0, -t70, 0;];
Cq = t4;
