% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:04
% EndTime: 2019-12-31 17:38:04
% DurationCPUTime: 0.26s
% Computational Cost: add. (396->64), mult. (929->97), div. (0->0), fcn. (810->6), ass. (0->46)
t40 = cos(qJ(5));
t58 = t40 * mrSges(6,2);
t38 = sin(qJ(5));
t59 = t38 * mrSges(6,1);
t25 = -t58 - t59;
t46 = -t59 / 0.2e1 - t58 / 0.2e1;
t67 = t25 / 0.2e1 + t46;
t66 = -t25 / 0.2e1 + t46;
t65 = -m(4) * qJ(3) - mrSges(4,3);
t34 = t40 ^ 2;
t64 = m(6) / 0.2e1;
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t20 = t39 * t35 + t41 * t36;
t61 = t35 * t20;
t21 = -t41 * t35 + t39 * t36;
t60 = t36 * t21;
t42 = -pkin(2) - pkin(3);
t57 = t36 * qJ(3) + t35 * t42;
t56 = t38 ^ 2 + t34;
t51 = m(6) * (-0.1e1 + t56);
t3 = t21 * t20 * t51;
t55 = t3 * qJD(1);
t14 = t67 * t36;
t54 = t14 * qJD(2);
t53 = t20 * t56;
t52 = t56 * t36;
t50 = t40 * mrSges(6,1) - t38 * mrSges(6,2);
t49 = -t35 * qJ(3) + t36 * t42;
t43 = (t35 * t53 + t60) * t64;
t44 = m(6) * (t21 * t52 + t61);
t2 = t43 - t44 / 0.2e1;
t23 = pkin(4) - t49;
t24 = -pkin(6) + t57;
t4 = (-m(6) * t24 + mrSges(6,3)) * t52 + (-m(5) * t57 - mrSges(5,2)) * t36 + (m(5) * t49 - m(6) * t23 - mrSges(5,1) - t50) * t35 + t65;
t48 = -t2 * qJD(1) - t4 * qJD(2);
t5 = t34 * Ifges(6,4) + t23 * t25 + (-Ifges(6,4) * t38 + (Ifges(6,1) - Ifges(6,2)) * t40) * t38;
t8 = t66 * t20;
t47 = -t8 * qJD(1) + t5 * qJD(2);
t45 = t50 * qJD(5);
t15 = t66 * t36;
t9 = t67 * t20;
t1 = m(4) * t39 + m(5) * (t60 + t61) + t44 / 0.2e1 + t43;
t6 = [t3 * qJD(2), t1 * qJD(3) + t9 * qJD(5) + t55 + (-t21 * mrSges(5,1) - t21 * t50 + (-t56 * mrSges(6,3) + mrSges(5,2)) * t20 + m(5) * (t57 * t20 + t49 * t21) + 0.2e1 * (-t23 * t21 + t24 * t53) * t64 + (-mrSges(3,2) - t65) * t41 + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * t39) * qJD(2), t1 * qJD(2), 0, t9 * qJD(2) - t21 * t45; -t2 * qJD(3) - t8 * qJD(5) - t55, -t4 * qJD(3) + t5 * qJD(5), t36 * qJD(3) * t35 * t51 + t15 * qJD(5) + t48, 0, t15 * qJD(3) + (-Ifges(6,5) * t40 + Ifges(6,6) * t38 - t50 * t24) * qJD(5) + t47; t2 * qJD(2), -t14 * qJD(5) - t48, 0, 0, -t35 * t45 - t54; 0, 0, 0, 0, t25 * qJD(5); t8 * qJD(2), t14 * qJD(3) - t47, t54, 0, 0;];
Cq = t6;
