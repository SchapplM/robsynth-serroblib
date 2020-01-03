% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:21
% EndTime: 2019-12-31 16:24:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (565->64), mult. (1468->102), div. (0->0), fcn. (1408->6), ass. (0->44)
t40 = sin(qJ(2));
t41 = cos(qJ(2));
t61 = t41 * t40;
t38 = sin(pkin(7));
t39 = cos(pkin(7));
t56 = sin(qJ(4));
t57 = cos(qJ(4));
t26 = -t56 * t38 + t57 * t39;
t27 = -t57 * t38 - t56 * t39;
t60 = t27 ^ 2;
t59 = m(5) / 0.2e1;
t58 = t40 / 0.2e1;
t55 = pkin(5) + qJ(3);
t36 = t38 ^ 2;
t37 = t39 ^ 2;
t54 = t36 + t37;
t19 = t27 * t40;
t20 = t27 * t41;
t21 = t26 * t40;
t22 = t26 * t41;
t6 = m(4) * (-0.1e1 + t54) * t61 + m(5) * (t19 * t20 + t21 * t22 - t61);
t53 = t6 * qJD(1);
t14 = -t27 * mrSges(5,1) + t26 * mrSges(5,2);
t52 = t14 * qJD(2);
t51 = m(4) * t58;
t50 = -t41 * t14 / 0.2e1;
t47 = t54 * mrSges(4,3);
t46 = t54 * qJ(3);
t33 = -t39 * pkin(3) - pkin(2);
t1 = -Ifges(5,4) * t60 + t33 * t14 + (Ifges(5,4) * t26 + (-Ifges(5,1) + Ifges(5,2)) * t27) * t26;
t43 = t20 * mrSges(5,1) / 0.2e1 - t22 * mrSges(5,2) / 0.2e1;
t2 = t50 - t43;
t45 = t2 * qJD(1) + t1 * qJD(2);
t28 = t55 * t38;
t29 = t55 * t39;
t15 = -t57 * t28 - t56 * t29;
t16 = -t56 * t28 + t57 * t29;
t5 = (t26 ^ 2 + t60) * mrSges(5,3) + t47 + m(5) * (t15 * t27 + t16 * t26) + m(4) * t46;
t42 = (t27 * t19 + t26 * t21) * t59;
t8 = t42 + (-m(5) / 0.2e1 + (t36 / 0.2e1 + t37 / 0.2e1 - 0.1e1 / 0.2e1) * m(4)) * t40;
t44 = t8 * qJD(1) + t5 * qJD(2);
t7 = m(5) * t58 + t54 * t51 + t42 + t51;
t3 = t50 + t43;
t4 = [t6 * qJD(2), t7 * qJD(3) + t3 * qJD(4) + t53 + ((t20 * t27 + t22 * t26) * mrSges(5,3) + (-mrSges(3,2) + t47) * t41 + (-t39 * mrSges(4,1) - t26 * mrSges(5,1) + t38 * mrSges(4,2) - t27 * mrSges(5,2) - mrSges(3,1)) * t40 + m(4) * (-t40 * pkin(2) + t41 * t46) + 0.2e1 * (t15 * t20 + t16 * t22 + t33 * t40) * t59) * qJD(2), t7 * qJD(2), t3 * qJD(2) + (-t21 * mrSges(5,1) - t19 * mrSges(5,2)) * qJD(4); t8 * qJD(3) + t2 * qJD(4) - t53, t5 * qJD(3) + t1 * qJD(4), t44, (-t16 * mrSges(5,1) - t15 * mrSges(5,2) + Ifges(5,5) * t26 + Ifges(5,6) * t27) * qJD(4) + t45; -t8 * qJD(2), t14 * qJD(4) - t44, 0, t52; -t2 * qJD(2), -t14 * qJD(3) - t45, -t52, 0;];
Cq = t4;
