% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:47
% EndTime: 2019-12-31 16:44:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (656->59), mult. (1386->74), div. (0->0), fcn. (1336->4), ass. (0->35)
t59 = -Ifges(5,5) + Ifges(4,4);
t36 = sin(pkin(6));
t37 = cos(pkin(6));
t38 = sin(qJ(3));
t54 = cos(qJ(3));
t28 = t38 * t36 - t37 * t54;
t30 = t36 * t54 + t38 * t37;
t43 = -t37 * pkin(2) - pkin(1);
t9 = t28 * pkin(3) - t30 * qJ(4) + t43;
t58 = m(5) * t9 + t28 * mrSges(5,1) - t30 * mrSges(5,3);
t33 = m(5) * qJ(4) + mrSges(5,3);
t56 = t28 ^ 2;
t55 = pkin(3) * t30;
t52 = pkin(5) + qJ(2);
t48 = t28 * qJ(4);
t15 = t48 + t55;
t23 = t28 * mrSges(5,3);
t39 = t30 * mrSges(4,1) - t28 * mrSges(4,2);
t1 = t43 * t39 + t9 * t23 + (t9 * mrSges(5,1) + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,3)) * t28 - t59 * t30) * t30 + t59 * t56 + t58 * t15;
t50 = t1 * qJD(1);
t31 = t52 * t37;
t40 = t52 * t36;
t17 = t38 * t31 + t40 * t54;
t18 = t31 * t54 - t38 * t40;
t2 = (t30 ^ 2 + t56) * (mrSges(4,3) + mrSges(5,2)) + (m(3) * qJ(2) + mrSges(3,3)) * (t36 ^ 2 + t37 ^ 2) + (m(4) + m(5)) * (t17 * t30 - t18 * t28);
t49 = t2 * qJD(1);
t3 = t30 * mrSges(5,1) + t23 + 0.2e1 * (t15 / 0.4e1 + t55 / 0.4e1 + t48 / 0.4e1) * m(5) + t39;
t47 = t3 * qJD(1);
t4 = t58 * t30;
t46 = t4 * qJD(1);
t10 = m(5) * t30;
t45 = t10 * qJD(1);
t44 = t33 * qJD(3);
t6 = m(5) * t18 - t28 * mrSges(5,2);
t5 = [t2 * qJD(2) + t1 * qJD(3) - t4 * qJD(4), t49, t6 * qJD(4) + t50 + ((-m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1)) * t18 + (mrSges(4,2) - t33) * t17 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t30 + (pkin(3) * mrSges(5,2) - Ifges(5,4) - Ifges(4,5)) * t28) * qJD(3), t6 * qJD(3) - t46; t3 * qJD(3) - t10 * qJD(4) - t49, 0, t47, -t45; -t3 * qJD(2) - t50, -t47, t33 * qJD(4), t44; t10 * qJD(2) + t46, t45, -t44, 0;];
Cq = t5;
