% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:50
% EndTime: 2019-12-05 15:04:51
% DurationCPUTime: 0.33s
% Computational Cost: add. (654->51), mult. (1942->92), div. (0->0), fcn. (2075->8), ass. (0->45)
t44 = cos(qJ(5));
t42 = t44 ^ 2;
t43 = sin(qJ(5));
t60 = -t43 ^ 2 - t42;
t61 = t44 * mrSges(6,2);
t62 = t43 * mrSges(6,1);
t32 = t61 + t62;
t47 = t61 / 0.2e1 + t62 / 0.2e1;
t71 = t47 + t32 / 0.2e1;
t58 = cos(pkin(9));
t39 = -t58 * pkin(3) - pkin(4);
t53 = t44 * mrSges(6,1) - t43 * mrSges(6,2);
t68 = m(5) * pkin(3);
t70 = qJD(3) * (-m(6) * t39 + t58 * t68 + mrSges(5,1) + t53);
t56 = sin(pkin(9));
t38 = t56 * pkin(3) + pkin(6);
t69 = (-t56 * t68 + mrSges(5,2) + (m(6) * t38 + mrSges(6,3)) * t60) * qJD(3);
t66 = cos(qJ(3));
t65 = sin(qJ(3));
t59 = cos(pkin(8));
t57 = sin(pkin(8));
t30 = -t56 * t66 - t58 * t65;
t54 = (-0.1e1 - t60) * t30;
t52 = t58 * t57;
t51 = t57 * t56;
t26 = t66 * t51 + t65 * t52;
t29 = t56 * t65 - t58 * t66;
t25 = t65 * t51 - t66 * t52;
t20 = t43 * t25 - t59 * t44;
t21 = -t44 * t25 - t59 * t43;
t48 = t20 * t43 - t21 * t44 - t25;
t2 = m(6) * (t26 * t54 + t48 * t29) / 0.2e1;
t4 = m(6) * t48 * t26;
t50 = t4 * qJD(1) + t2 * qJD(2);
t7 = m(6) * t29 * t54;
t49 = t2 * qJD(1) + t7 * qJD(2);
t46 = -t32 / 0.2e1 + t47;
t13 = t39 * t32 + t42 * Ifges(6,4) + (-Ifges(6,4) * t43 + (Ifges(6,1) - Ifges(6,2)) * t44) * t43;
t15 = t46 * t29;
t5 = t46 * t26;
t45 = t5 * qJD(1) + t15 * qJD(2) - t13 * qJD(3);
t16 = t71 * t29;
t6 = t71 * t26;
t1 = t2 * qJD(3);
t3 = [t4 * qJD(3), t1, t6 * qJD(5) + t50 + (-t66 * mrSges(4,1) + t65 * mrSges(4,2)) * qJD(3) * t57 + t25 * t70 + t26 * t69, 0, t6 * qJD(3) + (-t21 * mrSges(6,1) - t20 * mrSges(6,2)) * qJD(5); t1, t7 * qJD(3), (-t65 * mrSges(4,1) - t66 * mrSges(4,2)) * qJD(3) + t16 * qJD(5) + t49 + t30 * t70 + t29 * t69, 0, t53 * qJD(5) * t30 + t16 * qJD(3); -t5 * qJD(5) - t50, -t15 * qJD(5) - t49, t13 * qJD(5), 0, (Ifges(6,5) * t44 - Ifges(6,6) * t43 - t53 * t38) * qJD(5) - t45; 0, 0, 0, 0, -t32 * qJD(5); t5 * qJD(3), t15 * qJD(3), t45, 0, 0;];
Cq = t3;
