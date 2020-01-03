% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:29
% EndTime: 2019-12-31 17:01:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (335->57), mult. (787->87), div. (0->0), fcn. (562->6), ass. (0->43)
t38 = cos(qJ(4));
t63 = t38 / 0.2e1;
t36 = sin(qJ(4));
t62 = -Ifges(5,2) * t38 / 0.2e1 - Ifges(5,4) * t36 + Ifges(5,1) * t63;
t35 = cos(pkin(7));
t39 = cos(qJ(2));
t34 = sin(pkin(7));
t37 = sin(qJ(2));
t55 = t34 * t37;
t13 = (t35 * t39 - t55) * pkin(1);
t49 = t36 ^ 2 + t38 ^ 2;
t61 = (t49 * mrSges(5,3) - mrSges(4,2)) * t13 + (-t37 * mrSges(3,1) - t39 * mrSges(3,2)) * pkin(1);
t60 = -t13 / 0.2e1;
t51 = t38 * mrSges(5,2);
t53 = t36 * mrSges(5,1);
t19 = -t51 - t53;
t27 = t39 * pkin(1) + pkin(2);
t41 = -pkin(1) * t55 + t35 * t27;
t9 = -pkin(3) - t41;
t58 = t9 * t19;
t26 = -t35 * pkin(2) - pkin(3);
t56 = t26 * t19;
t54 = t35 * t37;
t43 = pkin(1) * t54 + t34 * t27;
t10 = pkin(6) + t43;
t12 = (t34 * t39 + t54) * pkin(1);
t18 = -t38 * mrSges(5,1) + t36 * mrSges(5,2);
t46 = t49 * t13;
t1 = (-mrSges(4,1) + t18) * t12 + m(5) * (t10 * t46 + t9 * t12) + m(4) * (-t41 * t12 + t43 * t13) + t61;
t48 = t1 * qJD(1);
t31 = Ifges(5,4) * t38;
t21 = -Ifges(5,2) * t36 + t31;
t22 = Ifges(5,1) * t36 + t31;
t44 = (t21 + t22) * t63 + t62 * t36;
t4 = t44 - t58;
t47 = t4 * qJD(1);
t45 = Ifges(5,5) * t38 - Ifges(5,6) * t36;
t2 = (t9 / 0.2e1 + t26 / 0.2e1) * t19 + (mrSges(5,2) * t60 - t22 / 0.2e1 - t21 / 0.2e1) * t38 + (mrSges(5,1) * t60 - t62) * t36;
t5 = t44 - t56;
t42 = t2 * qJD(1) - t5 * qJD(2);
t25 = t34 * pkin(2) + pkin(6);
t3 = -t58 / 0.2e1 - t56 / 0.2e1 + (-t51 / 0.2e1 - t53 / 0.2e1) * t13 + t44;
t6 = [t1 * qJD(2) + t4 * qJD(4), t3 * qJD(4) + t48 + (-t12 * mrSges(4,1) + t12 * t18 + m(5) * (t26 * t12 + t25 * t46) + m(4) * (-t12 * t35 + t13 * t34) * pkin(2) + t61) * qJD(2), 0, t47 + t3 * qJD(2) + (t18 * t10 + t45) * qJD(4); -t2 * qJD(4) - t48, t5 * qJD(4), 0, (t18 * t25 + t45) * qJD(4) - t42; 0, 0, 0, t19 * qJD(4); t2 * qJD(2) - t47, t42, 0, 0;];
Cq = t6;
