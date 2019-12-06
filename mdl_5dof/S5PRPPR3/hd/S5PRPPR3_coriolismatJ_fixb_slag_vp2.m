% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:21
% EndTime: 2019-12-05 15:26:22
% DurationCPUTime: 0.20s
% Computational Cost: add. (239->49), mult. (588->76), div. (0->0), fcn. (510->6), ass. (0->38)
t23 = sin(pkin(8));
t45 = t23 * pkin(2);
t33 = cos(pkin(8));
t44 = t33 * pkin(2);
t25 = sin(qJ(2));
t40 = cos(qJ(2));
t15 = t23 * t40 + t33 * t25;
t24 = sin(qJ(5));
t21 = t24 ^ 2;
t26 = cos(qJ(5));
t22 = t26 ^ 2;
t35 = t21 + t22;
t31 = t35 * t15;
t17 = t24 * mrSges(6,1) + t26 * mrSges(6,2);
t43 = mrSges(5,3) + t17;
t42 = m(6) / 0.4e1;
t41 = m(6) * t15;
t14 = t23 * t25 - t33 * t40;
t19 = qJ(4) + t45;
t38 = t19 * t14;
t37 = t24 * mrSges(6,2);
t36 = t26 * mrSges(6,1);
t1 = m(6) * (-t15 + t31) * t14;
t34 = t1 * qJD(1);
t32 = t17 * qJD(5);
t30 = -pkin(3) - t44;
t16 = -t36 + t37;
t3 = -t22 * Ifges(6,4) - t19 * t16 + (Ifges(6,4) * t24 + (-Ifges(6,1) + Ifges(6,2)) * t26) * t24;
t27 = -t37 / 0.2e1 + t36 / 0.2e1;
t4 = (t16 / 0.2e1 + t27) * t15;
t29 = t4 * qJD(1) - t3 * qJD(2);
t6 = 0.2e1 * (t21 / 0.4e1 + t22 / 0.4e1 - 0.1e1 / 0.4e1) * t41;
t8 = 0.4e1 * (t42 + m(5) / 0.4e1) * t19 + t43;
t28 = t6 * qJD(1) - t8 * qJD(2);
t18 = -pkin(6) + t30;
t5 = (-t16 / 0.2e1 + t27) * t15;
t2 = t41 / 0.2e1 + 0.2e1 * (m(5) / 0.2e1 + t35 * t42) * t15;
t7 = [t1 * qJD(2), t34 + (-t40 * mrSges(3,2) - t25 * mrSges(3,1) - m(5) * t38 + m(6) * (t18 * t31 - t38) - mrSges(6,3) * t31 + (-m(4) * t44 + m(5) * t30 - mrSges(4,1) + mrSges(5,2)) * t15 + (-m(4) * t45 + mrSges(4,2) - t43) * t14) * qJD(2) + t2 * qJD(4) + t5 * qJD(5), 0, t2 * qJD(2), t5 * qJD(2) - t14 * t32; -t6 * qJD(4) - t4 * qJD(5) - t34, t8 * qJD(4) + t3 * qJD(5), 0, -t28, (-Ifges(6,5) * t24 - Ifges(6,6) * t26 - t17 * t18) * qJD(5) - t29; 0, 0, 0, 0, t16 * qJD(5); t6 * qJD(2), t28, 0, 0, -t32; t4 * qJD(2), t29, 0, 0, 0;];
Cq = t7;
