% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:02
% EndTime: 2019-12-05 15:03:02
% DurationCPUTime: 0.18s
% Computational Cost: add. (217->44), mult. (540->70), div. (0->0), fcn. (482->6), ass. (0->34)
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t17 = t22 * mrSges(6,1) + t24 * mrSges(6,2);
t40 = -mrSges(5,3) - t17;
t20 = t24 ^ 2;
t39 = m(5) / 0.2e1;
t21 = sin(pkin(8));
t23 = sin(qJ(3));
t31 = cos(pkin(8));
t37 = cos(qJ(3));
t15 = t37 * t21 + t23 * t31;
t38 = m(6) * t15;
t36 = t22 * mrSges(6,2);
t35 = t24 * mrSges(6,1);
t19 = t22 ^ 2;
t34 = t19 + t20;
t14 = t23 * t21 - t37 * t31;
t33 = qJ(4) * t14;
t29 = t34 * t15;
t1 = m(6) * (-t15 + t29) * t14;
t32 = t1 * qJD(1);
t30 = t17 * qJD(5);
t16 = -t35 + t36;
t3 = -t20 * Ifges(6,4) - qJ(4) * t16 + (Ifges(6,4) * t22 + (-Ifges(6,1) + Ifges(6,2)) * t24) * t22;
t26 = -t36 / 0.2e1 + t35 / 0.2e1;
t4 = (t16 / 0.2e1 + t26) * t15;
t28 = t4 * qJD(1) - t3 * qJD(3);
t13 = (m(5) + m(6)) * qJ(4) - t40;
t6 = 0.2e1 * (t19 / 0.4e1 + t20 / 0.4e1 - 0.1e1 / 0.4e1) * t38;
t27 = t6 * qJD(1) - t13 * qJD(3);
t25 = -pkin(3) - pkin(6);
t5 = (-t16 / 0.2e1 + t26) * t15;
t2 = t38 / 0.2e1 + 0.2e1 * (t39 + m(6) * t34 / 0.4e1) * t15;
t7 = [t1 * qJD(3), 0, t2 * qJD(4) + t5 * qJD(5) + t32 + ((mrSges(4,2) + t40) * t14 + (-t34 * mrSges(6,3) - mrSges(4,1) + mrSges(5,2)) * t15 + 0.2e1 * (-pkin(3) * t15 - t33) * t39 + m(6) * (t25 * t29 - t33)) * qJD(3), t2 * qJD(3), t5 * qJD(3) - t14 * t30; 0, 0, 0, 0, t16 * qJD(5); -t6 * qJD(4) - t4 * qJD(5) - t32, 0, t13 * qJD(4) + t3 * qJD(5), -t27, ((-mrSges(6,2) * t25 - Ifges(6,6)) * t24 + (-mrSges(6,1) * t25 - Ifges(6,5)) * t22) * qJD(5) - t28; t6 * qJD(3), 0, t27, 0, -t30; t4 * qJD(3), 0, t28, 0, 0;];
Cq = t7;
