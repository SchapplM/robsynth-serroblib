% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->40), mult. (298->57), div. (0->0), fcn. (176->4), ass. (0->28)
t23 = sin(qJ(3));
t39 = m(6) * t23;
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t15 = t22 * mrSges(6,1) + t24 * mrSges(6,2);
t38 = mrSges(5,3) + t15;
t21 = t24 ^ 2;
t37 = m(6) / 0.2e1;
t25 = cos(qJ(3));
t20 = t22 ^ 2;
t33 = t20 + t21;
t36 = (0.1e1 - t33) * t25 * t39;
t35 = t22 * mrSges(6,2);
t34 = t24 * mrSges(6,1);
t28 = t34 - t35;
t1 = t21 * Ifges(6,4) - qJ(4) * t28 + (-Ifges(6,4) * t22 + (Ifges(6,1) - Ifges(6,2)) * t24) * t22;
t32 = t1 * qJD(3);
t31 = t15 * qJD(5);
t30 = qJD(2) * t36;
t29 = t33 * t23;
t6 = (0.1e1 / 0.2e1 - t20 / 0.2e1 - t21 / 0.2e1) * t39;
t8 = (m(5) + m(6)) * qJ(4) + t38;
t27 = t6 * qJD(2) + t8 * qJD(3);
t26 = -pkin(3) - pkin(6);
t19 = qJ(4) * t25;
t4 = t29 * t37 + (m(5) + t37) * t23;
t3 = (t28 / 0.2e1 - t35 / 0.2e1 + t34 / 0.2e1) * t23;
t2 = [0, 0, 0, 0, t28 * qJD(5); 0, qJD(3) * t36, t4 * qJD(4) + t3 * qJD(5) + t30 + ((-t33 * mrSges(6,3) - mrSges(4,1) + mrSges(5,2)) * t23 + 0.2e1 * (t26 * t29 + t19) * t37 + m(5) * (-pkin(3) * t23 + t19) + (-mrSges(4,2) + t38) * t25) * qJD(3), t4 * qJD(3), t3 * qJD(3) + t25 * t31; 0, t6 * qJD(4) - t30, t8 * qJD(4) - t1 * qJD(5), t27, -t32 + ((-mrSges(6,2) * t26 - Ifges(6,6)) * t24 + (-mrSges(6,1) * t26 - Ifges(6,5)) * t22) * qJD(5); 0, -t6 * qJD(3), -t27, 0, -t31; 0, 0, t32, 0, 0;];
Cq = t2;
