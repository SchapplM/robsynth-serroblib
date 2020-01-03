% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:14
% EndTime: 2019-12-31 16:30:15
% DurationCPUTime: 0.22s
% Computational Cost: add. (228->55), mult. (574->78), div. (0->0), fcn. (370->4), ass. (0->31)
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t59 = t34 ^ 2 + t36 ^ 2;
t41 = t34 * pkin(3) - t36 * qJ(4);
t39 = m(5) * t41;
t37 = cos(qJ(2));
t61 = t34 * t37;
t24 = -t36 * mrSges(5,1) - t34 * mrSges(5,3);
t60 = -t36 * mrSges(4,1) + t34 * mrSges(4,2) + t24;
t52 = t36 * mrSges(5,3);
t54 = t34 * mrSges(5,1);
t58 = -t39 + t52 - t54;
t51 = Ifges(4,4) - Ifges(5,5);
t50 = t59 * pkin(5) * t37;
t49 = m(5) * qJD(4);
t35 = sin(qJ(2));
t4 = 0.4e1 * (m(4) / 0.4e1 + m(5) / 0.4e1) * (-0.1e1 + t59) * t37 * t35;
t48 = t4 * qJD(1);
t42 = -t36 * pkin(3) - t34 * qJ(4);
t22 = -pkin(2) + t42;
t5 = (m(5) * t22 + t24) * t34;
t47 = t5 * qJD(2);
t46 = qJD(3) * t36;
t31 = m(5) * qJ(4) + mrSges(5,3);
t45 = t31 * qJD(3);
t1 = -t41 * t24 + (pkin(2) * mrSges(4,2) - t51 * t36) * t36 + (pkin(2) * mrSges(4,1) + t51 * t34 + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,3)) * t36) * t34 + t58 * t22;
t40 = t1 * qJD(2);
t38 = (m(5) * t42 + t60) * qJD(3);
t23 = (m(5) * pkin(5) + mrSges(5,2)) * t36;
t3 = (-t54 / 0.2e1 + t52 / 0.2e1 - t36 * mrSges(4,2) - t34 * mrSges(4,1) - t39 / 0.2e1 + t58 / 0.2e1) * t37;
t2 = [t4 * qJD(2), t49 * t61 + t3 * qJD(3) + t48 + ((-mrSges(3,1) + t60) * t35 + m(4) * (-t35 * pkin(2) + t50) + m(5) * (t22 * t35 + t50) + (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t59) * t37) * qJD(2), t3 * qJD(2) + (t36 * t49 + t38) * t35, (qJD(2) * t61 + t35 * t46) * m(5); -t48, -t1 * qJD(3) - t5 * qJD(4), t23 * qJD(4) + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t46 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * qJD(3) * t34 + pkin(5) * t38 - t40, t23 * qJD(3) - t47; 0, t40, t31 * qJD(4), t45; 0, t47, -t45, 0;];
Cq = t2;
