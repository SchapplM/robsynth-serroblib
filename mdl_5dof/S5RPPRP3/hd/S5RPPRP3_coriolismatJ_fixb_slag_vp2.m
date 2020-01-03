% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:45
% EndTime: 2019-12-31 17:50:45
% DurationCPUTime: 0.18s
% Computational Cost: add. (358->57), mult. (618->68), div. (0->0), fcn. (373->4), ass. (0->35)
t26 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t54 = -qJ(5) + t26;
t34 = sin(qJ(4));
t32 = t34 ^ 2;
t53 = m(6) * pkin(4);
t35 = cos(qJ(4));
t52 = pkin(4) * t35;
t51 = t35 * mrSges(5,1);
t50 = mrSges(5,1) + mrSges(6,1);
t49 = mrSges(5,2) + mrSges(6,2);
t48 = Ifges(6,4) + Ifges(5,4);
t33 = t35 ^ 2;
t47 = t32 + t33;
t27 = sin(pkin(7)) * pkin(1) + qJ(3);
t38 = t34 * pkin(4) + t27;
t37 = m(6) * t38;
t31 = t35 * mrSges(6,1);
t39 = -t34 * mrSges(6,2) + t31;
t1 = -t27 * t51 - t48 * t32 - t38 * t39 + (-mrSges(6,2) * t52 - pkin(4) * t37 + t48 * t35) * t35 + (t27 * mrSges(5,2) - mrSges(6,1) * t52 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,2)) * t35) * t34;
t46 = t1 * qJD(1);
t3 = mrSges(4,3) + t49 * t35 + t50 * t34 + t37 + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t27;
t45 = t3 * qJD(1);
t11 = t54 * t34;
t12 = t54 * t35;
t4 = m(6) * (-t34 * t11 - t35 * t12) + t47 * mrSges(6,3);
t44 = t4 * qJD(1);
t43 = qJD(4) * t34;
t42 = qJD(4) * t35;
t28 = m(6) * t52;
t13 = -t28 - t39;
t41 = t13 * qJD(1);
t23 = (-t32 / 0.2e1 - t33 / 0.2e1 - 0.1e1 / 0.2e1) * m(6);
t40 = t23 * qJD(1);
t22 = (-t47 / 0.2e1 + 0.1e1 / 0.2e1) * m(6);
t2 = [t3 * qJD(3) - t1 * qJD(4) + t4 * qJD(5), 0, t22 * qJD(5) + t45, -t46 + (-t12 * mrSges(6,2) + (-mrSges(6,1) - t53) * t11) * qJD(4) + (-mrSges(5,2) * t26 - Ifges(5,6) - Ifges(6,6)) * t42 + (-mrSges(5,1) * t26 + mrSges(6,3) * pkin(4) - Ifges(5,5) - Ifges(6,5)) * t43, t22 * qJD(3) + t44; 0, 0, 0, (t49 * t34 - t28 - t31 - t51) * qJD(4), 0; t23 * qJD(5) - t45, 0, 0, -t49 * t42 + (-t50 - t53) * t43, t40; t13 * qJD(5) + t46, 0, 0, 0, t41; -t23 * qJD(3) - t13 * qJD(4) - t44, 0, -t40, -t41, 0;];
Cq = t2;
