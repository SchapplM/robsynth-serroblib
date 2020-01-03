% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:42
% EndTime: 2019-12-31 18:10:43
% DurationCPUTime: 0.26s
% Computational Cost: add. (425->83), mult. (769->96), div. (0->0), fcn. (500->4), ass. (0->44)
t21 = sin(pkin(7)) * pkin(1) + pkin(6);
t55 = -qJ(5) + t21;
t54 = Ifges(4,4) - Ifges(6,4) - Ifges(5,5);
t53 = m(6) / 0.2e1;
t52 = m(5) + m(6);
t51 = pkin(3) + pkin(4);
t29 = sin(qJ(3));
t50 = m(6) * t29;
t30 = cos(qJ(3));
t49 = t30 * pkin(3);
t48 = t51 * t29;
t47 = -mrSges(5,2) + mrSges(6,3);
t46 = mrSges(6,2) + mrSges(5,3);
t35 = -cos(pkin(7)) * pkin(1) - pkin(2);
t43 = t29 * qJ(4);
t32 = -t35 + t43;
t12 = -t32 - t49;
t41 = t30 * qJ(4);
t18 = -t29 * pkin(3) + t41;
t15 = -t29 * pkin(4) + t18;
t23 = t29 * mrSges(6,1);
t16 = -t30 * mrSges(5,1) - t29 * mrSges(5,3);
t34 = -m(5) * t12 - t16;
t9 = t51 * t30 + t32;
t36 = m(6) * t9 + t30 * mrSges(6,1) + t29 * mrSges(6,2);
t1 = t36 * t15 + t34 * t18 - t9 * t23 + (t35 * mrSges(4,2) + t9 * mrSges(6,2) - t12 * mrSges(5,3) + t54 * t30) * t30 + (t35 * mrSges(4,1) + t12 * mrSges(5,1) + (-Ifges(6,2) + Ifges(5,1) - Ifges(5,3) + Ifges(4,1) - Ifges(4,2) + Ifges(6,1)) * t30 - t54 * t29) * t29;
t45 = t1 * qJD(1);
t2 = (t34 + t36) * t29;
t44 = t2 * qJD(1);
t13 = t55 * t29;
t14 = t55 * t30;
t3 = m(6) * (-t13 * t29 - t14 * t30) + (t29 ^ 2 + t30 ^ 2) * mrSges(6,3);
t42 = t3 * qJD(1);
t5 = t30 * mrSges(6,2) - t23 + 0.2e1 * (t15 / 0.4e1 + t41 / 0.4e1 - t48 / 0.4e1) * m(6);
t40 = t5 * qJD(1);
t39 = qJD(3) * t30;
t20 = t52 * qJ(4) + t46;
t38 = t20 * qJD(3);
t37 = qJD(1) * t50;
t33 = t41 - t48;
t19 = t52 * t29;
t8 = -m(6) * t33 / 0.2e1 + t15 * t53;
t6 = m(6) * t14 + (m(5) * t21 - t47) * t30;
t4 = [t1 * qJD(3) + t2 * qJD(4) + t3 * qJD(5), 0, t45 + t6 * qJD(4) + t8 * qJD(5) + (-pkin(3) * mrSges(5,2) + mrSges(6,3) * t51 + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t39 + (-t14 * mrSges(6,1) - t13 * mrSges(6,2) + m(6) * (-qJ(4) * t13 - t14 * t51) + (t47 * qJ(4) - Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t29 + (t29 * mrSges(4,2) - t30 * mrSges(4,1) + m(5) * (-t43 - t49) + t16) * t21) * qJD(3), t6 * qJD(3) + t44, t8 * qJD(3) + t42; 0, 0, t19 * qJD(4) + (-mrSges(4,2) + t46) * t39 + (m(5) * t18 + 0.2e1 * t33 * t53 - t23 + (-mrSges(4,1) - mrSges(5,1)) * t29) * qJD(3), t19 * qJD(3), 0; -t5 * qJD(5) - t45, 0, t20 * qJD(4), t38, -t40; -qJD(5) * t50 - t44, 0, -t38, 0, -t37; t5 * qJD(3) + qJD(4) * t50 - t42, 0, t40, t37, 0;];
Cq = t4;
