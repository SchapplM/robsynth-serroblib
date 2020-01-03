% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:52
% EndTime: 2019-12-31 16:28:53
% DurationCPUTime: 0.26s
% Computational Cost: add. (278->63), mult. (737->80), div. (0->0), fcn. (468->4), ass. (0->31)
t35 = sin(qJ(3));
t33 = t35 ^ 2;
t37 = cos(qJ(3));
t34 = t37 ^ 2;
t50 = t33 + t34;
t55 = t37 * mrSges(5,2);
t51 = -qJ(4) - pkin(5);
t25 = t51 * t35;
t26 = t51 * t37;
t59 = m(5) * (-t25 * t35 - t26 * t37);
t58 = m(5) * pkin(3);
t45 = mrSges(5,1) + t58;
t36 = sin(qJ(2));
t56 = m(5) * t36;
t52 = Ifges(5,4) + Ifges(4,4);
t38 = cos(qJ(2));
t4 = 0.4e1 * (m(4) / 0.4e1 + m(5) / 0.4e1) * (-0.1e1 + t50) * t38 * t36;
t49 = t4 * qJD(1);
t13 = t35 * t45 + t55;
t48 = t13 * qJD(2);
t44 = t37 * pkin(3) + pkin(2);
t43 = -t37 * mrSges(5,1) + t35 * mrSges(5,2);
t42 = t35 * mrSges(5,1) + t55;
t1 = t44 * t42 + (pkin(2) * mrSges(4,1) - pkin(3) * t43 + t35 * t52 + t44 * t58) * t35 + (pkin(2) * mrSges(4,2) + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,2)) * t35 - t52 * t37) * t37;
t40 = t1 * qJD(2);
t10 = (t33 / 0.2e1 + t34 / 0.2e1 - 0.1e1 / 0.2e1) * t56;
t7 = mrSges(5,3) * t50 + t59;
t39 = t10 * qJD(1) + t7 * qJD(2);
t9 = (t50 + 0.1e1) * t56 / 0.2e1;
t3 = (-t42 / 0.2e1 - t55 / 0.2e1 - t37 * mrSges(4,2) + (-t58 - mrSges(5,1) / 0.2e1 - mrSges(4,1)) * t35) * t38;
t2 = [t4 * qJD(2), t3 * qJD(3) + t9 * qJD(4) + t49 + ((-m(4) * pkin(2) - m(5) * t44 - t37 * mrSges(4,1) + t35 * mrSges(4,2) - mrSges(3,1) + t43) * t36 + (-mrSges(3,2) + t59 + (m(4) * pkin(5) + mrSges(4,3) + mrSges(5,3)) * t50) * t38) * qJD(2), t3 * qJD(2) + ((mrSges(4,2) + mrSges(5,2)) * t35 + (-mrSges(4,1) - t45) * t37) * qJD(3) * t36, t9 * qJD(2); t10 * qJD(4) - t49, -t1 * qJD(3) + t7 * qJD(4), -t40 + (-t25 * mrSges(5,2) + (mrSges(4,2) * pkin(5) - Ifges(4,6) - Ifges(5,6)) * t35 + (-mrSges(4,1) * pkin(5) - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t37 + t45 * t26) * qJD(3), t39; 0, -t13 * qJD(4) + t40, 0, -t48; -t10 * qJD(2), t13 * qJD(3) - t39, t48, 0;];
Cq = t2;
