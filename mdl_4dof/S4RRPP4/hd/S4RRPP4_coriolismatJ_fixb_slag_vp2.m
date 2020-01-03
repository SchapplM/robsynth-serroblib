% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:54
% EndTime: 2019-12-31 16:58:54
% DurationCPUTime: 0.29s
% Computational Cost: add. (320->80), mult. (607->88), div. (0->0), fcn. (379->2), ass. (0->35)
t35 = cos(qJ(2));
t51 = t35 * qJ(3);
t63 = m(5) * t51;
t62 = pkin(5) - qJ(4);
t60 = pkin(2) + pkin(3);
t34 = sin(qJ(2));
t59 = m(5) * t34;
t58 = t35 * pkin(2);
t57 = m(5) * qJ(3);
t56 = t35 * mrSges(4,1);
t55 = -mrSges(4,2) + mrSges(5,3);
t52 = t34 * qJ(3);
t37 = t60 * t35 + pkin(1) + t52;
t40 = t34 * pkin(2) - t51;
t41 = -t52 - t58;
t42 = t35 * mrSges(5,1) + t34 * mrSges(5,2);
t43 = -t34 * mrSges(4,3) - t56;
t44 = -t34 * mrSges(5,1) + t35 * mrSges(5,2);
t46 = m(5) * t60;
t47 = Ifges(3,4) - Ifges(5,4) - Ifges(4,5);
t1 = -t40 * t43 - (-t60 * t34 + t51) * t42 - t37 * t44 + ((-t60 * t57 - t47) * t35 + (mrSges(3,2) - t57) * pkin(1)) * t35 + (t47 * t34 + (mrSges(3,1) + t46) * pkin(1) + (t60 ^ 2 * m(5) - Ifges(3,1) - Ifges(4,1) - Ifges(5,1) + Ifges(3,2) + Ifges(5,2) + Ifges(4,3)) * t35 + (t46 * t34 - t63) * qJ(3)) * t34 + (m(4) * t40 + t34 * mrSges(4,1) - t35 * mrSges(4,3)) * (pkin(1) - t41);
t54 = t1 * qJD(1);
t3 = -t37 * t59 + (-t42 + m(4) * (-pkin(1) - t58) - t56 + (-m(4) * qJ(3) - mrSges(4,3)) * t34) * t34;
t53 = t3 * qJD(1);
t4 = -t63 + 0.2e1 * (t60 / 0.4e1 + pkin(2) / 0.4e1 + pkin(3) / 0.4e1) * t59 - t44;
t50 = t4 * qJD(1);
t29 = t62 * t34;
t30 = t62 * t35;
t5 = m(5) * (-t29 * t34 - t30 * t35) + (t34 ^ 2 + t35 ^ 2) * mrSges(5,3);
t49 = t5 * qJD(1);
t31 = mrSges(5,2) + mrSges(4,3) + (m(4) + m(5)) * qJ(3);
t48 = t31 * qJD(2);
t45 = qJD(1) * t59;
t7 = m(5) * t30 + (m(4) * pkin(5) - t55) * t35;
t2 = [-t1 * qJD(2) - t3 * qJD(3) + t5 * qJD(4), t7 * qJD(3) - t54 + (m(5) * (-qJ(3) * t29 - t30 * t60) - t30 * mrSges(5,1) - t29 * mrSges(5,2) + (-pkin(2) * mrSges(4,2) + mrSges(5,3) * t60 + Ifges(4,4) + Ifges(3,5) - Ifges(5,5)) * t35 + (t55 * qJ(3) - Ifges(3,6) + Ifges(4,6) - Ifges(5,6)) * t34 + (m(4) * t41 - t35 * mrSges(3,1) + t34 * mrSges(3,2) + t43) * pkin(5)) * qJD(2), t7 * qJD(2) - t53, t49; t4 * qJD(4) + t54, t31 * qJD(3), t48, t50; -qJD(4) * t59 + t53, -t48, 0, -t45; -t4 * qJD(2) + qJD(3) * t59 - t49, -t50, t45, 0;];
Cq = t2;
