% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:06
% EndTime: 2019-12-05 15:35:08
% DurationCPUTime: 0.35s
% Computational Cost: add. (417->67), mult. (955->92), div. (0->0), fcn. (856->6), ass. (0->42)
t38 = sin(qJ(4));
t39 = cos(qJ(4));
t78 = -t38 ^ 2 - t39 ^ 2;
t52 = t39 * qJ(5);
t64 = t38 * pkin(4);
t42 = -t52 + t64;
t77 = -m(6) * t42 - (mrSges(5,2) - mrSges(6,3)) * t39;
t76 = mrSges(5,1) + mrSges(6,1);
t25 = -t39 * mrSges(6,1) - t38 * mrSges(6,3);
t43 = t39 * pkin(4) + t38 * qJ(5);
t44 = -t39 * mrSges(5,1) + t38 * mrSges(5,2);
t73 = qJD(4) * (m(6) * t43 - t25 - t44);
t49 = sin(pkin(8));
t72 = t49 * pkin(2);
t50 = cos(pkin(8));
t71 = t50 * pkin(2);
t70 = Ifges(5,4) - Ifges(6,5);
t31 = -pkin(3) - t71;
t22 = -t43 + t31;
t45 = m(6) * t22 + t25;
t62 = sin(qJ(2));
t63 = cos(qJ(2));
t23 = t49 * t62 - t50 * t63;
t66 = m(6) * t23;
t54 = m(6) * qJD(5);
t24 = -t49 * t63 - t50 * t62;
t1 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (-0.1e1 - t78) * t24 * t23;
t53 = t1 * qJD(1);
t9 = t45 * t38;
t51 = t9 * qJD(2);
t48 = qJD(4) * t38;
t47 = qJD(4) * t39;
t32 = m(6) * qJ(5) + mrSges(6,3);
t46 = t32 * qJD(4);
t2 = t45 * t42 + (t31 * mrSges(5,2) - t22 * mrSges(6,3) + t70 * t39) * t39 + (t31 * mrSges(5,1) + t22 * mrSges(6,1) + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t39 - t70 * t38) * t38;
t3 = (t64 / 0.2e1 - t52 / 0.2e1 - t42 / 0.2e1) * t66;
t41 = t3 * qJD(1) - t2 * qJD(2);
t30 = pkin(6) + t72;
t21 = (m(6) * t30 + mrSges(6,2)) * t39;
t8 = t38 * t66;
t4 = (t76 * t38 - t77) * t23;
t5 = [t1 * qJD(2), t4 * qJD(4) - t8 * qJD(5) + t53 + (-t62 * mrSges(3,1) - t63 * mrSges(3,2) + (m(4) * t71 - m(5) * t31 + mrSges(4,1) - t44 - t45) * t24 + (-m(4) * t72 + mrSges(4,2) + (mrSges(5,3) + mrSges(6,2) + (m(5) + m(6)) * t30) * t78) * t23) * qJD(2), 0, t4 * qJD(2) + (-t39 * t54 + t73) * t24, -m(6) * t24 * t47 - t8 * qJD(2); -t3 * qJD(4) - t53, t2 * qJD(4) - t9 * qJD(5), 0, t21 * qJD(5) + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t47 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t48 - t30 * t73 - t41, t21 * qJD(4) - t51; 0, 0, 0, t77 * qJD(4) + (-t76 * qJD(4) + t54) * t38, m(6) * t48; t3 * qJD(2), t41, 0, t32 * qJD(5), t46; 0, t51, 0, -t46, 0;];
Cq = t5;
