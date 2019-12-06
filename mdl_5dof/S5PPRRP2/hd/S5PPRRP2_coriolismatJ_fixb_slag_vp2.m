% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:35
% EndTime: 2019-12-05 15:08:36
% DurationCPUTime: 0.30s
% Computational Cost: add. (385->63), mult. (915->93), div. (0->0), fcn. (818->6), ass. (0->41)
t37 = cos(qJ(4));
t69 = (mrSges(5,2) - mrSges(6,3)) * t37;
t68 = mrSges(5,1) + mrSges(6,1);
t35 = sin(qJ(4));
t52 = -t35 ^ 2 - t37 ^ 2;
t66 = Ifges(5,4) - Ifges(6,5);
t25 = -t37 * mrSges(6,1) - t35 * mrSges(6,3);
t65 = -t37 * mrSges(5,1) + t35 * mrSges(5,2) + t25;
t64 = m(6) / 0.2e1;
t34 = sin(pkin(8));
t36 = sin(qJ(3));
t47 = cos(pkin(8));
t58 = cos(qJ(3));
t21 = t36 * t34 - t58 * t47;
t62 = m(6) * t21;
t49 = t37 * qJ(5);
t59 = t35 * pkin(4);
t40 = -t49 + t59;
t61 = m(6) * t40;
t53 = t52 * pkin(6) * t21;
t51 = m(6) * qJD(5);
t22 = t58 * t34 + t36 * t47;
t1 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (0.1e1 + t52) * t22 * t21;
t50 = t1 * qJD(1);
t41 = -t37 * pkin(4) - t35 * qJ(5);
t23 = -pkin(3) + t41;
t43 = m(6) * t23 + t25;
t9 = t43 * t35;
t48 = t9 * qJD(3);
t46 = qJD(4) * t35;
t45 = qJD(4) * t37;
t28 = m(6) * qJ(5) + mrSges(6,3);
t44 = t28 * qJD(4);
t2 = (t59 / 0.2e1 - t49 / 0.2e1 - t40 / 0.2e1) * t62;
t4 = t43 * t40 + (-pkin(3) * mrSges(5,2) - t23 * mrSges(6,3) + t66 * t37) * t37 + (-pkin(3) * mrSges(5,1) + t23 * mrSges(6,1) + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t37 - t66 * t35) * t35;
t39 = t2 * qJD(1) - t4 * qJD(3);
t38 = (m(6) * t41 + t65) * qJD(4);
t24 = (m(6) * pkin(6) + mrSges(6,2)) * t37;
t8 = t35 * t62;
t3 = (t61 / 0.2e1 + t40 * t64 + t69 + t68 * t35) * t21;
t5 = [t1 * qJD(3), 0, t3 * qJD(4) - t8 * qJD(5) + t50 + ((-mrSges(4,1) + t65) * t22 + m(5) * (-pkin(3) * t22 + t53) + 0.2e1 * (t23 * t22 + t53) * t64 + (mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t52) * t21) * qJD(3), t3 * qJD(3) + (t37 * t51 + t38) * t22, m(6) * t22 * t45 - t8 * qJD(3); 0, 0, 0, (-t61 - t69) * qJD(4) + (-t68 * qJD(4) + t51) * t35, m(6) * t46; -t2 * qJD(4) - t50, 0, t4 * qJD(4) - t9 * qJD(5), t24 * qJD(5) + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t45 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t46 + pkin(6) * t38 - t39, t24 * qJD(4) - t48; t2 * qJD(3), 0, t39, t28 * qJD(5), t44; 0, 0, t48, -t44, 0;];
Cq = t5;
