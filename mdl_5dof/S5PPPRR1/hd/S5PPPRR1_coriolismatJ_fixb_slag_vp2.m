% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:54
% EndTime: 2019-12-05 14:57:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (619->48), mult. (1885->81), div. (0->0), fcn. (2022->8), ass. (0->40)
t46 = cos(qJ(5));
t60 = t46 * mrSges(6,2);
t44 = sin(qJ(5));
t62 = t44 * mrSges(6,1);
t36 = t60 + t62;
t49 = t60 / 0.2e1 + t62 / 0.2e1;
t67 = t49 + t36 / 0.2e1;
t42 = sin(pkin(9));
t45 = sin(qJ(4));
t57 = cos(pkin(9));
t64 = cos(qJ(4));
t66 = -t45 * t42 + t64 * t57;
t41 = t46 ^ 2;
t34 = t64 * t42 + t45 * t57;
t35 = -t46 * mrSges(6,1) + t44 * mrSges(6,2);
t63 = t34 * t35;
t59 = -t44 ^ 2 - t41;
t58 = cos(pkin(8));
t55 = t59 * t66;
t54 = (0.1e1 + t59) * t34;
t43 = sin(pkin(8));
t25 = t34 * t43;
t26 = t66 * t43;
t20 = -t44 * t26 - t58 * t46;
t21 = t46 * t26 - t58 * t44;
t50 = t20 * t44 - t21 * t46 + t26;
t2 = m(6) * (t25 * t54 - t50 * t66) / 0.2e1;
t4 = m(6) * t50 * t25;
t52 = t4 * qJD(1) + t2 * qJD(2);
t7 = m(6) * t66 * t54;
t51 = t2 * qJD(1) - t7 * qJD(2);
t48 = -t36 / 0.2e1 + t49;
t14 = -Ifges(6,4) * t41 + pkin(4) * t36 + (Ifges(6,4) * t44 + (-Ifges(6,1) + Ifges(6,2)) * t46) * t44;
t15 = t48 * t66;
t5 = t48 * t25;
t47 = t5 * qJD(1) - t15 * qJD(2) + t14 * qJD(4);
t16 = t67 * t66;
t6 = t67 * t25;
t1 = t2 * qJD(4);
t3 = [t4 * qJD(4), t1, 0, ((-m(6) * pkin(4) - mrSges(5,1) + t35) * t26 + (mrSges(5,2) + (m(6) * pkin(6) + mrSges(6,3)) * t59) * t25) * qJD(4) + t6 * qJD(5) + t52, t6 * qJD(4) + (-t21 * mrSges(6,1) - t20 * mrSges(6,2)) * qJD(5); t1, -t7 * qJD(4), 0, (-t66 * mrSges(5,2) - t34 * mrSges(5,1) + t63 + m(6) * (-pkin(4) * t34 - pkin(6) * t55) - mrSges(6,3) * t55) * qJD(4) - t16 * qJD(5) + t51, -t16 * qJD(4) + qJD(5) * t63; 0, 0, 0, 0, -t36 * qJD(5); -t5 * qJD(5) - t52, t15 * qJD(5) - t51, 0, -t14 * qJD(5), (Ifges(6,5) * t46 - Ifges(6,6) * t44 + t35 * pkin(6)) * qJD(5) - t47; t5 * qJD(4), -t15 * qJD(4), 0, t47, 0;];
Cq = t3;
