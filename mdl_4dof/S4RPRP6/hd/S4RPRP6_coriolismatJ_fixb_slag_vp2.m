% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:00
% EndTime: 2019-12-31 16:46:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (250->52), mult. (453->61), div. (0->0), fcn. (261->2), ass. (0->31)
t50 = -m(5) * pkin(3) - mrSges(5,1);
t31 = -pkin(1) - pkin(5);
t49 = -qJ(4) + t31;
t29 = sin(qJ(3));
t27 = t29 ^ 2;
t30 = cos(qJ(3));
t47 = pkin(3) * t30;
t46 = t29 * mrSges(5,1);
t45 = mrSges(5,2) + mrSges(4,2);
t44 = Ifges(5,4) + Ifges(4,4);
t28 = t30 ^ 2;
t43 = t27 + t28;
t34 = t29 * pkin(3) + qJ(2);
t33 = m(5) * t34;
t35 = t30 * mrSges(5,1) - t29 * mrSges(5,2);
t1 = -t44 * t27 - t34 * t35 - t46 * t47 + (-qJ(2) * mrSges(4,1) - mrSges(5,2) * t47 - pkin(3) * t33 + t44 * t30) * t30 + (qJ(2) * mrSges(4,2) + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,2)) * t30) * t29;
t42 = t1 * qJD(1);
t22 = t49 * t29;
t23 = t49 * t30;
t3 = m(5) * (-t29 * t22 - t30 * t23) + t43 * mrSges(5,3);
t41 = t3 * qJD(1);
t5 = t29 * mrSges(4,1) + t45 * t30 + mrSges(3,3) + t33 + t46 + (m(3) + m(4)) * qJ(2);
t40 = t5 * qJD(1);
t8 = -m(5) * t47 - t35;
t39 = t8 * qJD(1);
t38 = qJD(3) * t29;
t37 = qJD(3) * t30;
t18 = (-t27 / 0.2e1 - t28 / 0.2e1 - 0.1e1 / 0.2e1) * m(5);
t36 = t18 * qJD(1);
t17 = (-t43 / 0.2e1 + 0.1e1 / 0.2e1) * m(5);
t2 = [t5 * qJD(2) - t1 * qJD(3) + t3 * qJD(4), t17 * qJD(4) + t40, -t42 + (-t23 * mrSges(5,2) + t50 * t22) * qJD(3) + (-mrSges(4,2) * t31 - Ifges(4,6) - Ifges(5,6)) * t37 + (-mrSges(4,1) * t31 + mrSges(5,3) * pkin(3) - Ifges(4,5) - Ifges(5,5)) * t38, t17 * qJD(2) + t41; t18 * qJD(4) - t40, 0, -t45 * t37 + (-mrSges(4,1) + t50) * t38, t36; t8 * qJD(4) + t42, 0, 0, t39; -t18 * qJD(2) - t8 * qJD(3) - t41, -t36, -t39, 0;];
Cq = t2;
