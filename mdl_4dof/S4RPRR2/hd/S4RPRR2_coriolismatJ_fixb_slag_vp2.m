% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:08
% EndTime: 2019-12-31 16:48:09
% DurationCPUTime: 0.20s
% Computational Cost: add. (359->42), mult. (741->59), div. (0->0), fcn. (530->6), ass. (0->35)
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t53 = -Ifges(5,4) * t25 + (Ifges(5,1) - Ifges(5,2)) * t27;
t54 = t53 * t25;
t20 = cos(pkin(7)) * pkin(1) + pkin(2);
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t45 = pkin(1) * sin(pkin(7));
t11 = t28 * t20 - t26 * t45;
t12 = t26 * t20 + t28 * t45;
t18 = -t27 * mrSges(5,1) + t25 * mrSges(5,2);
t23 = t27 ^ 2;
t36 = (t25 ^ 2 + t23) * t11;
t50 = -t11 * mrSges(4,2) + (-mrSges(4,1) + t18) * t12 + mrSges(5,3) * t36;
t47 = -t11 / 0.2e1;
t39 = t27 * mrSges(5,2);
t40 = t25 * mrSges(5,1);
t19 = -t39 - t40;
t44 = pkin(3) * t19;
t42 = Ifges(5,4) * t27;
t10 = pkin(6) + t12;
t31 = -pkin(3) - t11;
t1 = m(5) * (t10 * t36 + t31 * t12) + t50;
t38 = t1 * qJD(1);
t29 = -Ifges(5,4) * t23 - t54;
t7 = t31 * t19;
t4 = t7 + t29;
t37 = t4 * qJD(1);
t35 = Ifges(5,5) * t27 - Ifges(5,6) * t25;
t34 = t7 / 0.2e1 - t44 / 0.2e1;
t2 = (mrSges(5,2) * t47 - t42) * t27 + (mrSges(5,1) * t47 - t53) * t25 + t34;
t5 = t29 - t44;
t33 = t2 * qJD(1) + t5 * qJD(3);
t3 = t42 * t27 + (-t39 / 0.2e1 - t40 / 0.2e1) * t11 - t34 + t54;
t6 = [t1 * qJD(3) - t4 * qJD(4), 0, t38 + (m(5) * (-pkin(3) * t12 + pkin(6) * t36) + t50) * qJD(3) + t3 * qJD(4), -t37 + t3 * qJD(3) + (t18 * t10 + t35) * qJD(4); 0, 0, 0, t19 * qJD(4); -t2 * qJD(4) - t38, 0, -t5 * qJD(4), (t18 * pkin(6) + t35) * qJD(4) - t33; t2 * qJD(3) + t37, 0, t33, 0;];
Cq = t6;
