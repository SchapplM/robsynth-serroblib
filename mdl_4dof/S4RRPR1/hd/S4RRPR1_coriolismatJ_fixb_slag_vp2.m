% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:36
% EndTime: 2019-01-31 13:16:36
% DurationCPUTime: 0.19s
% Computational Cost: add. (720->46), mult. (1504->71), div. (0->0), fcn. (1230->6), ass. (0->40)
t32 = sin(pkin(7));
t37 = cos(qJ(2));
t33 = cos(pkin(7));
t35 = sin(qJ(2));
t50 = t33 * t35;
t28 = (-t32 * t37 - t50) * pkin(1);
t51 = t32 * t35;
t29 = (t33 * t37 - t51) * pkin(1);
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t41 = t36 * t28 - t34 * t29;
t18 = t41 * mrSges(5,1);
t19 = t34 * t28 + t36 * t29;
t52 = t19 * mrSges(5,2);
t55 = t28 * mrSges(4,1) - t29 * mrSges(4,2) - (t35 * mrSges(3,1) + t37 * mrSges(3,2)) * pkin(1) + t18 - t52;
t53 = pkin(2) * t32;
t31 = t37 * pkin(1) + pkin(2);
t39 = -pkin(1) * t51 + t33 * t31;
t23 = pkin(3) + t39;
t25 = pkin(1) * t50 + t32 * t31;
t14 = t36 * t23 - t34 * t25;
t15 = t34 * t23 + t36 * t25;
t3 = -m(4) * (t25 * t29 + t39 * t28) - m(5) * (t14 * t41 + t15 * t19) - t55;
t47 = t3 * qJD(1);
t13 = t15 * mrSges(5,1);
t4 = t14 * mrSges(5,2) + t13;
t46 = t4 * qJD(1);
t45 = t4 * qJD(4);
t30 = t33 * pkin(2) + pkin(3);
t26 = t34 * t30 + t36 * t53;
t22 = t26 * mrSges(5,1);
t24 = t36 * t30 - t34 * t53;
t8 = t24 * mrSges(5,2) + t22;
t44 = t8 * qJD(4);
t43 = -t13 / 0.2e1 - t22 / 0.2e1;
t42 = -t14 / 0.2e1 - t24 / 0.2e1;
t1 = -t18 / 0.2e1 + (t19 / 0.2e1 + t42) * mrSges(5,2) + t43;
t40 = -t1 * qJD(1) + t8 * qJD(2);
t2 = -t52 / 0.2e1 + t18 / 0.2e1 + t42 * mrSges(5,2) + t43;
t5 = [-t3 * qJD(2) - t45, t2 * qJD(4) - t47 + (m(5) * (t26 * t19 + t24 * t41) + m(4) * (t28 * t33 + t29 * t32) * pkin(2) + t55) * qJD(2), 0, t2 * qJD(2) - t45 - t46; t1 * qJD(4) + t47, -t44, 0, -t40 - t44; 0, 0, 0, 0; -t1 * qJD(2) + t46, t40, 0, 0;];
Cq  = t5;
