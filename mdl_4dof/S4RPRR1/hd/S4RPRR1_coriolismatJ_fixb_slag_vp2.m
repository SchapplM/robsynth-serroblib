% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR1
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
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:40
% EndTime: 2019-03-08 18:31:40
% DurationCPUTime: 0.11s
% Computational Cost: add. (572->35), mult. (1112->53), div. (0->0), fcn. (894->6), ass. (0->32)
t40 = -pkin(3) / 0.2e1;
t39 = pkin(1) * sin(pkin(7));
t22 = cos(pkin(7)) * pkin(1) + pkin(2);
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t15 = t27 * t22 - t25 * t39;
t24 = sin(qJ(4));
t16 = t25 * t22 + t27 * t39;
t26 = cos(qJ(4));
t35 = t26 * t16;
t12 = -t24 * t15 - t35;
t38 = t12 * mrSges(5,1);
t36 = t24 * t16;
t13 = t26 * t15 - t36;
t37 = t13 * mrSges(5,2);
t14 = pkin(3) + t15;
t10 = t26 * t14 - t36;
t11 = t24 * t14 + t35;
t28 = -t16 * mrSges(4,1) - t15 * mrSges(4,2) - t37 + t38;
t1 = m(5) * (t10 * t12 + t11 * t13) + t28;
t34 = t1 * qJD(1);
t4 = -t11 * mrSges(5,1) - t10 * mrSges(5,2);
t33 = t4 * qJD(1);
t32 = t4 * qJD(4);
t31 = t24 * t40 - t11 / 0.2e1;
t30 = t26 * t40 - t10 / 0.2e1;
t18 = (t24 * mrSges(5,1) + t26 * mrSges(5,2)) * pkin(3);
t2 = (t13 / 0.2e1 + t30) * mrSges(5,2) + (-t12 / 0.2e1 + t31) * mrSges(5,1);
t29 = -t2 * qJD(1) + t18 * qJD(3);
t17 = t18 * qJD(4);
t3 = -t37 / 0.2e1 + t38 / 0.2e1 + t30 * mrSges(5,2) + t31 * mrSges(5,1);
t5 = [t1 * qJD(3) + t32, 0, t34 + (m(5) * (t12 * t26 + t13 * t24) * pkin(3) + t28) * qJD(3) + t3 * qJD(4), t3 * qJD(3) + t32 + t33; 0, 0, 0, 0; t2 * qJD(4) - t34, 0, -t17, -t17 - t29; -t2 * qJD(3) - t33, 0, t29, 0;];
Cq  = t5;
