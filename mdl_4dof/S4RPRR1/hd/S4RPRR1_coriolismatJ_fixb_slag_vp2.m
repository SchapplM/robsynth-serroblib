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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:34
% EndTime: 2018-11-14 13:50:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (572->35), mult. (1112->52), div. (0->0), fcn. (894->6), ass. (0->32)
t43 = -pkin(3) / 0.2e1;
t42 = pkin(1) * sin(pkin(7));
t25 = cos(pkin(7)) * pkin(1) + pkin(2);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t18 = t30 * t25 - t28 * t42;
t29 = cos(qJ(4));
t19 = t28 * t25 + t30 * t42;
t27 = sin(qJ(4));
t40 = t27 * t19;
t15 = t29 * t18 - t40;
t41 = t15 * mrSges(5,2);
t39 = t29 * t19;
t17 = pkin(3) + t18;
t13 = t29 * t17 - t40;
t14 = t27 * t17 + t39;
t34 = -t27 * t18 - t39;
t12 = t34 * mrSges(5,1);
t31 = -t19 * mrSges(4,1) - t18 * mrSges(4,2) + t12 - t41;
t1 = -m(5) * (t13 * t34 + t14 * t15) - t31;
t38 = t1 * qJD(1);
t4 = -t14 * mrSges(5,1) - t13 * mrSges(5,2);
t37 = t4 * qJD(1);
t36 = t4 * qJD(4);
t35 = t29 * t43 - t13 / 0.2e1;
t32 = (t27 * t43 - t14 / 0.2e1) * mrSges(5,1);
t2 = -t12 / 0.2e1 + t32 + (t15 / 0.2e1 + t35) * mrSges(5,2);
t21 = (t27 * mrSges(5,1) + t29 * mrSges(5,2)) * pkin(3);
t33 = -t2 * qJD(1) + t21 * qJD(3);
t20 = t21 * qJD(4);
t3 = -t41 / 0.2e1 + t12 / 0.2e1 + t35 * mrSges(5,2) + t32;
t5 = [-t1 * qJD(3) + t36, 0, -t38 + (m(5) * (t27 * t15 + t29 * t34) * pkin(3) + t31) * qJD(3) + t3 * qJD(4), t3 * qJD(3) + t36 + t37; 0, 0, 0, 0; t2 * qJD(4) + t38, 0, -t20, -t20 - t33; -t2 * qJD(3) - t37, 0, t33, 0;];
Cq  = t5;
