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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:19
% EndTime: 2018-11-14 13:53:19
% DurationCPUTime: 0.17s
% Computational Cost: add. (720->46), mult. (1504->71), div. (0->0), fcn. (1230->6), ass. (0->40)
t30 = sin(pkin(7));
t35 = cos(qJ(2));
t31 = cos(pkin(7));
t33 = sin(qJ(2));
t48 = t31 * t33;
t26 = (-t30 * t35 - t48) * pkin(1);
t49 = t30 * t33;
t27 = (t31 * t35 - t49) * pkin(1);
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t39 = t34 * t26 - t32 * t27;
t16 = t39 * mrSges(5,1);
t17 = t32 * t26 + t34 * t27;
t50 = t17 * mrSges(5,2);
t54 = t26 * mrSges(4,1) - t27 * mrSges(4,2) - (t33 * mrSges(3,1) + t35 * mrSges(3,2)) * pkin(1) + t16 - t50;
t52 = pkin(2) * t30;
t29 = t35 * pkin(1) + pkin(2);
t37 = -pkin(1) * t49 + t31 * t29;
t21 = pkin(3) + t37;
t23 = pkin(1) * t48 + t30 * t29;
t13 = t32 * t21 + t34 * t23;
t51 = t13 * mrSges(5,1);
t12 = t34 * t21 - t32 * t23;
t3 = -m(5) * (t12 * t39 + t13 * t17) - m(4) * (t23 * t27 + t37 * t26) - t54;
t45 = t3 * qJD(1);
t4 = -t12 * mrSges(5,2) - t51;
t44 = t4 * qJD(1);
t43 = t4 * qJD(4);
t28 = t31 * pkin(2) + pkin(3);
t24 = t32 * t28 + t34 * t52;
t20 = t24 * mrSges(5,1);
t22 = t34 * t28 - t32 * t52;
t7 = t22 * mrSges(5,2) + t20;
t42 = t7 * qJD(4);
t41 = -t22 / 0.2e1 - t12 / 0.2e1;
t40 = -t20 / 0.2e1 - t51 / 0.2e1;
t1 = -t16 / 0.2e1 + (t17 / 0.2e1 + t41) * mrSges(5,2) + t40;
t38 = -t1 * qJD(1) + t7 * qJD(2);
t2 = -t50 / 0.2e1 + t16 / 0.2e1 + t41 * mrSges(5,2) + t40;
t5 = [-t3 * qJD(2) + t43, t2 * qJD(4) - t45 + (m(5) * (t24 * t17 + t22 * t39) + m(4) * (t26 * t31 + t27 * t30) * pkin(2) + t54) * qJD(2), 0, t2 * qJD(2) + t43 + t44; t1 * qJD(4) + t45, -t42, 0, -t38 - t42; 0, 0, 0, 0; -t1 * qJD(2) - t44, t38, 0, 0;];
Cq  = t5;
