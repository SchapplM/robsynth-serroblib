% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.10s
% Computational Cost: add. (214->32), mult. (554->51), div. (0->0), fcn. (354->4), ass. (0->31)
t39 = -pkin(2) / 0.2e1;
t25 = cos(qJ(3));
t38 = t25 * pkin(1);
t24 = cos(qJ(4));
t22 = sin(qJ(4));
t23 = sin(qJ(3));
t35 = t22 * t23;
t16 = (t24 * t25 - t35) * pkin(1);
t37 = t16 * mrSges(5,2);
t36 = t22 * mrSges(5,1);
t34 = t23 * t24;
t21 = pkin(2) + t38;
t13 = -pkin(1) * t35 + t24 * t21;
t14 = pkin(1) * t34 + t22 * t21;
t27 = (-t22 * t25 - t34) * pkin(1);
t15 = mrSges(5,1) * t27;
t26 = -t23 * pkin(1) * mrSges(4,1) - mrSges(4,2) * t38 + t15 - t37;
t3 = -m(5) * (t13 * t27 + t14 * t16) - t26;
t33 = t3 * qJD(2);
t10 = t14 * mrSges(5,1);
t4 = t13 * mrSges(5,2) + t10;
t32 = t4 * qJD(2);
t31 = t4 * qJD(4);
t30 = t24 * t39 - t13 / 0.2e1;
t29 = -t10 / 0.2e1 + t36 * t39;
t1 = -t15 / 0.2e1 + (t16 / 0.2e1 + t30) * mrSges(5,2) + t29;
t18 = (t24 * mrSges(5,2) + t36) * pkin(2);
t28 = -t1 * qJD(2) + t18 * qJD(3);
t17 = t18 * qJD(4);
t2 = -t37 / 0.2e1 + t15 / 0.2e1 + t30 * mrSges(5,2) + t29;
t5 = [0, 0, 0, 0; 0, -t3 * qJD(3) - t31, -t33 + (m(5) * (t22 * t16 + t24 * t27) * pkin(2) + t26) * qJD(3) + t2 * qJD(4), t2 * qJD(3) - t31 - t32; 0, t1 * qJD(4) + t33, -t17, -t17 - t28; 0, -t1 * qJD(3) + t32, t28, 0;];
Cq  = t5;
