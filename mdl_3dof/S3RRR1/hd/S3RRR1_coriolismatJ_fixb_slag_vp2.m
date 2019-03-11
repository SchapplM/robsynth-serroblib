% Calculate matrix of centrifugal and coriolis load on the joints for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [3x3]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S3RRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (214->32), mult. (554->52), div. (0->0), fcn. (354->4), ass. (0->29)
t36 = -pkin(2) / 0.2e1;
t23 = cos(qJ(2));
t35 = t23 * pkin(1);
t22 = cos(qJ(3));
t20 = sin(qJ(3));
t21 = sin(qJ(2));
t33 = t20 * t21;
t14 = (t22 * t23 - t33) * pkin(1);
t34 = t14 * mrSges(4,2);
t32 = t21 * t22;
t19 = pkin(2) + t35;
t11 = -pkin(1) * t33 + t22 * t19;
t12 = pkin(1) * t32 + t20 * t19;
t25 = (-t20 * t23 - t32) * pkin(1);
t13 = mrSges(4,1) * t25;
t24 = -t21 * pkin(1) * mrSges(3,1) - mrSges(3,2) * t35 + t13 - t34;
t3 = -m(4) * (t11 * t25 + t12 * t14) - t24;
t31 = t3 * qJD(1);
t4 = -t12 * mrSges(4,1) - t11 * mrSges(4,2);
t30 = t4 * qJD(1);
t29 = t4 * qJD(3);
t28 = t22 * t36 - t11 / 0.2e1;
t26 = (t20 * t36 - t12 / 0.2e1) * mrSges(4,1);
t1 = -t13 / 0.2e1 + t26 + (t14 / 0.2e1 + t28) * mrSges(4,2);
t16 = (t20 * mrSges(4,1) + t22 * mrSges(4,2)) * pkin(2);
t27 = -t1 * qJD(1) + t16 * qJD(2);
t15 = t16 * qJD(3);
t2 = -t34 / 0.2e1 + t13 / 0.2e1 + t28 * mrSges(4,2) + t26;
t5 = [-t3 * qJD(2) + t29, -t31 + (m(4) * (t20 * t14 + t22 * t25) * pkin(2) + t24) * qJD(2) + t2 * qJD(3), t2 * qJD(2) + t29 + t30; t1 * qJD(3) + t31, -t15, -t15 - t27; -t1 * qJD(2) - t30, t27, 0;];
Cq  = t5;
