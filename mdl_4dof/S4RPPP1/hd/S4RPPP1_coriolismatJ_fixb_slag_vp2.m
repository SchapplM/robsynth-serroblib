% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:20
% EndTime: 2019-03-08 18:26:21
% DurationCPUTime: 0.22s
% Computational Cost: add. (203->56), mult. (580->78), div. (0->0), fcn. (487->4), ass. (0->33)
t18 = cos(pkin(4));
t16 = sin(pkin(4));
t17 = cos(pkin(6));
t34 = t16 * t17;
t38 = (mrSges(4,1) + mrSges(5,1)) * t34 + (mrSges(5,2) + mrSges(4,3)) * t18;
t37 = m(5) * t18;
t36 = pkin(1) * t17;
t15 = sin(pkin(6));
t35 = t16 * t15;
t32 = t18 * t15;
t30 = qJ(2) * t16;
t31 = pkin(1) * t32 + t17 * t30;
t12 = t15 * t30;
t23 = -pkin(2) - t36;
t19 = t12 + (-qJ(4) + t23) * t18;
t6 = -t18 * qJ(3) - t31;
t5 = pkin(3) * t34 - t6;
t9 = mrSges(5,1) * t35 - t18 * mrSges(5,3);
t1 = (-t18 * mrSges(3,2) + mrSges(3,3) * t34 + t38) * t34 + (t9 + (mrSges(4,1) + mrSges(3,3)) * t35 + (-mrSges(3,1) + mrSges(4,2)) * t18) * t35 + (m(4) * (-t6 * t17 + (t23 * t18 + t12) * t15) + m(5) * (t5 * t17 + (pkin(3) * t35 + t19) * t15) + m(3) * (t31 * t17 - (t18 * t36 - t12) * t15)) * t16;
t29 = t1 * qJD(1);
t22 = -qJ(3) * t15 - pkin(1);
t20 = -pkin(2) * t17 + t22;
t7 = (-mrSges(5,2) * t15 - mrSges(5,3) * t17) * t16;
t2 = -t7 * t35 + (-m(4) * t20 - mrSges(4,2) * t17 + mrSges(4,3) * t15 - m(5) * ((-pkin(2) - qJ(4)) * t17 + t22)) * t16 ^ 2 * t15 + (-m(4) * t6 + m(5) * t5 + t38) * t18;
t28 = t2 * qJD(1);
t3 = t7 * t34 - m(5) * (-t19 * t18 + (-pkin(3) * t32 - (-qJ(4) * t17 + t20) * t34) * t16) + t18 * t9;
t27 = t3 * qJD(1);
t8 = (m(4) + m(5)) * t35;
t26 = t8 * qJD(1);
t25 = m(5) * t34;
t24 = qJD(1) * t37;
t21 = qJD(1) * t25;
t4 = [t1 * qJD(2) + t2 * qJD(3) - t3 * qJD(4), t29, t28, -t27; -t8 * qJD(3) - qJD(4) * t25 - t29, 0, -t26, -t21; t8 * qJD(2) - qJD(4) * t37 - t28, t26, 0, -t24; t27 + (qJD(2) * t34 + qJD(3) * t18) * m(5), t21, t24, 0;];
Cq  = t4;
