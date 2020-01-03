% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (122->38), mult. (234->46), div. (0->0), fcn. (156->2), ass. (0->17)
t19 = Ifges(4,4) - Ifges(5,5);
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t13 = -t12 * pkin(3) - t11 * qJ(4);
t3 = -pkin(2) + t13;
t5 = -t12 * mrSges(5,1) - t11 * mrSges(5,3);
t14 = m(5) * t3 + t5;
t6 = -t11 * pkin(3) + t12 * qJ(4);
t1 = -t14 * t6 + (-pkin(2) * mrSges(4,2) - t3 * mrSges(5,3) + t19 * t12) * t12 + (-pkin(2) * mrSges(4,1) + t3 * mrSges(5,1) + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t12 - t19 * t11) * t11;
t18 = t1 * qJD(2);
t2 = t14 * t11;
t17 = t2 * qJD(2);
t7 = m(5) * qJ(4) + mrSges(5,3);
t16 = t7 * qJD(3);
t15 = qJD(3) * t11;
t4 = (m(5) * pkin(5) + mrSges(5,2)) * t12;
t8 = [0, 0, (m(5) * t6 + (-mrSges(4,2) + mrSges(5,3)) * t12) * qJD(3) + ((-mrSges(4,1) - mrSges(5,1)) * qJD(3) + m(5) * qJD(4)) * t11, m(5) * t15; 0, t1 * qJD(3) - t2 * qJD(4), t18 + t4 * qJD(4) + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t15 + ((-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t12 + (m(5) * t13 - t12 * mrSges(4,1) + t11 * mrSges(4,2) + t5) * pkin(5)) * qJD(3), t4 * qJD(3) - t17; 0, -t18, t7 * qJD(4), t16; 0, t17, -t16, 0;];
Cq = t8;
