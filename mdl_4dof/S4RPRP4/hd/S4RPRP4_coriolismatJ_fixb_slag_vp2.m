% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:39
% EndTime: 2019-12-31 16:43:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (142->40), mult. (254->48), div. (0->0), fcn. (176->4), ass. (0->19)
t22 = Ifges(5,5) - Ifges(4,4);
t13 = sin(qJ(3));
t14 = cos(qJ(3));
t16 = -cos(pkin(6)) * pkin(1) - pkin(2);
t15 = -t14 * pkin(3) - t13 * qJ(4);
t4 = t15 + t16;
t5 = -t14 * mrSges(5,1) - t13 * mrSges(5,3);
t17 = m(5) * t4 + t5;
t6 = -t13 * pkin(3) + t14 * qJ(4);
t1 = -t17 * t6 + (t16 * mrSges(4,2) - t4 * mrSges(5,3) - t22 * t14) * t14 + (t4 * mrSges(5,1) + t16 * mrSges(4,1) + (Ifges(5,1) - Ifges(5,3) + Ifges(4,1) - Ifges(4,2)) * t14 + t22 * t13) * t13;
t21 = t1 * qJD(1);
t2 = t17 * t13;
t20 = t2 * qJD(1);
t8 = m(5) * qJ(4) + mrSges(5,3);
t19 = t8 * qJD(3);
t18 = qJD(3) * t13;
t7 = sin(pkin(6)) * pkin(1) + pkin(5);
t3 = (m(5) * t7 + mrSges(5,2)) * t14;
t9 = [t1 * qJD(3) - t2 * qJD(4), 0, t21 + t3 * qJD(4) + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t18 + ((-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t14 + (m(5) * t15 - t14 * mrSges(4,1) + t13 * mrSges(4,2) + t5) * t7) * qJD(3), t3 * qJD(3) - t20; 0, 0, (m(5) * t6 + (-mrSges(4,2) + mrSges(5,3)) * t14) * qJD(3) + ((-mrSges(4,1) - mrSges(5,1)) * qJD(3) + m(5) * qJD(4)) * t13, m(5) * t18; -t21, 0, t8 * qJD(4), t19; t20, 0, -t19, 0;];
Cq = t9;
