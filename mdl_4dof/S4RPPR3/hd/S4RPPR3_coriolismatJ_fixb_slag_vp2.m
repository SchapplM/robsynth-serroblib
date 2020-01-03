% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:51
% EndTime: 2019-12-31 16:37:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (505->28), mult. (990->42), div. (0->0), fcn. (964->6), ass. (0->21)
t20 = sin(pkin(7));
t21 = cos(pkin(7));
t22 = sin(qJ(4));
t23 = cos(qJ(4));
t16 = t23 * t20 + t22 * t21;
t33 = t16 ^ 2;
t17 = sin(pkin(6)) * pkin(1) + qJ(3);
t31 = pkin(5) + t17;
t15 = -t22 * t20 + t23 * t21;
t11 = t31 * t20;
t12 = t31 * t21;
t8 = -t23 * t11 - t22 * t12;
t9 = -t22 * t11 + t23 * t12;
t4 = (t15 ^ 2 + t33) * mrSges(5,3) + m(5) * (t9 * t15 - t8 * t16) + (m(4) * t17 + mrSges(4,3)) * (t20 ^ 2 + t21 ^ 2);
t29 = t4 * qJD(1);
t10 = t16 * mrSges(5,1) + t15 * mrSges(5,2);
t28 = t10 * qJD(1);
t27 = t10 * qJD(4);
t1 = -t33 * Ifges(5,4) + (-cos(pkin(6)) * pkin(1) - pkin(2) - t21 * pkin(3)) * t10 + (Ifges(5,4) * t15 + (Ifges(5,1) - Ifges(5,2)) * t16) * t15;
t24 = t1 * qJD(1);
t2 = [t4 * qJD(3) + t1 * qJD(4), 0, t29, (-t9 * mrSges(5,1) - t8 * mrSges(5,2) + Ifges(5,5) * t15 - Ifges(5,6) * t16) * qJD(4) + t24; 0, 0, 0, -t27; t27 - t29, 0, 0, t28; -t10 * qJD(3) - t24, 0, -t28, 0;];
Cq = t2;
