% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (150->34), mult. (287->51), div. (0->0), fcn. (190->4), ass. (0->21)
t25 = cos(qJ(4));
t20 = t25 ^ 2;
t24 = sin(qJ(4));
t32 = t24 ^ 2 + t20;
t36 = -pkin(1) - pkin(2);
t34 = t24 * mrSges(5,1);
t33 = t25 * mrSges(5,2);
t21 = sin(pkin(6));
t22 = cos(pkin(6));
t14 = t22 * qJ(2) + t21 * t36 - pkin(5);
t26 = t21 * qJ(2) - t22 * t36 + pkin(3);
t28 = t25 * mrSges(5,1) - t24 * mrSges(5,2);
t1 = -mrSges(3,3) + t32 * t22 * mrSges(5,3) + (-m(5) * t14 * t32 - mrSges(4,2)) * t22 + (-m(5) * t26 - mrSges(4,1) - t28) * t21 + (-m(4) * (t21 ^ 2 + t22 ^ 2) - m(3)) * qJ(2);
t31 = t1 * qJD(1);
t27 = -t33 - t34;
t2 = -t20 * Ifges(5,4) - t26 * t27 + (Ifges(5,4) * t24 + (-Ifges(5,1) + Ifges(5,2)) * t25) * t24;
t30 = t2 * qJD(1);
t5 = t22 * t27;
t29 = t5 * qJD(1);
t6 = -t5 / 0.2e1 + (-t33 / 0.2e1 - t34 / 0.2e1) * t22;
t3 = [-t1 * qJD(2) - t2 * qJD(4), -t31 + t6 * qJD(4) + m(5) * (-0.1e1 + t32) * t22 * qJD(2) * t21, 0, -t30 + t6 * qJD(2) + (-Ifges(5,5) * t25 + Ifges(5,6) * t24 - t28 * t14) * qJD(4); -t5 * qJD(4) + t31, 0, 0, -t28 * qJD(4) * t21 - t29; 0, 0, 0, t27 * qJD(4); t5 * qJD(2) + t30, t29, 0, 0;];
Cq = t3;
