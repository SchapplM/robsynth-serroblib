% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP1
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
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:37
% EndTime: 2019-03-08 18:29:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (181->26), mult. (374->33), div. (0->0), fcn. (228->4), ass. (0->19)
t15 = cos(pkin(6)) * pkin(1) + pkin(2);
t18 = sin(qJ(3));
t13 = t18 * t15;
t19 = cos(qJ(3));
t25 = pkin(1) * sin(pkin(6));
t14 = t19 * t25;
t12 = t14 + t13;
t22 = qJ(4) + t12;
t28 = m(5) * t22 + mrSges(5,3);
t29 = t28 * qJD(4);
t11 = t19 * t15 - t18 * t25;
t20 = (-mrSges(4,1) - mrSges(5,1)) * t12 + (-mrSges(4,2) + mrSges(5,3)) * t11;
t1 = -m(5) * (t22 * t11 + (-pkin(3) - t11) * t12) - t20;
t24 = t1 * qJD(1);
t23 = t28 * qJD(1);
t16 = m(5) * qJ(4) + mrSges(5,3);
t4 = -mrSges(5,3) + 0.2e1 * (t12 / 0.4e1 - t14 / 0.4e1 - t13 / 0.4e1 - qJ(4) / 0.2e1) * m(5);
t21 = t4 * qJD(1) - t16 * qJD(3);
t2 = [-t1 * qJD(3) + t29, 0, -t24 + (m(5) * (-pkin(3) * t12 + qJ(4) * t11) + t20) * qJD(3) + t29, qJD(3) * t28 + t23; 0, 0, 0, 0; -t4 * qJD(4) + t24, 0, t16 * qJD(4), -t21; t4 * qJD(3) - t23, 0, t21, 0;];
Cq  = t2;
