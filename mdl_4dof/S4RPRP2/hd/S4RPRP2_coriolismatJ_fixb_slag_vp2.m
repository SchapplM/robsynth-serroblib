% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:42
% EndTime: 2019-03-08 18:30:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (305->31), mult. (433->35), div. (0->0), fcn. (265->2), ass. (0->22)
t27 = mrSges(4,1) + mrSges(5,1);
t17 = sin(qJ(3));
t18 = cos(qJ(3));
t19 = -pkin(1) - pkin(2);
t13 = -t17 * qJ(2) + t18 * t19;
t12 = -pkin(3) + t13;
t32 = m(5) * (-t12 + t13);
t26 = mrSges(4,2) + mrSges(5,2);
t21 = t26 * t18;
t14 = t18 * qJ(2) + t17 * t19;
t31 = t26 * t13 + t27 * t14;
t30 = m(5) * pkin(3);
t28 = t14 * t18;
t1 = t14 * t32 + t31;
t24 = t1 * qJD(1);
t20 = t17 * t32 / 0.2e1;
t2 = t21 + (t30 / 0.2e1 + t27) * t17 + t20;
t23 = t2 * qJD(1);
t4 = m(3) * qJ(2) + mrSges(3,3) + t21 + t27 * t17 + m(4) * (-t13 * t17 + t28) + m(5) * (-t12 * t17 + t28);
t22 = t4 * qJD(1);
t3 = -t17 * t30 / 0.2e1 + t20;
t5 = [t4 * qJD(2) + t1 * qJD(3), t3 * qJD(3) + t22, t24 + t3 * qJD(2) + (-t14 * t30 - t31) * qJD(3), 0; t2 * qJD(3) - t22, 0, t23 + (-t21 + (-t27 - t30) * t17) * qJD(3), 0; -t2 * qJD(2) - t24, -t23, 0, 0; 0, 0, 0, 0;];
Cq  = t5;
