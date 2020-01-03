% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (71->24), mult. (204->44), div. (0->0), fcn. (132->4), ass. (0->20)
t14 = cos(qJ(4));
t11 = t14 ^ 2;
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t12 = sin(qJ(4));
t22 = t12 ^ 2 + t11;
t25 = m(5) * (-0.1e1 + t22) * t15 * t13;
t24 = t12 * mrSges(5,1);
t23 = t14 * mrSges(5,2);
t21 = qJD(2) * t25;
t20 = t22 * t15;
t19 = -t14 * mrSges(5,1) + t12 * mrSges(5,2);
t6 = t23 + t24;
t1 = -pkin(3) * t6 + t11 * Ifges(5,4) + (-Ifges(5,4) * t12 + (Ifges(5,1) - Ifges(5,2)) * t14) * t12;
t17 = -t23 / 0.2e1 - t24 / 0.2e1;
t2 = (t6 / 0.2e1 + t17) * t15;
t18 = t2 * qJD(2) - t1 * qJD(3);
t16 = t13 * t19;
t3 = (-t6 / 0.2e1 + t17) * t15;
t4 = [0, 0, 0, t6 * qJD(4); 0, qJD(3) * t25, t21 + (-t15 * mrSges(4,2) - t13 * mrSges(4,1) + t16 + m(5) * (-pkin(3) * t13 + pkin(5) * t20) + mrSges(5,3) * t20) * qJD(3) + t3 * qJD(4), t3 * qJD(3) + qJD(4) * t16; 0, -t2 * qJD(4) - t21, t1 * qJD(4), (Ifges(5,5) * t14 - Ifges(5,6) * t12 + t19 * pkin(5)) * qJD(4) - t18; 0, t2 * qJD(3), t18, 0;];
Cq = t4;
