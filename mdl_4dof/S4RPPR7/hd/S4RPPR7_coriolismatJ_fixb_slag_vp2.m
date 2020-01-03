% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR7
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (595->43), mult. (1061->59), div. (0->0), fcn. (996->4), ass. (0->33)
t28 = sin(pkin(6));
t29 = cos(pkin(6));
t31 = cos(qJ(4));
t45 = sin(qJ(4));
t20 = -t45 * t28 + t31 * t29;
t48 = t20 ^ 2;
t47 = -m(5) / 0.2e1;
t30 = -pkin(1) - qJ(3);
t46 = -pkin(5) + t30;
t44 = t28 ^ 2 + t29 ^ 2;
t19 = -t31 * t28 - t45 * t29;
t13 = t20 * mrSges(5,1) + t19 * mrSges(5,2);
t24 = t28 * pkin(3) + qJ(2);
t1 = -t48 * Ifges(5,4) + t24 * t13 + (Ifges(5,4) * t19 + (Ifges(5,1) - Ifges(5,2)) * t20) * t19;
t43 = t1 * qJD(1);
t22 = t46 * t28;
t23 = t46 * t29;
t11 = -t45 * t22 + t31 * t23;
t12 = t31 * t22 + t45 * t23;
t36 = m(4) * t44;
t37 = t19 ^ 2 + t48;
t5 = t37 * mrSges(5,3) + t44 * mrSges(4,3) + m(5) * (-t11 * t20 + t12 * t19) - t30 * t36;
t42 = t5 * qJD(1);
t34 = t37 * t47 - t36 / 0.2e1;
t38 = -m(4) / 0.2e1 + t47;
t7 = t34 + t38;
t41 = t7 * qJD(1);
t35 = t19 * mrSges(5,1) - t20 * mrSges(5,2);
t9 = m(5) * t24 + mrSges(3,3) + t29 * mrSges(4,2) + t28 * mrSges(4,1) + (m(4) + m(3)) * qJ(2) - t35;
t40 = t9 * qJD(1);
t39 = t13 * qJD(1);
t6 = t34 - t38;
t2 = [t9 * qJD(2) + t5 * qJD(3) + t1 * qJD(4), t6 * qJD(3) + t40, t6 * qJD(2) + t42, t43 + (-t12 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,5) * t19 - Ifges(5,6) * t20) * qJD(4); t7 * qJD(3) - t40, 0, t41, t35 * qJD(4); -t7 * qJD(2) + t13 * qJD(4) - t42, -t41, 0, t39; -t13 * qJD(3) - t43, 0, -t39, 0;];
Cq = t2;
