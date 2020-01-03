% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:48
% EndTime: 2019-12-31 17:36:48
% DurationCPUTime: 0.22s
% Computational Cost: add. (609->48), mult. (1382->72), div. (0->0), fcn. (1270->4), ass. (0->27)
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t45 = cos(qJ(5));
t61 = sin(qJ(5));
t32 = -t43 * t61 - t44 * t45;
t33 = t43 * t45 - t44 * t61;
t19 = t33 * mrSges(6,1) + t32 * mrSges(6,2);
t66 = t19 * qJD(5);
t64 = -t61 * t32 / 0.2e1 + t45 * t33 / 0.2e1;
t63 = t33 ^ 2;
t60 = -pkin(6) + qJ(3);
t30 = t43 * qJ(4) + pkin(2) + (pkin(3) + pkin(4)) * t44;
t7 = (-t32 * mrSges(6,1) + t33 * mrSges(6,2) + t44 * mrSges(5,1) + m(6) * t30 - m(5) * (-t44 * pkin(3) - pkin(2)) + (m(5) * qJ(4) + mrSges(5,3)) * t43) * t43;
t56 = t7 * qJD(2);
t55 = t19 * qJD(2);
t13 = m(5) * t43 + (t43 / 0.2e1 + t64) * m(6);
t54 = t13 * qJD(2);
t1 = -t63 * Ifges(6,4) + t30 * t19 + (Ifges(6,4) * t32 + (Ifges(6,1) - Ifges(6,2)) * t33) * t32;
t48 = t1 * qJD(2);
t35 = t60 * t43;
t36 = t60 * t44;
t20 = t45 * t35 - t61 * t36;
t21 = t61 * t35 + t45 * t36;
t3 = (-t32 ^ 2 - t63) * mrSges(6,3) + m(6) * (t20 * t33 - t21 * t32) + (mrSges(4,3) + mrSges(5,2) + 0.4e1 * (m(4) / 0.4e1 + m(5) / 0.4e1) * qJ(3)) * (t43 ^ 2 + t44 ^ 2);
t47 = t3 * qJD(2);
t14 = (-t43 / 0.2e1 + t64) * m(6);
t2 = [0, 0, 0, 0, -t66; 0, t3 * qJD(3) + t7 * qJD(4) + t1 * qJD(5), t14 * qJD(4) + t47, t14 * qJD(3) + t56, (-t21 * mrSges(6,1) - t20 * mrSges(6,2) + Ifges(6,5) * t32 - Ifges(6,6) * t33) * qJD(5) + t48; 0, -t13 * qJD(4) - t47 - t66, 0, -t54, -t55; 0, t13 * qJD(3) - t56, t54, 0, (-t61 * mrSges(6,1) - t45 * mrSges(6,2)) * qJD(5); 0, qJD(3) * t19 - t48, t55, 0, 0;];
Cq = t2;
