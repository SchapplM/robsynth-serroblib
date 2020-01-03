% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:35
% EndTime: 2019-12-31 16:31:36
% DurationCPUTime: 0.20s
% Computational Cost: add. (201->40), mult. (483->60), div. (0->0), fcn. (290->4), ass. (0->32)
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t50 = -Ifges(5,4) * t23 + (Ifges(5,1) - Ifges(5,2)) * t25;
t51 = t50 * t23;
t22 = t25 ^ 2;
t47 = t23 ^ 2 + t22;
t26 = cos(qJ(3));
t45 = t47 * t26;
t35 = t25 * mrSges(5,2);
t36 = t23 * mrSges(5,1);
t13 = -t35 - t36;
t43 = pkin(3) * t13;
t42 = t26 * pkin(2);
t40 = Ifges(5,4) * t25;
t17 = -pkin(3) - t42;
t39 = t17 * t13;
t24 = sin(qJ(3));
t16 = t24 * pkin(2) + pkin(6);
t12 = -t25 * mrSges(5,1) + t23 * mrSges(5,2);
t27 = (t47 * mrSges(5,3) - mrSges(4,2)) * t26 + (t12 - mrSges(4,1)) * t24;
t3 = (m(5) * (t45 * t16 + t17 * t24) + t27) * pkin(2);
t34 = t3 * qJD(2);
t31 = t40 * t25 + t51;
t4 = t31 - t39;
t33 = t4 * qJD(2);
t32 = -t42 / 0.2e1;
t30 = Ifges(5,5) * t25 - Ifges(5,6) * t23;
t1 = (mrSges(5,2) * t32 - t40) * t25 + (-pkin(3) / 0.2e1 + t17 / 0.2e1) * t13 + (mrSges(5,1) * t32 - t50) * t23;
t5 = -Ifges(5,4) * t22 - t43 - t51;
t29 = t1 * qJD(2) + t5 * qJD(3);
t2 = t43 / 0.2e1 - t39 / 0.2e1 + (-t35 / 0.2e1 - t36 / 0.2e1) * t42 + t31;
t6 = [0, 0, 0, t13 * qJD(4); 0, t3 * qJD(3) + t4 * qJD(4), t34 + t2 * qJD(4) + (m(5) * (-pkin(3) * t24 + t45 * pkin(6)) + t27) * qJD(3) * pkin(2), t33 + t2 * qJD(3) + (t12 * t16 + t30) * qJD(4); 0, -t1 * qJD(4) - t34, -t5 * qJD(4), (t12 * pkin(6) + t30) * qJD(4) - t29; 0, t1 * qJD(3) - t33, t29, 0;];
Cq = t6;
