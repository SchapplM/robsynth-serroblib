% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:38
% EndTime: 2019-12-31 16:40:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (416->48), mult. (916->74), div. (0->0), fcn. (818->4), ass. (0->25)
t33 = sin(pkin(6));
t34 = cos(pkin(6));
t35 = cos(qJ(4));
t45 = sin(qJ(4));
t22 = -t33 * t45 - t34 * t35;
t23 = t33 * t35 - t34 * t45;
t46 = -t45 * t22 / 0.2e1 + t35 * t23 / 0.2e1;
t43 = -pkin(5) + qJ(2);
t20 = t33 * qJ(3) + pkin(1) + (pkin(2) + pkin(3)) * t34;
t1 = (t20 * mrSges(5,1) - Ifges(5,4) * t23) * t23 + (t20 * mrSges(5,2) + Ifges(5,4) * t22 + (Ifges(5,1) - Ifges(5,2)) * t23) * t22;
t41 = t1 * qJD(1);
t25 = t43 * t33;
t26 = t43 * t34;
t12 = t35 * t25 - t45 * t26;
t13 = t45 * t25 + t35 * t26;
t2 = (-t22 ^ 2 - t23 ^ 2) * mrSges(5,3) + m(5) * (t12 * t23 - t13 * t22) + (mrSges(3,3) + mrSges(4,2) + 0.4e1 * (m(4) / 0.4e1 + m(3) / 0.4e1) * qJ(2)) * (t33 ^ 2 + t34 ^ 2);
t40 = t2 * qJD(1);
t5 = (-t22 * mrSges(5,1) + t23 * mrSges(5,2) + m(5) * t20 - m(4) * (-t34 * pkin(2) - pkin(1)) + t34 * mrSges(4,1) + (m(4) * qJ(3) + mrSges(4,3)) * t33) * t33;
t39 = t5 * qJD(1);
t8 = -t23 * mrSges(5,1) - t22 * mrSges(5,2);
t38 = t8 * qJD(1);
t9 = m(4) * t33 + (t33 / 0.2e1 + t46) * m(5);
t37 = t9 * qJD(1);
t10 = (-t33 / 0.2e1 + t46) * m(5);
t3 = [t2 * qJD(2) + t5 * qJD(3) + t1 * qJD(4), t10 * qJD(3) + t40, t10 * qJD(2) + t39, t41 + (-t13 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,5) * t22 - Ifges(5,6) * t23) * qJD(4); -t9 * qJD(3) + t8 * qJD(4) - t40, 0, -t37, t38; t9 * qJD(2) - t39, t37, 0, (-t45 * mrSges(5,1) - t35 * mrSges(5,2)) * qJD(4); -t8 * qJD(2) - t41, -t38, 0, 0;];
Cq = t3;
