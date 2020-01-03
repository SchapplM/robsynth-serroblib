% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:41
% DurationCPUTime: 0.15s
% Computational Cost: add. (143->32), mult. (473->59), div. (0->0), fcn. (274->4), ass. (0->24)
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t37 = t22 ^ 2 + t24 ^ 2;
t43 = (t37 * mrSges(5,3) - mrSges(4,2)) * qJD(3);
t42 = t22 * mrSges(5,1) + t24 * mrSges(5,2);
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t3 = (-0.1e1 + t37) * (-t23 ^ 2 + t25 ^ 2);
t40 = -t3 / 0.2e1;
t35 = m(5) * qJD(3);
t1 = (-pkin(3) * mrSges(5,2) + Ifges(5,4) * t24) * t24 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) * t22 + (Ifges(5,1) - Ifges(5,2)) * t24) * t22;
t34 = t1 * qJD(3);
t32 = t37 * t25;
t8 = (-t25 + t32) * t23;
t33 = t8 * t35;
t11 = -t24 * mrSges(5,1) + t22 * mrSges(5,2);
t31 = (-mrSges(4,1) + t11) * qJD(3);
t29 = qJD(1) * t40 - t8 * qJD(2);
t28 = t8 * qJD(1) + qJD(2) * t40;
t26 = t11 * qJD(4);
t7 = t42 * t25;
t6 = t42 * t23;
t2 = t3 * t35 / 0.2e1;
t4 = [-t33, t2, t6 * qJD(4) + t25 * t31 - t23 * t43 + ((-t37 * t23 * pkin(5) - pkin(3) * t25) * qJD(3) - t28) * m(5), t6 * qJD(3) + t25 * t26; t2, t33, -t7 * qJD(4) + t23 * t31 + t25 * t43 + ((-t23 * pkin(3) + pkin(5) * t32) * qJD(3) - t29) * m(5), -t7 * qJD(3) + t23 * t26; t28 * m(5), t29 * m(5), t1 * qJD(4), t34 + (Ifges(5,5) * t24 - Ifges(5,6) * t22 + t11 * pkin(5)) * qJD(4); 0, 0, -t34, 0;];
Cq = t4;
