% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:45
% DurationCPUTime: 0.12s
% Computational Cost: add. (166->35), mult. (369->43), div. (0->0), fcn. (234->2), ass. (0->17)
t36 = m(5) * pkin(3);
t23 = sin(qJ(3));
t34 = pkin(3) * t23;
t33 = Ifges(5,4) + Ifges(4,4);
t32 = -qJ(4) - pkin(5);
t24 = cos(qJ(3));
t31 = t23 * mrSges(5,1) + t24 * mrSges(5,2);
t18 = t32 * t23;
t19 = t32 * t24;
t3 = m(5) * (-t18 * t23 - t19 * t24) + (t23 ^ 2 + t24 ^ 2) * mrSges(5,3);
t30 = t3 * qJD(2);
t7 = -m(5) * t34 - t31;
t29 = t7 * qJD(2);
t26 = t24 * pkin(3) + pkin(2);
t1 = t26 * t31 + (pkin(2) * mrSges(4,1) - mrSges(5,2) * t34 + t33 * t23 + t26 * t36) * t23 + (mrSges(5,1) * t34 + pkin(2) * mrSges(4,2) + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,2)) * t23 - t33 * t24) * t24;
t25 = t1 * qJD(2);
t2 = [0, 0, (-t23 * mrSges(4,1) - t24 * mrSges(4,2) + t7) * qJD(3), 0; 0, -t1 * qJD(3) + t3 * qJD(4), -t25 + (-t18 * mrSges(5,2) + (mrSges(4,2) * pkin(5) - Ifges(4,6) - Ifges(5,6)) * t23 + (-mrSges(4,1) * pkin(5) - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t24 - (-mrSges(5,1) - t36) * t19) * qJD(3), t30; 0, t7 * qJD(4) + t25, 0, t29; 0, -t7 * qJD(3) - t30, -t29, 0;];
Cq = t2;
