% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP3
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (216->36), mult. (419->44), div. (0->0), fcn. (284->4), ass. (0->21)
t26 = cos(qJ(3));
t41 = t26 ^ 2;
t40 = m(5) * pkin(3);
t25 = sin(qJ(3));
t38 = pkin(3) * t25;
t37 = Ifges(5,4) + Ifges(4,4);
t36 = t25 * mrSges(5,1) + t26 * mrSges(5,2);
t21 = sin(pkin(6)) * pkin(1) + pkin(5);
t34 = qJ(4) + t21;
t10 = t34 * t26;
t9 = t34 * t25;
t3 = m(5) * (t10 * t26 + t9 * t25) + (t25 ^ 2 + t41) * mrSges(5,3);
t35 = t3 * qJD(1);
t11 = -m(5) * t38 - t36;
t33 = t11 * qJD(1);
t30 = -cos(pkin(6)) * pkin(1) - pkin(2);
t29 = t25 * mrSges(4,1) + t26 * mrSges(4,2);
t27 = t26 * pkin(3) - t30;
t1 = t27 * t36 - t30 * t29 - t37 * t41 + (-mrSges(5,2) * t38 + t37 * t25 + t27 * t40) * t25 + (mrSges(5,1) * t38 + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,2)) * t25) * t26;
t28 = t1 * qJD(1);
t2 = [-t1 * qJD(3) + t3 * qJD(4), 0, -t28 + (t9 * mrSges(5,2) + (mrSges(4,2) * t21 - Ifges(4,6) - Ifges(5,6)) * t25 + (-mrSges(4,1) * t21 - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t26 + (-mrSges(5,1) - t40) * t10) * qJD(3), t35; 0, 0, (-t29 + t11) * qJD(3), 0; t11 * qJD(4) + t28, 0, 0, t33; -t11 * qJD(3) - t35, 0, -t33, 0;];
Cq = t2;
