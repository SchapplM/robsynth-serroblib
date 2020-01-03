% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:06
% EndTime: 2019-12-31 16:49:07
% DurationCPUTime: 0.33s
% Computational Cost: add. (1181->51), mult. (2249->73), div. (0->0), fcn. (2079->6), ass. (0->29)
t50 = sin(qJ(3));
t44 = sin(pkin(7)) * pkin(1) + pkin(5);
t71 = pkin(6) + t44;
t33 = t71 * t50;
t52 = cos(qJ(3));
t34 = t71 * t52;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t26 = -t49 * t33 + t51 * t34;
t40 = -t49 * t50 + t51 * t52;
t41 = -t49 * t52 - t51 * t50;
t57 = -t51 * t33 - t49 * t34;
t4 = -t26 * mrSges(5,1) - t57 * mrSges(5,2) + Ifges(5,5) * t40 + Ifges(5,6) * t41;
t28 = -t41 * mrSges(5,1) + t40 * mrSges(5,2);
t77 = t28 * qJD(4);
t72 = t50 * pkin(3);
t59 = -cos(pkin(7)) * pkin(1) - pkin(2);
t43 = -t52 * pkin(3) + t59;
t66 = t43 * t28;
t65 = t50 * mrSges(4,1);
t56 = Ifges(5,4) * t40 + (-Ifges(5,1) + Ifges(5,2)) * t41;
t1 = m(5) * t43 * t72 + t66 + t59 * t65 - t50 ^ 2 * Ifges(4,4) + (-mrSges(5,2) * t72 - Ifges(5,4) * t41) * t41 + (t59 * mrSges(4,2) + Ifges(4,4) * t52 + (Ifges(4,1) - Ifges(4,2)) * t50) * t52 + (-mrSges(5,1) * t72 + t56) * t40;
t55 = t1 * qJD(1);
t2 = -Ifges(5,4) * t41 ^ 2 + t40 * t56 + t66;
t54 = t2 * qJD(1);
t42 = (mrSges(5,1) * t49 + mrSges(5,2) * t51) * pkin(3);
t53 = t42 * qJD(3);
t39 = t42 * qJD(4);
t3 = [t1 * qJD(3) + t2 * qJD(4), 0, t4 * qJD(4) + t55 + (Ifges(4,5) * t52 - Ifges(4,6) * t50 + (m(5) * (-t26 * t51 + t49 * t57) + (-t51 * t40 + t49 * t41) * mrSges(5,3)) * pkin(3) + (-mrSges(4,1) * t52 + mrSges(4,2) * t50) * t44 + t4) * qJD(3), t54 + (qJD(4) + qJD(3)) * t4; 0, 0, -t77 + (-t52 * mrSges(4,2) - t28 - t65 + (t40 * t49 + t41 * t51) * pkin(3) * m(5)) * qJD(3), -t28 * qJD(3) - t77; -t55, 0, -t39, -t39 - t53; -t54, 0, t53, 0;];
Cq = t3;
