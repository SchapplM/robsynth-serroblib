% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:06
% EndTime: 2019-12-31 16:55:07
% DurationCPUTime: 0.33s
% Computational Cost: add. (1274->54), mult. (2187->68), div. (0->0), fcn. (2026->4), ass. (0->30)
t52 = sin(qJ(3));
t55 = -pkin(1) - pkin(5);
t72 = t52 * t55;
t46 = -t52 * pkin(6) + t72;
t54 = cos(qJ(3));
t47 = (-pkin(6) + t55) * t54;
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t30 = t53 * t46 + t51 * t47;
t43 = t51 * t52 - t53 * t54;
t44 = -t51 * t54 - t53 * t52;
t60 = -t51 * t46 + t53 * t47;
t3 = -t30 * mrSges(5,1) - t60 * mrSges(5,2) + Ifges(5,5) * t44 + Ifges(5,6) * t43;
t61 = t44 * mrSges(5,1) + t43 * mrSges(5,2);
t88 = qJD(4) * t61;
t48 = t52 * pkin(3) + qJ(2);
t87 = m(5) * t48;
t84 = -t51 * t43 + t53 * t44;
t2 = (-Ifges(5,4) * t43 + (-Ifges(5,1) + Ifges(5,2)) * t44) * t43 + Ifges(5,4) * t44 ^ 2 + t48 * (-t43 * mrSges(5,1) + t44 * mrSges(5,2));
t70 = t54 * mrSges(4,2);
t66 = t2 * qJD(1);
t57 = -t52 * mrSges(4,1) + t61 - t70;
t15 = t87 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) - t57;
t64 = t15 * qJD(1);
t1 = (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t52) * t52 + t2 + (qJ(2) * mrSges(4,1) + (-Ifges(4,1) + Ifges(4,2)) * t52 + (-t61 + t87) * pkin(3) - Ifges(4,4) * t54) * t54;
t59 = t1 * qJD(1);
t45 = (mrSges(5,1) * t51 + mrSges(5,2) * t53) * pkin(3);
t56 = t45 * qJD(3);
t42 = t45 * qJD(4);
t4 = [t15 * qJD(2) + t1 * qJD(3) + t2 * qJD(4), t64, t3 * qJD(4) + t59 + (-mrSges(4,1) * t72 - Ifges(4,5) * t52 - Ifges(4,6) * t54 - t55 * t70 + (m(5) * (-t30 * t53 + t51 * t60) - t84 * mrSges(5,3)) * pkin(3) + t3) * qJD(3), t66 + (qJD(4) + qJD(3)) * t3; -t64, 0, (m(5) * t84 * pkin(3) + t57) * qJD(3) + t88, qJD(3) * t61 + t88; -t59, 0, -t42, -t42 - t56; -t66, 0, t56, 0;];
Cq = t4;
