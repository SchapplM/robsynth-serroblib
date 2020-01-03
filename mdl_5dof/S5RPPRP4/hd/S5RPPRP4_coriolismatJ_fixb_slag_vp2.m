% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:04
% EndTime: 2019-12-31 17:52:05
% DurationCPUTime: 0.28s
% Computational Cost: add. (594->73), mult. (1049->89), div. (0->0), fcn. (717->4), ass. (0->41)
t40 = sin(qJ(4));
t35 = t40 ^ 2;
t41 = cos(qJ(4));
t36 = t41 ^ 2;
t56 = t36 + t35;
t72 = -t56 + 0.1e1;
t37 = sin(pkin(7));
t42 = -pkin(1) - pkin(2);
t38 = cos(pkin(7));
t66 = t38 * qJ(2);
t24 = t37 * t42 - pkin(6) + t66;
t65 = qJ(5) - t24;
t15 = t65 * t40;
t16 = t65 * t41;
t71 = m(6) * (t15 * t40 + t16 * t41);
t58 = mrSges(5,2) + mrSges(6,2);
t70 = t58 * t40;
t69 = -mrSges(5,1) - mrSges(6,1);
t48 = m(6) * pkin(4) + mrSges(6,1);
t46 = -mrSges(5,1) - t48;
t68 = t46 * t40 - t58 * t41;
t62 = m(6) * t37;
t57 = Ifges(6,4) + Ifges(5,4);
t31 = t37 * qJ(2);
t45 = t38 * t42 - pkin(3) - t31;
t18 = t41 * pkin(4) - t45;
t1 = (t45 * mrSges(5,2) - t18 * mrSges(6,2) + t57 * t41) * t41 + (t45 * mrSges(5,1) + (pkin(4) * mrSges(6,2) - t57) * t40 - t48 * t18 + (-pkin(4) * mrSges(6,1) + Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,2)) * t41) * t40;
t55 = t1 * qJD(1);
t2 = -m(3) * qJ(2) - mrSges(3,3) + (-m(4) * t31 + m(5) * t45 - m(6) * t18 + t69 * t41 - mrSges(4,1) + t70) * t37 + (-m(4) * t66 - mrSges(4,2) + t71 + (-m(5) * t24 + mrSges(5,3) + mrSges(6,3)) * t56) * t38;
t54 = t2 * qJD(1);
t3 = t68 * t38;
t53 = t3 * qJD(1);
t7 = t56 * mrSges(6,3) + t71;
t52 = t7 * qJD(1);
t14 = (-t35 / 0.2e1 - t36 / 0.2e1 - 0.1e1 / 0.2e1) * t62;
t51 = t14 * qJD(1);
t19 = t41 * mrSges(6,2) + t48 * t40;
t50 = t19 * qJD(1);
t13 = t72 * t62 / 0.2e1;
t4 = (-mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1 - t69 / 0.2e1) * t40 * t38;
t5 = [-t2 * qJD(2) + t1 * qJD(4) + t7 * qJD(5), -t54 - 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t72 * t38 * t37 * qJD(2) + t4 * qJD(4) + t13 * qJD(5), 0, t4 * qJD(2) + t55 + (-t15 * mrSges(6,2) + (mrSges(5,2) * t24 + Ifges(5,6) + Ifges(6,6)) * t40 + (-mrSges(5,1) * t24 + mrSges(6,3) * pkin(4) - Ifges(5,5) - Ifges(6,5)) * t41 + t48 * t16) * qJD(4), t13 * qJD(2) + t52; -t3 * qJD(4) + t14 * qJD(5) + t54, 0, 0, -t53 + (t46 * t41 + t70) * qJD(4) * t37, t51; 0, 0, 0, t68 * qJD(4), 0; t3 * qJD(2) + t19 * qJD(5) - t55, t53, 0, 0, t50; -t14 * qJD(2) - t19 * qJD(4) - t52, -t51, 0, -t50, 0;];
Cq = t5;
