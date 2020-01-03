% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:29
% EndTime: 2019-12-31 17:03:29
% DurationCPUTime: 0.30s
% Computational Cost: add. (327->74), mult. (747->105), div. (0->0), fcn. (384->4), ass. (0->50)
t39 = sin(qJ(4));
t78 = -t39 / 0.2e1;
t37 = t39 ^ 2;
t41 = cos(qJ(4));
t38 = t41 ^ 2;
t57 = t37 + t38;
t77 = t57 * mrSges(5,3) + mrSges(3,1);
t58 = t39 * mrSges(5,1) + t41 * mrSges(5,2);
t50 = mrSges(4,3) + t58;
t76 = Ifges(5,2) * t39 / 0.2e1 - Ifges(5,4) * t41 + Ifges(5,1) * t78;
t74 = m(5) / 0.4e1 + m(4) / 0.4e1;
t13 = -(-m(5) - m(4)) * qJ(3) + t50;
t73 = m(4) / 0.2e1;
t71 = m(5) / 0.2e1;
t40 = sin(qJ(2));
t67 = t40 * pkin(1);
t42 = cos(qJ(2));
t66 = t42 * pkin(1);
t65 = Ifges(5,4) * t39;
t60 = t41 * mrSges(5,1);
t61 = t39 * mrSges(5,2);
t20 = t60 - t61;
t31 = qJ(3) + t67;
t63 = t31 * t20;
t62 = t31 * t42;
t59 = t42 * mrSges(3,2);
t56 = qJ(3) * t20;
t49 = -pkin(2) - t66;
t29 = -pkin(6) + t49;
t48 = t57 * t40;
t51 = mrSges(4,2) * t67 + t50 * t66;
t1 = -t51 + (-m(5) * (t29 * t48 + t62) - m(4) * (t49 * t40 + t62) + t59) * pkin(1) + t77 * t67;
t55 = t1 * qJD(1);
t21 = -Ifges(5,2) * t41 - t65;
t24 = Ifges(5,1) * t41 - t65;
t46 = (t21 + t24) * t78 + t76 * t41;
t4 = t46 + t63;
t54 = t4 * qJD(1);
t9 = 0.4e1 * t74 * t31 + t50;
t53 = t9 * qJD(1);
t2 = (-qJ(3) / 0.2e1 - t31 / 0.2e1) * t20 + (mrSges(5,1) * t67 / 0.2e1 - t76) * t41 + (-mrSges(5,2) * t67 / 0.2e1 + t24 / 0.2e1 + t21 / 0.2e1) * t39;
t5 = t46 + t56;
t45 = t2 * qJD(1) - t5 * qJD(2);
t6 = (t37 / 0.2e1 + t38 / 0.2e1 - 0.1e1 / 0.2e1) * m(5) * t67 - t13;
t44 = t6 * qJD(1) - t13 * qJD(2);
t43 = -pkin(2) - pkin(6);
t32 = qJ(3) * t66;
t7 = (t57 * t71 + t73) * t67 + t50 + t74 * (0.4e1 * qJ(3) + 0.2e1 * t67);
t3 = t56 / 0.2e1 + t63 / 0.2e1 + (-t61 / 0.2e1 + t60 / 0.2e1) * t67 + t46;
t8 = [-t1 * qJD(2) + t9 * qJD(3) + t4 * qJD(4), t7 * qJD(3) + t3 * qJD(4) - t55 + (t51 + (-t77 * t40 - t59) * pkin(1) + 0.2e1 * (t43 * pkin(1) * t48 + t32) * t71 + 0.2e1 * (-pkin(2) * t67 + t32) * t73) * qJD(2), t7 * qJD(2) + t53, t54 + t3 * qJD(2) + (-Ifges(5,5) * t39 - Ifges(5,6) * t41 - t58 * t29) * qJD(4); -t6 * qJD(3) - t2 * qJD(4) + t55, t13 * qJD(3) + t5 * qJD(4), -t44, ((-mrSges(5,2) * t43 - Ifges(5,6)) * t41 + (-mrSges(5,1) * t43 - Ifges(5,5)) * t39) * qJD(4) - t45; t6 * qJD(2) - t53, t44, 0, -t58 * qJD(4); t2 * qJD(2) - t54, t45, 0, 0;];
Cq = t8;
