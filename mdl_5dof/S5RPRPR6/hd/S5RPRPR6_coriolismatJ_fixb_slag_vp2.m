% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:44
% EndTime: 2019-12-31 18:17:44
% DurationCPUTime: 0.30s
% Computational Cost: add. (602->72), mult. (1212->94), div. (0->0), fcn. (800->6), ass. (0->51)
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t65 = t42 * mrSges(6,1) + t44 * mrSges(6,2);
t59 = mrSges(5,3) + t65;
t81 = Ifges(6,1) - Ifges(6,2);
t43 = sin(qJ(3));
t54 = cos(pkin(8)) * pkin(1) + pkin(2);
t72 = cos(qJ(3));
t73 = pkin(1) * sin(pkin(8));
t25 = t43 * t73 - t72 * t54;
t26 = t43 * t54 + t72 * t73;
t40 = t44 ^ 2;
t64 = t42 ^ 2 + t40;
t55 = t64 * t26;
t80 = -(mrSges(4,1) - mrSges(5,2)) * t26 - mrSges(6,3) * t55 + (mrSges(4,2) - t59) * t25;
t74 = -t42 / 0.2e1;
t76 = m(5) / 0.4e1;
t75 = m(6) / 0.4e1;
t71 = Ifges(6,4) * t42;
t70 = Ifges(6,4) * t44;
t24 = qJ(4) + t26;
t69 = t24 * t25;
t67 = t42 * mrSges(6,2);
t66 = t44 * mrSges(6,1);
t63 = qJ(4) * t25;
t52 = -pkin(3) + t25;
t21 = -pkin(7) + t52;
t1 = -m(5) * (t52 * t26 - t69) - m(6) * (t21 * t55 - t69) - t80;
t62 = t1 * qJD(1);
t53 = t66 - t67;
t14 = t24 * t53;
t49 = t81 * t44 - t71;
t47 = t40 * Ifges(6,4) + t49 * t42;
t7 = -t14 + t47;
t61 = t7 * qJD(1);
t11 = 0.4e1 * (t76 + t75) * t24 + t59;
t60 = t11 * qJD(1);
t33 = qJ(4) * t53;
t56 = -t33 / 0.2e1 - t14 / 0.2e1;
t10 = -t33 + t47;
t2 = (t26 * mrSges(6,1) / 0.2e1 + t70) * t44 + (-t26 * mrSges(6,2) / 0.2e1 + t49) * t42 + t56;
t51 = t2 * qJD(1) + t10 * qJD(3);
t27 = (m(5) + m(6)) * qJ(4) + t59;
t46 = -t59 + (-m(5) / 0.2e1 - m(6) / 0.2e1) * (0.2e1 * qJ(4) + t26);
t48 = 0.2e1 * (t64 * t75 + t76) * t26;
t4 = t48 + t46;
t50 = t4 * qJD(1) - t27 * qJD(3);
t45 = -pkin(3) - pkin(7);
t5 = t48 - t46;
t3 = -0.2e1 * t71 * t74 + (-t67 / 0.2e1 + t66 / 0.2e1) * t26 - t56 + (0.2e1 * t81 * t74 - t70) * t44;
t6 = [-t1 * qJD(3) + t11 * qJD(4) - t7 * qJD(5), 0, t5 * qJD(4) + t3 * qJD(5) - t62 + (m(5) * (-pkin(3) * t26 - t63) + m(6) * (t45 * t55 - t63) + t80) * qJD(3), t5 * qJD(3) + t60, -t61 + t3 * qJD(3) + (-Ifges(6,5) * t42 - Ifges(6,6) * t44 - t65 * t21) * qJD(5); 0, 0, 0, 0, -t53 * qJD(5); -t4 * qJD(4) - t2 * qJD(5) + t62, 0, t27 * qJD(4) - t10 * qJD(5), -t50, ((-mrSges(6,2) * t45 - Ifges(6,6)) * t44 + (-mrSges(6,1) * t45 - Ifges(6,5)) * t42) * qJD(5) - t51; t4 * qJD(3) - t60, 0, t50, 0, -t65 * qJD(5); t2 * qJD(3) + t61, 0, t51, 0, 0;];
Cq = t6;
