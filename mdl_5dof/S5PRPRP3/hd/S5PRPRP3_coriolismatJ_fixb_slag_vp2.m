% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:42
% EndTime: 2019-12-05 15:32:43
% DurationCPUTime: 0.35s
% Computational Cost: add. (528->72), mult. (1212->98), div. (0->0), fcn. (1036->6), ass. (0->44)
t45 = sin(qJ(4));
t43 = t45 ^ 2;
t46 = cos(qJ(4));
t44 = t46 ^ 2;
t66 = t44 + t43;
t62 = sin(pkin(8));
t82 = t62 * pkin(2);
t42 = pkin(6) + t82;
t64 = qJ(5) + t42;
t25 = t64 * t45;
t26 = t64 * t46;
t83 = m(6) * (t25 * t45 + t26 * t46);
t63 = cos(pkin(8));
t81 = t63 * pkin(2);
t78 = m(6) * pkin(4);
t55 = mrSges(6,1) + t78;
t73 = sin(qJ(2));
t74 = cos(qJ(2));
t38 = -t62 * t74 - t63 * t73;
t75 = m(6) * t38;
t69 = t46 * mrSges(6,2);
t68 = -mrSges(5,2) - mrSges(6,2);
t67 = Ifges(5,4) + Ifges(6,4);
t37 = t62 * t73 - t63 * t74;
t2 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (-0.1e1 + t66) * t38 * t37;
t65 = t2 * qJD(1);
t61 = qJD(4) * t45;
t60 = qJD(4) * t46;
t27 = t55 * t45 + t69;
t59 = t27 * qJD(2);
t54 = mrSges(5,1) + t55;
t53 = pkin(3) + t81;
t52 = t45 * mrSges(5,1) + t46 * mrSges(5,2);
t51 = t46 * mrSges(6,1) - t45 * mrSges(6,2);
t50 = t45 * mrSges(6,1) + t69;
t10 = t66 * mrSges(6,3) + t83;
t8 = 0.2e1 * (-0.1e1 / 0.4e1 + t43 / 0.4e1 + t44 / 0.4e1) * t75;
t49 = -t8 * qJD(1) + t10 * qJD(2);
t48 = t46 * pkin(4) + t53;
t1 = -t67 * t44 + t48 * t50 + t53 * t52 + (pkin(4) * t51 + t48 * t78 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t46 + t67 * t45) * t45;
t47 = t1 * qJD(2);
t9 = (-t66 / 0.2e1 - 0.1e1 / 0.2e1) * t75;
t4 = (t52 / 0.2e1 + t50 / 0.2e1 + (mrSges(5,2) / 0.2e1 + mrSges(6,2) / 0.2e1) * t46 + (t78 + mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1) * t45) * t37;
t3 = [t2 * qJD(2), t4 * qJD(4) + t9 * qJD(5) + t65 + (-t73 * mrSges(3,1) - t74 * mrSges(3,2) + (m(4) * t81 + m(5) * t53 + m(6) * t48 + t46 * mrSges(5,1) - t45 * mrSges(5,2) + mrSges(4,1) + t51) * t38 + (-m(4) * t82 + mrSges(4,2) - t83 + (-m(5) * t42 - mrSges(5,3) - mrSges(6,3)) * t66) * t37) * qJD(2), 0, t4 * qJD(2) + (t45 * t68 + t46 * t54) * qJD(4) * t38, t9 * qJD(2); -t8 * qJD(5) - t65, -t1 * qJD(4) + t10 * qJD(5), 0, (t25 * mrSges(6,2) - t55 * t26) * qJD(4) + (mrSges(5,2) * t42 - Ifges(5,6) - Ifges(6,6)) * t61 + (-mrSges(5,1) * t42 - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t60 - t47, t49; 0, 0, 0, -t54 * t61 + t60 * t68, 0; 0, -t27 * qJD(5) + t47, 0, 0, -t59; t8 * qJD(2), t27 * qJD(4) - t49, 0, t59, 0;];
Cq = t3;
