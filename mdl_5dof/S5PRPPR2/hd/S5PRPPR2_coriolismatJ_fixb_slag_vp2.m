% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:02
% EndTime: 2019-12-05 15:24:04
% DurationCPUTime: 0.56s
% Computational Cost: add. (1347->72), mult. (3069->120), div. (0->0), fcn. (3333->8), ass. (0->52)
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t79 = sin(qJ(2));
t81 = cos(qJ(2));
t41 = t72 * t79 - t73 * t81;
t54 = sin(pkin(9));
t55 = cos(pkin(9));
t78 = sin(qJ(5));
t80 = cos(qJ(5));
t44 = -t80 * t54 - t78 * t55;
t25 = t44 * t41;
t86 = t78 * t54 - t80 * t55;
t27 = t86 * t41;
t89 = m(6) * (t25 * t86 - t44 * t27);
t75 = t54 ^ 2 + t55 ^ 2;
t66 = m(5) * t75;
t88 = t72 * pkin(2);
t87 = t73 * pkin(2);
t85 = t44 ^ 2;
t42 = -t72 * t81 - t73 * t79;
t84 = 0.2e1 * t42;
t46 = qJ(4) + t88;
t83 = pkin(6) + t46;
t33 = t41 * t42;
t32 = -t44 * mrSges(6,1) - mrSges(6,2) * t86;
t71 = t32 * qJD(2);
t70 = t32 * qJD(5);
t69 = -m(6) / 0.4e1 - m(5) / 0.4e1;
t68 = qJD(1) * t89;
t67 = t41 * t32 / 0.2e1;
t61 = -pkin(3) - t87;
t37 = t83 * t54;
t38 = t83 * t55;
t30 = -t80 * t37 - t78 * t38;
t31 = -t78 * t37 + t80 * t38;
t10 = (t86 ^ 2 + t85) * mrSges(6,3) + m(6) * (t30 * t44 - t31 * t86) + (m(5) * t46 + mrSges(5,3)) * t75;
t26 = t44 * t42;
t28 = t86 * t42;
t58 = m(6) * (-t44 * t26 - t28 * t86);
t11 = -t58 / 0.2e1 + (t66 / 0.4e1 + t69) * t84;
t60 = -t11 * qJD(1) + t10 * qJD(2);
t59 = t25 * mrSges(6,1) / 0.2e1 + t27 * mrSges(6,2) / 0.2e1;
t4 = m(6) * (t26 * t25 + t28 * t27 - t33) + m(5) * (t75 - 0.1e1) * t33;
t57 = -t4 * qJD(1) - qJD(3) * t89 / 0.2e1;
t2 = t67 + t59;
t45 = -t55 * pkin(4) + t61;
t3 = -Ifges(6,4) * t85 + t45 * t32 - (-Ifges(6,4) * t86 + (-Ifges(6,1) + Ifges(6,2)) * t44) * t86;
t56 = -t2 * qJD(1) - t3 * qJD(2);
t12 = t58 / 0.2e1 - t42 * t66 / 0.2e1 + t69 * t84;
t5 = qJD(2) * t89 / 0.2e1;
t1 = t67 - t59;
t6 = [t4 * qJD(2), (m(6) * (-t30 * t25 + t31 * t27) - t81 * mrSges(3,2) - t79 * mrSges(3,1) + (m(4) * t87 - m(5) * t61 - m(6) * t45 + t55 * mrSges(5,1) - mrSges(6,1) * t86 - t54 * mrSges(5,2) + t44 * mrSges(6,2) + mrSges(4,1)) * t42 + (-t25 * t44 - t27 * t86) * mrSges(6,3) + (-m(4) * t88 - t75 * mrSges(5,3) - t46 * t66 + mrSges(4,2)) * t41) * qJD(2) + t12 * qJD(4) + t1 * qJD(5) - t57, t5, t12 * qJD(2), t1 * qJD(2) + (-t28 * mrSges(6,1) + t26 * mrSges(6,2)) * qJD(5); -t11 * qJD(4) + t2 * qJD(5) + t57, t10 * qJD(4) + t3 * qJD(5), -t68 / 0.2e1, t60, (-t31 * mrSges(6,1) - t30 * mrSges(6,2) - Ifges(6,5) * t86 + Ifges(6,6) * t44) * qJD(5) - t56; t5, t68 / 0.2e1, 0, 0, -t70; t11 * qJD(2), -t60 + t70, 0, 0, t71; -t2 * qJD(2), -t32 * qJD(4) + t56, 0, -t71, 0;];
Cq = t6;
