% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:42
% EndTime: 2019-12-05 15:42:44
% DurationCPUTime: 0.82s
% Computational Cost: add. (3361->84), mult. (6728->122), div. (0->0), fcn. (7546->6), ass. (0->52)
t73 = sin(pkin(9));
t74 = cos(pkin(9));
t76 = sin(qJ(4));
t78 = cos(qJ(4));
t65 = -t73 * t76 + t74 * t78;
t66 = t73 * t78 + t74 * t76;
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t56 = t65 * t75 + t66 * t77;
t88 = t65 * t77 - t66 * t75;
t139 = t56 * mrSges(6,1) + t88 * mrSges(6,2);
t138 = t139 * qJD(5);
t99 = pkin(6) + qJ(3);
t69 = t99 * t73;
t70 = t99 * t74;
t58 = -t69 * t76 + t70 * t78;
t38 = pkin(7) * t65 + t58;
t57 = -t69 * t78 - t70 * t76;
t82 = -pkin(7) * t66 + t57;
t123 = -t38 * t75 + t77 * t82;
t29 = t38 * t77 + t75 * t82;
t4 = -t29 * mrSges(6,1) - mrSges(6,2) * t123 + Ifges(6,5) * t88 - Ifges(6,6) * t56;
t136 = t4 * qJD(5);
t115 = Ifges(6,4) * t56;
t117 = t56 / 0.2e1;
t118 = -t56 / 0.2e1;
t119 = t88 / 0.2e1;
t90 = -pkin(3) * t74 - pkin(2);
t59 = -pkin(4) * t65 + t90;
t3 = t59 * t139 + (0.2e1 * Ifges(6,4) * t88 + (Ifges(6,1) - Ifges(6,2)) * t56) * t119 + (Ifges(6,2) * t88 + t115) * t118 + (Ifges(6,1) * t88 - t115) * t117;
t122 = t66 ^ 2;
t120 = m(6) * pkin(4);
t116 = t66 * pkin(4);
t106 = t56 * t77;
t103 = t88 * t75;
t94 = t3 * qJD(2);
t62 = t66 * mrSges(5,1);
t81 = -t65 * mrSges(5,2) - t139 - t62;
t14 = (t66 / 0.2e1 + t106 / 0.2e1 - t103 / 0.2e1) * t120 - t81;
t92 = t14 * qJD(2);
t91 = t139 * qJD(2);
t1 = -t122 * Ifges(5,4) + t90 * t62 + (Ifges(5,4) * t65 + t90 * mrSges(5,2) + (Ifges(5,1) - Ifges(5,2)) * t66) * t65 + (m(6) * t59 - mrSges(6,1) * t88 + mrSges(6,2) * t56) * t116 + t3;
t85 = t1 * qJD(2);
t7 = (t56 ^ 2 + t88 ^ 2) * mrSges(6,3) + (t65 ^ 2 + t122) * mrSges(5,3) + m(6) * (-t123 * t56 + t29 * t88) + m(5) * (-t57 * t66 + t58 * t65) + (m(4) * qJ(3) + mrSges(4,3)) * (t73 ^ 2 + t74 ^ 2);
t83 = t7 * qJD(2);
t80 = (t103 - t106) * t120;
t5 = (t118 + t117) * Ifges(6,6) + (t119 - t88 / 0.2e1) * Ifges(6,5);
t68 = (mrSges(6,1) * t75 + mrSges(6,2) * t77) * pkin(4);
t79 = -qJD(2) * t5 + qJD(4) * t68;
t67 = t68 * qJD(5);
t21 = m(6) * t116 / 0.2e1 + t80 / 0.2e1;
t2 = [0, 0, 0, (t80 + t81) * qJD(4) - t138, -qJD(4) * t139 - t138; 0, qJD(3) * t7 + qJD(4) * t1 + qJD(5) * t3, qJD(4) * t21 + t83, t21 * qJD(3) + t136 + t85 + (-t58 * mrSges(5,1) - t57 * mrSges(5,2) + Ifges(5,5) * t65 - Ifges(5,6) * t66 + (m(6) * (t123 * t75 - t29 * t77) + (-t56 * t75 - t77 * t88) * mrSges(6,3)) * pkin(4) + t4) * qJD(4), qJD(4) * t4 + t136 + t94; 0, qJD(4) * t14 + t138 - t83, 0, t92, t91; 0, -qJD(3) * t14 + qJD(5) * t5 - t85, -t92, -t67, -t67 - t79; 0, -qJD(3) * t139 - qJD(4) * t5 - t94, -t91, t79, 0;];
Cq = t2;
