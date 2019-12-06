% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:55
% EndTime: 2019-12-05 15:00:56
% DurationCPUTime: 0.44s
% Computational Cost: add. (1176->67), mult. (2890->116), div. (0->0), fcn. (3156->8), ass. (0->49)
t52 = sin(pkin(8));
t55 = sin(qJ(3));
t71 = cos(pkin(8));
t82 = cos(qJ(3));
t42 = t55 * t52 - t82 * t71;
t51 = sin(pkin(9));
t53 = cos(pkin(9));
t54 = sin(qJ(5));
t56 = cos(qJ(5));
t43 = t56 * t51 + t54 * t53;
t27 = t43 * t42;
t88 = -t54 * t51 + t56 * t53;
t29 = t88 * t42;
t90 = m(6) * (t27 * t88 - t43 * t29);
t73 = t51 ^ 2 + t53 ^ 2;
t89 = t73 * mrSges(5,3);
t86 = t43 ^ 2;
t85 = m(6) / 0.2e1;
t32 = t43 * mrSges(6,1) + mrSges(6,2) * t88;
t81 = t42 * t32;
t77 = pkin(6) + qJ(4);
t70 = t32 * qJD(3);
t69 = t32 * qJD(5);
t68 = qJD(1) * t90;
t67 = m(5) * t73;
t64 = t73 * t42;
t44 = t82 * t52 + t55 * t71;
t28 = t43 * t44;
t30 = t88 * t44;
t63 = (t43 * t28 + t30 * t88) * t85;
t45 = t77 * t51;
t46 = t77 * t53;
t33 = -t56 * t45 - t54 * t46;
t34 = -t54 * t45 + t56 * t46;
t13 = (t88 ^ 2 + t86) * mrSges(6,3) + t89 + m(6) * (-t33 * t43 + t34 * t88) + qJ(4) * t67;
t9 = t63 + 0.2e1 * (t67 / 0.4e1 - m(5) / 0.4e1 - m(6) / 0.4e1) * t44;
t61 = t9 * qJD(1) + t13 * qJD(3);
t60 = t27 * mrSges(6,1) / 0.2e1 + t29 * mrSges(6,2) / 0.2e1;
t35 = t42 * t44;
t4 = m(5) * (-t44 * t64 + t35) + m(6) * (-t28 * t27 - t30 * t29 + t35);
t58 = -t4 * qJD(1) - qJD(2) * t90 / 0.2e1;
t1 = -t81 / 0.2e1 + t60;
t48 = -t53 * pkin(4) - pkin(3);
t3 = -t86 * Ifges(6,4) + t48 * t32 + (Ifges(6,4) * t88 + (Ifges(6,1) - Ifges(6,2)) * t43) * t88;
t57 = -t1 * qJD(1) + t3 * qJD(3);
t8 = t63 + (m(5) + m(6) + t67) * t44 / 0.2e1;
t5 = qJD(3) * t90 / 0.2e1;
t2 = t81 / 0.2e1 + t60;
t6 = [t4 * qJD(3), t5, t8 * qJD(4) + t2 * qJD(5) - t58 + ((-t27 * t43 - t29 * t88) * mrSges(6,3) + (-t53 * mrSges(5,1) - mrSges(6,1) * t88 + t51 * mrSges(5,2) + t43 * mrSges(6,2) - mrSges(4,1)) * t44 + (mrSges(4,2) - t89) * t42 + m(5) * (-pkin(3) * t44 - qJ(4) * t64) + 0.2e1 * (t33 * t27 - t34 * t29 + t48 * t44) * t85) * qJD(3), t8 * qJD(3), t2 * qJD(3) + (-t30 * mrSges(6,1) + t28 * mrSges(6,2)) * qJD(5); t5, 0, t68 / 0.2e1, 0, -t69; t9 * qJD(4) - t1 * qJD(5) + t58, -t68 / 0.2e1, t13 * qJD(4) + t3 * qJD(5), t61, (-t34 * mrSges(6,1) - t33 * mrSges(6,2) + Ifges(6,5) * t88 - Ifges(6,6) * t43) * qJD(5) + t57; -t9 * qJD(3), 0, -t61 + t69, 0, t70; t1 * qJD(3), 0, -t32 * qJD(4) - t57, -t70, 0;];
Cq = t6;
