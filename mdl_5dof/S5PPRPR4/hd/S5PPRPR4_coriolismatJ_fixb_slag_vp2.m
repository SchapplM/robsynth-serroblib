% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:18
% EndTime: 2019-12-31 17:32:18
% DurationCPUTime: 0.34s
% Computational Cost: add. (866->66), mult. (2252->111), div. (0->0), fcn. (2208->6), ass. (0->48)
t51 = sin(pkin(8));
t52 = cos(pkin(8));
t53 = sin(qJ(5));
t79 = cos(qJ(5));
t41 = t79 * t51 + t53 * t52;
t55 = cos(qJ(3));
t31 = t41 * t55;
t58 = -t53 * t51 + t79 * t52;
t33 = t58 * t55;
t87 = m(6) * (t31 * t58 - t33 * t41);
t54 = sin(qJ(3));
t85 = t55 * t54;
t84 = t41 ^ 2;
t83 = m(6) / 0.2e1;
t82 = t54 / 0.2e1;
t77 = pkin(6) + qJ(4);
t49 = t51 ^ 2;
t50 = t52 ^ 2;
t76 = t49 + t50;
t24 = t41 * mrSges(6,1) + mrSges(6,2) * t58;
t74 = t24 * qJD(3);
t73 = t24 * qJD(5);
t72 = qJD(2) * t87;
t71 = m(5) * t82;
t69 = -t55 * t24 / 0.2e1;
t65 = t76 * mrSges(5,3);
t63 = t76 * qJ(4);
t46 = -t52 * pkin(4) - pkin(3);
t1 = -t84 * Ifges(6,4) + t46 * t24 - (-Ifges(6,4) * t58 + (-Ifges(6,1) + Ifges(6,2)) * t41) * t58;
t59 = t31 * mrSges(6,1) / 0.2e1 + t33 * mrSges(6,2) / 0.2e1;
t3 = t69 + t59;
t62 = t3 * qJD(2) + t1 * qJD(3);
t30 = t41 * t54;
t32 = t58 * t54;
t61 = (t30 * t41 + t32 * t58) * t83;
t11 = m(5) * (-0.1e1 + t76) * t85 + m(6) * (t30 * t31 + t32 * t33 - t85);
t57 = -t11 * qJD(2) - qJD(1) * t87 / 0.2e1;
t13 = t61 + (-m(6) / 0.2e1 + (t49 / 0.2e1 + t50 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t54;
t42 = t77 * t51;
t43 = t77 * t52;
t25 = -t79 * t42 - t53 * t43;
t26 = -t53 * t42 + t79 * t43;
t9 = (t58 ^ 2 + t84) * mrSges(6,3) + t65 + m(6) * (-t25 * t41 + t26 * t58) + m(5) * t63;
t56 = t13 * qJD(2) + t9 * qJD(3);
t12 = m(6) * t82 + t76 * t71 + t61 + t71;
t4 = qJD(3) * t87 / 0.2e1;
t2 = t69 - t59;
t5 = [0, t4, t72 / 0.2e1, 0, t73; t4, t11 * qJD(3), t12 * qJD(4) + t2 * qJD(5) - t57 + ((t31 * t41 + t33 * t58) * mrSges(6,3) + (-mrSges(4,2) + t65) * t55 + (-t52 * mrSges(5,1) - mrSges(6,1) * t58 + t51 * mrSges(5,2) + t41 * mrSges(6,2) - mrSges(4,1)) * t54 + 0.2e1 * (-t25 * t31 + t26 * t33 + t46 * t54) * t83 + m(5) * (-pkin(3) * t54 + t55 * t63)) * qJD(3), t12 * qJD(3), t2 * qJD(3) + (-t32 * mrSges(6,1) + t30 * mrSges(6,2)) * qJD(5); -t72 / 0.2e1, t13 * qJD(4) + t3 * qJD(5) + t57, t9 * qJD(4) + t1 * qJD(5), t56, (-t26 * mrSges(6,1) - t25 * mrSges(6,2) + Ifges(6,5) * t58 - Ifges(6,6) * t41) * qJD(5) + t62; 0, -t13 * qJD(3), -t56 + t73, 0, t74; 0, -t3 * qJD(3), -t24 * qJD(4) - t62, -t74, 0;];
Cq = t5;
