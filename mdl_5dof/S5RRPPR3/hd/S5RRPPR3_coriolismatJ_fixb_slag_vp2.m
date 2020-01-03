% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:30
% EndTime: 2019-12-31 19:26:30
% DurationCPUTime: 0.34s
% Computational Cost: add. (564->87), mult. (1308->115), div. (0->0), fcn. (852->6), ass. (0->60)
t47 = sin(qJ(5));
t94 = -t47 / 0.2e1;
t49 = cos(qJ(5));
t93 = Ifges(6,2) * t47 / 0.2e1 - Ifges(6,4) * t49 + Ifges(6,1) * t94;
t45 = sin(pkin(8));
t50 = cos(qJ(2));
t46 = cos(pkin(8));
t48 = sin(qJ(2));
t77 = t46 * t48;
t26 = (t45 * t50 + t77) * pkin(1);
t37 = t45 * t48 * pkin(1);
t85 = t50 * pkin(1);
t27 = t46 * t85 - t37;
t71 = t47 ^ 2 + t49 ^ 2;
t62 = t71 * t26;
t72 = t47 * mrSges(6,1) + t49 * mrSges(6,2);
t92 = (t48 * mrSges(3,1) + t50 * mrSges(3,2)) * pkin(1) - t27 * t72 + (mrSges(4,1) - mrSges(5,2)) * t26 + mrSges(6,3) * t62;
t87 = m(6) / 0.4e1;
t88 = m(5) / 0.4e1;
t91 = 0.4e1 * t88 + 0.4e1 * t87;
t84 = Ifges(6,4) * t47;
t40 = pkin(2) + t85;
t67 = pkin(1) * t77;
t54 = t45 * t40 + t67;
t22 = qJ(4) + t54;
t82 = t22 * t27;
t74 = t49 * mrSges(6,1);
t76 = t47 * mrSges(6,2);
t32 = t74 - t76;
t81 = t22 * t32;
t39 = t45 * pkin(2) + qJ(4);
t79 = t39 * t27;
t78 = t39 * t32;
t61 = t46 * t40 - t37;
t59 = -pkin(3) - t61;
t21 = -pkin(7) + t59;
t1 = (mrSges(4,2) - mrSges(5,3)) * t27 - m(6) * (t21 * t62 + t82) - m(4) * (-t61 * t26 + t54 * t27) - m(5) * (t59 * t26 + t82) + t92;
t70 = t1 * qJD(1);
t33 = -Ifges(6,2) * t49 - t84;
t36 = Ifges(6,1) * t49 - t84;
t60 = (t33 + t36) * t94 + t93 * t49;
t7 = t60 + t81;
t69 = t7 * qJD(1);
t64 = mrSges(5,3) + t72;
t11 = t22 * t91 + t64;
t68 = t11 * qJD(1);
t63 = -t46 * pkin(2) - pkin(3);
t57 = -Ifges(6,5) * t47 - Ifges(6,6) * t49;
t10 = t60 + t78;
t2 = (-t22 / 0.2e1 - t39 / 0.2e1) * t32 + (t26 * mrSges(6,1) / 0.2e1 - t93) * t49 + (-t26 * mrSges(6,2) / 0.2e1 + t33 / 0.2e1 + t36 / 0.2e1) * t47;
t56 = t2 * qJD(1) - t10 * qJD(2);
t13 = t39 * t91 + t64;
t51 = -t64 + (-m(6) / 0.2e1 - m(5) / 0.2e1) * (t67 + 0.2e1 * qJ(4) + (pkin(2) + t40) * t45);
t53 = 0.2e1 * (t71 * t87 + t88) * t26;
t4 = t53 + t51;
t55 = t4 * qJD(1) - t13 * qJD(2);
t38 = -pkin(7) + t63;
t5 = t53 - t51;
t3 = t81 / 0.2e1 + t78 / 0.2e1 + (-t76 / 0.2e1 + t74 / 0.2e1) * t26 + t60;
t6 = [-t1 * qJD(2) + t11 * qJD(4) + t7 * qJD(5), t5 * qJD(4) + t3 * qJD(5) - t70 + (-t27 * mrSges(4,2) + t27 * mrSges(5,3) + m(6) * (t38 * t62 + t79) + m(4) * (-t46 * t26 + t45 * t27) * pkin(2) + m(5) * (t63 * t26 + t79) - t92) * qJD(2), 0, t5 * qJD(2) + t68, t69 + t3 * qJD(2) + (-t21 * t72 + t57) * qJD(5); -t4 * qJD(4) - t2 * qJD(5) + t70, t13 * qJD(4) + t10 * qJD(5), 0, -t55, (-t38 * t72 + t57) * qJD(5) - t56; 0, 0, 0, 0, -t32 * qJD(5); t4 * qJD(2) - t68, t55, 0, 0, -t72 * qJD(5); t2 * qJD(2) - t69, t56, 0, 0, 0;];
Cq = t6;
