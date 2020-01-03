% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:22
% EndTime: 2019-12-31 17:46:23
% DurationCPUTime: 0.47s
% Computational Cost: add. (1475->78), mult. (2725->116), div. (0->0), fcn. (2658->6), ass. (0->52)
t101 = -m(6) / 0.2e1;
t63 = sin(pkin(8));
t65 = cos(pkin(8));
t68 = sin(qJ(5));
t69 = cos(qJ(5));
t103 = -t68 * t63 + t69 * t65;
t53 = t69 * t63 + t68 * t65;
t66 = cos(pkin(7));
t44 = t53 * t66;
t46 = t103 * t66;
t18 = -t103 * t44 + t53 * t46;
t107 = m(6) * t18 * qJD(1);
t64 = sin(pkin(7));
t97 = -pkin(1) - pkin(2);
t54 = t66 * qJ(2) + t64 * t97 - qJ(4);
t61 = t63 ^ 2;
t62 = t65 ^ 2;
t91 = t61 + t62;
t82 = m(5) * t91;
t106 = t54 * t82;
t100 = m(6) / 0.2e1;
t104 = t91 * mrSges(5,3);
t102 = t103 ^ 2;
t96 = pkin(6) - t54;
t37 = -t53 * mrSges(6,1) - mrSges(6,2) * t103;
t71 = t64 * qJ(2) - t66 * t97 + pkin(3);
t70 = t65 * pkin(4) + t71;
t1 = -Ifges(6,4) * t102 - t70 * t37 + (-(Ifges(6,1) - Ifges(6,2)) * t103 + Ifges(6,4) * t53) * t53;
t87 = t1 * qJD(1);
t43 = t53 * t64;
t45 = t103 * t64;
t77 = (-t103 * t45 - t43 * t53) * t100;
t15 = t77 + (t101 + (-t61 / 0.2e1 - t62 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t64;
t86 = t15 * qJD(1);
t85 = t37 * qJD(1);
t84 = t37 * qJD(5);
t83 = -t66 * t37 / 0.2e1;
t75 = t44 * mrSges(6,1) / 0.2e1 + t46 * mrSges(6,2) / 0.2e1;
t3 = t83 + t75;
t74 = t3 * qJD(1);
t41 = t96 * t63;
t42 = t96 * t65;
t23 = t69 * t41 + t68 * t42;
t24 = t68 * t41 - t69 * t42;
t5 = -m(6) * (-t23 * t44 + t24 * t46) - mrSges(3,3) + t66 * t104 + (-mrSges(4,2) - t106) * t66 + (-m(5) * t71 - m(6) * t70 - t65 * mrSges(5,1) - mrSges(6,1) * t103 + t63 * mrSges(5,2) + t53 * mrSges(6,2) - mrSges(4,1)) * t64 + (-m(4) * (t64 ^ 2 + t66 ^ 2) - m(3)) * qJ(2) + (t103 * t46 + t44 * t53) * mrSges(6,3);
t73 = t18 * qJD(3) * t101 + t5 * qJD(1);
t9 = (t53 ^ 2 + t102) * mrSges(6,3) + t104 + m(6) * (-t103 * t24 + t23 * t53) - t106;
t72 = t9 * qJD(1);
t14 = t77 + (m(5) + m(6) - t82) * t64 / 0.2e1;
t4 = t18 * qJD(2) * t100;
t2 = t83 - t75;
t6 = [-t5 * qJD(2) + t9 * qJD(4) - t1 * qJD(5), 0.2e1 * ((t43 * t44 + t45 * t46) * t100 + (t101 + m(5) * (-0.1e1 + t91) / 0.2e1) * t66 * t64) * qJD(2) + t14 * qJD(4) + t2 * qJD(5) - t73, t4, t14 * qJD(2) + t72, -t87 + t2 * qJD(2) + (-t24 * mrSges(6,1) - t23 * mrSges(6,2) - Ifges(6,5) * t103 + Ifges(6,6) * t53) * qJD(5); t15 * qJD(4) + t3 * qJD(5) + t73, 0, -t107 / 0.2e1, t86, (-t45 * mrSges(6,1) + t43 * mrSges(6,2)) * qJD(5) + t74; t4, t107 / 0.2e1, 0, 0, t84; -t15 * qJD(2) - t72 + t84, -t86, 0, 0, t85; -t3 * qJD(2) - t37 * qJD(4) + t87, -t74, 0, -t85, 0;];
Cq = t6;
