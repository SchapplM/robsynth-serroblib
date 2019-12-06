% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:06
% EndTime: 2019-12-05 16:19:10
% DurationCPUTime: 1.21s
% Computational Cost: add. (3885->107), mult. (7717->152), div. (0->0), fcn. (8501->6), ass. (0->59)
t113 = cos(pkin(9));
t89 = sin(pkin(9));
t91 = sin(qJ(3));
t93 = cos(qJ(3));
t74 = t113 * t93 - t89 * t91;
t75 = -t113 * t91 - t89 * t93;
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t103 = t92 * t74 + t75 * t90;
t120 = -qJ(4) - pkin(6);
t81 = t120 * t91;
t82 = t120 * t93;
t144 = t113 * t81 + t89 * t82;
t146 = pkin(7) * t75 + t144;
t62 = t113 * t82 - t89 * t81;
t43 = pkin(7) * t74 - t62;
t151 = t146 * t92 - t43 * t90;
t33 = t146 * t90 + t43 * t92;
t61 = t74 * t90 - t75 * t92;
t4 = -t33 * mrSges(6,1) - t151 * mrSges(6,2) + Ifges(6,5) * t103 - Ifges(6,6) * t61;
t163 = t4 * qJD(5);
t160 = t61 * mrSges(6,1) + t103 * mrSges(6,2);
t161 = qJD(5) * t160;
t158 = m(5) * pkin(3);
t155 = -t75 * mrSges(5,1) + t74 * mrSges(5,2);
t132 = Ifges(6,4) * t61;
t136 = t61 / 0.2e1;
t137 = -t61 / 0.2e1;
t138 = t103 / 0.2e1;
t87 = -pkin(3) * t93 - pkin(2);
t67 = -t74 * pkin(4) + t87;
t3 = t67 * t160 + (0.2e1 * Ifges(6,4) * t103 + (Ifges(6,1) - Ifges(6,2)) * t61) * t138 + (Ifges(6,2) * t103 + t132) * t137 + (Ifges(6,1) * t103 - t132) * t136;
t133 = t91 * pkin(3);
t108 = m(5) * t133;
t139 = m(6) / 0.2e1;
t135 = pkin(3) * t89;
t114 = t3 * qJD(2);
t107 = (t113 * t75 + t74 * t89) * t158 / 0.2e1;
t105 = t113 * pkin(3);
t86 = t105 + pkin(4);
t69 = -t90 * t135 + t86 * t92;
t70 = t92 * t135 + t86 * t90;
t25 = t103 * t70 - t61 * t69;
t68 = -pkin(4) * t75 + t133;
t95 = -t160 - t155;
t13 = -t108 / 0.2e1 + t107 + 0.2e1 * (t25 / 0.4e1 - t68 / 0.4e1) * m(6) + t95;
t112 = t13 * qJD(2);
t110 = t160 * qJD(2);
t23 = -mrSges(6,1) * t70 - mrSges(6,2) * t69;
t109 = t23 * qJD(5);
t1 = (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t91) * t91 + (-mrSges(5,2) * t133 - Ifges(5,4) * t75) * t75 + (-mrSges(4,2) * pkin(2) + Ifges(4,4) * t93 + (Ifges(4,1) - Ifges(4,2)) * t91) * t93 + (-mrSges(5,1) * t133 + Ifges(5,4) * t74 + (-Ifges(5,1) + Ifges(5,2)) * t75) * t74 + (t108 + t155) * t87 + t3 + (m(6) * t67 - mrSges(6,1) * t103 + t61 * mrSges(6,2)) * t68;
t99 = t1 * qJD(2);
t7 = (t103 ^ 2 + t61 ^ 2) * mrSges(6,3) + (t74 ^ 2 + t75 ^ 2) * mrSges(5,3) + m(6) * (t103 * t33 - t151 * t61) + m(5) * (t144 * t75 - t62 * t74);
t97 = qJD(2) * t7;
t96 = t25 * t139 + t107;
t5 = (t137 + t136) * Ifges(6,6) + (t138 - t103 / 0.2e1) * Ifges(6,5);
t94 = t5 * qJD(2) + t23 * qJD(3);
t17 = t68 * t139 + t108 / 0.2e1 + t96;
t2 = [0, 0, -t161 + (-t91 * mrSges(4,1) - t93 * mrSges(4,2) + t95 + 0.2e1 * t96) * qJD(3), 0, -qJD(3) * t160 - t161; 0, qJD(3) * t1 + qJD(4) * t7 + qJD(5) * t3, (Ifges(5,5) * t74 + Ifges(5,6) * t75 + Ifges(4,5) * t93 - Ifges(4,6) * t91 + m(6) * (t151 * t70 - t33 * t69) + t62 * mrSges(5,1) - t144 * mrSges(5,2) + (t113 * t62 + t144 * t89) * t158 + (-t93 * mrSges(4,1) + t91 * mrSges(4,2)) * pkin(6) + (-t103 * t69 - t61 * t70) * mrSges(6,3) + (-t74 * t105 + t75 * t135) * mrSges(5,3) + t4) * qJD(3) + t17 * qJD(4) + t163 + t99, qJD(3) * t17 + t97, t4 * qJD(3) + t114 + t163; 0, qJD(4) * t13 + qJD(5) * t5 - t99, t109, t112, t94 + t109; 0, -qJD(3) * t13 + t161 - t97, -t112, 0, t110; 0, -t5 * qJD(3) - qJD(4) * t160 - t114, -t94, -t110, 0;];
Cq = t2;
