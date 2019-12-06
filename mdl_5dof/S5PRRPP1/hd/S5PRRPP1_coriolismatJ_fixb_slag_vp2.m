% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:14
% EndTime: 2019-12-05 16:06:17
% DurationCPUTime: 0.71s
% Computational Cost: add. (1233->87), mult. (2659->107), div. (0->0), fcn. (2518->4), ass. (0->50)
t92 = sin(qJ(3));
t80 = t92 * pkin(3);
t104 = m(5) * t80;
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t93 = cos(qJ(3));
t41 = t84 * t92 - t85 * t93;
t43 = -t84 * t93 - t85 * t92;
t54 = -t93 * pkin(3) - pkin(2);
t21 = t41 * pkin(4) + t43 * qJ(5) + t54;
t103 = m(6) * t21 + mrSges(6,1) * t41 + mrSges(6,3) * t43;
t102 = -t43 * mrSges(5,1) - t41 * mrSges(5,2);
t77 = t84 * pkin(3);
t48 = t77 + qJ(5);
t47 = m(6) * t48 + mrSges(6,3);
t99 = -Ifges(6,5) + Ifges(5,4);
t23 = m(6) * t43;
t98 = qJD(5) * t23;
t79 = t92 * pkin(6);
t62 = -t92 * qJ(4) - t79;
t57 = t93 * pkin(6);
t87 = t93 * qJ(4) + t57;
t30 = -t85 * t62 + t84 * t87;
t97 = t84 * t62 + t85 * t87;
t95 = m(6) / 0.2e1;
t94 = m(5) * pkin(3);
t91 = t41 * mrSges(6,2);
t60 = (-t84 * t41 + t85 * t43) * t94;
t78 = t85 * pkin(3);
t53 = -t78 - pkin(4);
t66 = m(6) * (-t48 * t41 - t53 * t43);
t58 = t66 / 0.2e1 + t60 / 0.2e1;
t24 = -t43 * pkin(4) + t41 * qJ(5) + t80;
t61 = t24 * t95 + t104 / 0.2e1;
t37 = t41 * mrSges(6,3);
t65 = -t43 * mrSges(6,1) + t102 + t37;
t5 = -t58 + t61 + t65;
t86 = t5 * qJD(2);
t11 = t103 * t43;
t83 = qJD(2) * t11;
t81 = t23 * qJD(2);
t63 = t92 * mrSges(4,1) + t93 * mrSges(4,2);
t1 = -pkin(2) * t63 + t21 * t37 + (-t21 * mrSges(6,1) - mrSges(5,2) * t80 - t99 * t43) * t43 + (mrSges(5,1) * t80 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t43 + t99 * t41) * t41 + (-Ifges(4,2) + Ifges(4,1)) * t93 * t92 + (-t92 ^ 2 + t93 ^ 2) * Ifges(4,4) + (t102 + t104) * t54 + t103 * t24;
t70 = t1 * qJD(2);
t4 = (t41 ^ 2 + t43 ^ 2) * (mrSges(5,3) + mrSges(6,2)) + (m(5) + m(6)) * (-t30 * t43 - t97 * t41);
t68 = qJD(2) * t4;
t67 = qJD(3) * t47;
t12 = 0.2e1 * t97 * t95 - t91;
t8 = t58 + t61;
t2 = [0, 0, (t60 + t66 - t63 - t65) * qJD(3) - t98, 0, -qJD(3) * t23; 0, qJD(3) * t1 + qJD(4) * t4 + qJD(5) * t11, t8 * qJD(4) + t12 * qJD(5) + t70 + (-mrSges(4,1) * t57 + mrSges(4,2) * t79 + Ifges(4,5) * t93 - Ifges(4,6) * t92 - t53 * t91 - (t84 * t94 - mrSges(5,2) + t47) * t30 + (m(6) * t53 - t85 * t94 - mrSges(5,1) - mrSges(6,1)) * t97 + (mrSges(5,3) * t78 - Ifges(6,4) - Ifges(5,5)) * t41 + (t48 * mrSges(6,2) + mrSges(5,3) * t77 + Ifges(5,6) - Ifges(6,6)) * t43) * qJD(3), qJD(3) * t8 + t68, qJD(3) * t12 + t83; 0, -qJD(4) * t5 - t70, t47 * qJD(5), -t86, t67; 0, qJD(3) * t5 - t68 + t98, t86, 0, t81; 0, -qJD(4) * t23 - t83, -t67, -t81, 0;];
Cq = t2;
