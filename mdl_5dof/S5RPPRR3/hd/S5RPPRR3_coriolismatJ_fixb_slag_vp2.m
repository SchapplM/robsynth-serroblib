% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:33
% EndTime: 2020-01-03 11:27:37
% DurationCPUTime: 0.94s
% Computational Cost: add. (4031->85), mult. (7398->121), div. (0->0), fcn. (8216->8), ass. (0->53)
t74 = sin(pkin(9));
t75 = cos(pkin(9));
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t67 = -t74 * t78 + t75 * t80;
t68 = t74 * t80 + t75 * t78;
t77 = sin(qJ(5));
t79 = cos(qJ(5));
t58 = t67 * t77 + t68 * t79;
t91 = t67 * t79 - t68 * t77;
t141 = t58 * mrSges(6,1) + t91 * mrSges(6,2);
t140 = t141 * qJD(5);
t71 = sin(pkin(8)) * pkin(1) + qJ(3);
t118 = pkin(6) + t71;
t63 = t118 * t74;
t64 = t118 * t75;
t42 = -t63 * t78 + t64 * t80;
t38 = pkin(7) * t67 + t42;
t41 = -t63 * t80 - t64 * t78;
t84 = -pkin(7) * t68 + t41;
t125 = -t38 * t77 + t79 * t84;
t28 = t38 * t79 + t77 * t84;
t4 = -t28 * mrSges(6,1) - mrSges(6,2) * t125 + Ifges(6,5) * t91 - Ifges(6,6) * t58;
t138 = t4 * qJD(5);
t117 = Ifges(6,4) * t58;
t120 = t58 / 0.2e1;
t121 = -t58 / 0.2e1;
t122 = t91 / 0.2e1;
t85 = -cos(pkin(8)) * pkin(1) - pkin(2) - t75 * pkin(3);
t59 = -pkin(4) * t67 + t85;
t3 = t59 * t141 + (0.2e1 * Ifges(6,4) * t91 + (Ifges(6,1) - Ifges(6,2)) * t58) * t122 + (Ifges(6,2) * t91 + t117) * t121 + (Ifges(6,1) * t91 - t117) * t120;
t124 = t67 ^ 2;
t123 = m(6) * pkin(4);
t119 = t68 / 0.2e1;
t108 = t58 * t79;
t105 = t91 * t77;
t97 = t3 * qJD(1);
t93 = mrSges(5,1) * t68 + t67 * mrSges(5,2);
t83 = -t141 - t93;
t16 = (t119 + t108 / 0.2e1 - t105 / 0.2e1) * t123 - t83;
t95 = t16 * qJD(1);
t94 = t141 * qJD(1);
t1 = t85 * t93 + t124 * Ifges(5,4) + (-Ifges(5,4) * t68 + (-Ifges(5,2) + Ifges(5,1)) * t67 + (m(6) * t59 - mrSges(6,1) * t91 + t58 * mrSges(6,2)) * pkin(4)) * t68 + t3;
t88 = t1 * qJD(1);
t7 = (t58 ^ 2 + t91 ^ 2) * mrSges(6,3) + (t68 ^ 2 + t124) * mrSges(5,3) + m(6) * (-t125 * t58 + t28 * t91) + m(5) * (-t41 * t68 + t42 * t67) + (m(4) * t71 + mrSges(4,3)) * (t74 ^ 2 + t75 ^ 2);
t86 = t7 * qJD(1);
t82 = (t105 - t108) * t123;
t5 = (t121 + t120) * Ifges(6,6) + (t122 - t91 / 0.2e1) * Ifges(6,5);
t70 = (mrSges(6,1) * t77 + mrSges(6,2) * t79) * pkin(4);
t81 = -qJD(1) * t5 + qJD(4) * t70;
t69 = t70 * qJD(5);
t29 = t119 * t123 + t82 / 0.2e1;
t2 = [qJD(3) * t7 + qJD(4) * t1 + qJD(5) * t3, 0, qJD(4) * t29 + t86, t29 * qJD(3) + t138 + t88 + (-t42 * mrSges(5,1) - t41 * mrSges(5,2) + Ifges(5,5) * t67 - Ifges(5,6) * t68 + (m(6) * (t125 * t77 - t28 * t79) + (-t58 * t77 - t79 * t91) * mrSges(6,3)) * pkin(4) + t4) * qJD(4), qJD(4) * t4 + t138 + t97; 0, 0, 0, (t82 + t83) * qJD(4) - t140, -qJD(4) * t141 - t140; qJD(4) * t16 + t140 - t86, 0, 0, t95, t94; -qJD(3) * t16 + qJD(5) * t5 - t88, 0, -t95, -t69, -t69 - t81; -qJD(3) * t141 - qJD(4) * t5 - t97, 0, -t94, t81, 0;];
Cq = t2;
