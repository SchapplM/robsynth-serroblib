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
% m [6x1]
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:14:10
% EndTime: 2022-01-23 09:14:12
% DurationCPUTime: 0.83s
% Computational Cost: add. (4031->84), mult. (7398->122), div. (0->0), fcn. (8216->8), ass. (0->51)
t75 = sin(pkin(9));
t76 = cos(pkin(9));
t79 = sin(qJ(4));
t81 = cos(qJ(4));
t68 = -t79 * t75 + t81 * t76;
t69 = t81 * t75 + t79 * t76;
t78 = sin(qJ(5));
t80 = cos(qJ(5));
t59 = t78 * t68 + t80 * t69;
t91 = t80 * t68 - t78 * t69;
t140 = t59 * mrSges(6,1) + t91 * mrSges(6,2);
t141 = t140 * qJD(5);
t72 = sin(pkin(8)) * pkin(1) + qJ(3);
t116 = pkin(6) + t72;
t64 = t116 * t75;
t65 = t116 * t76;
t43 = -t79 * t64 + t81 * t65;
t39 = t68 * pkin(7) + t43;
t42 = -t81 * t64 - t79 * t65;
t84 = -t69 * pkin(7) + t42;
t124 = -t78 * t39 + t80 * t84;
t28 = t80 * t39 + t78 * t84;
t4 = -t28 * mrSges(6,1) - t124 * mrSges(6,2) + Ifges(6,5) * t91 - Ifges(6,6) * t59;
t139 = t4 * qJD(5);
t115 = Ifges(6,4) * t59;
t119 = t59 / 0.2e1;
t120 = -t59 / 0.2e1;
t121 = t91 / 0.2e1;
t85 = -cos(pkin(8)) * pkin(1) - pkin(2) - t76 * pkin(3);
t60 = -t68 * pkin(4) + t85;
t3 = t60 * t140 + (0.2e1 * Ifges(6,4) * t91 + (Ifges(6,1) - Ifges(6,2)) * t59) * t121 + (Ifges(6,2) * t91 + t115) * t120 + (Ifges(6,1) * t91 - t115) * t119;
t123 = t69 ^ 2;
t32 = (-t59 * t80 + t78 * t91) * pkin(4);
t118 = m(6) * t32;
t117 = t69 * pkin(4);
t97 = t3 * qJD(1);
t93 = t69 * mrSges(5,1) + t68 * mrSges(5,2);
t83 = -t140 - t93;
t16 = (t32 / 0.2e1 - t117 / 0.2e1) * m(6) + t83;
t95 = t16 * qJD(1);
t94 = t140 * qJD(1);
t1 = t85 * t93 - t123 * Ifges(5,4) + (Ifges(5,4) * t68 + (Ifges(5,1) - Ifges(5,2)) * t69) * t68 + (m(6) * t60 - mrSges(6,1) * t91 + t59 * mrSges(6,2)) * t117 + t3;
t88 = t1 * qJD(1);
t7 = (t59 ^ 2 + t91 ^ 2) * mrSges(6,3) + (t68 ^ 2 + t123) * mrSges(5,3) + m(6) * (-t124 * t59 + t28 * t91) + m(5) * (-t42 * t69 + t43 * t68) + (m(4) * t72 + mrSges(4,3)) * (t75 ^ 2 + t76 ^ 2);
t86 = t7 * qJD(1);
t5 = (t120 + t119) * Ifges(6,6) + (t121 - t91 / 0.2e1) * Ifges(6,5);
t71 = (mrSges(6,1) * t78 + mrSges(6,2) * t80) * pkin(4);
t82 = -t5 * qJD(1) + t71 * qJD(4);
t70 = t71 * qJD(5);
t29 = m(6) * t117 / 0.2e1 + t118 / 0.2e1;
t2 = [t7 * qJD(3) + t1 * qJD(4) + t3 * qJD(5), 0, t29 * qJD(4) + t86, t29 * qJD(3) + t139 + t88 + (-t43 * mrSges(5,1) - t42 * mrSges(5,2) + Ifges(5,5) * t68 - Ifges(5,6) * t69 + (m(6) * (t124 * t78 - t28 * t80) + (-t59 * t78 - t80 * t91) * mrSges(6,3)) * pkin(4) + t4) * qJD(4), t4 * qJD(4) + t139 + t97; 0, 0, 0, (t83 + t118) * qJD(4) - t141, -qJD(4) * t140 - t141; -t16 * qJD(4) + t141 - t86, 0, 0, -t95, t94; t16 * qJD(3) + t5 * qJD(5) - t88, 0, t95, -t70, -t70 - t82; -qJD(3) * t140 - t5 * qJD(4) - t97, 0, -t94, t82, 0;];
Cq = t2;
