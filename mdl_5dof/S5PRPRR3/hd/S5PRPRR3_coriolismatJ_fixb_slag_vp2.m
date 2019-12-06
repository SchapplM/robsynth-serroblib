% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR3
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:41
% EndTime: 2019-12-05 15:46:43
% DurationCPUTime: 0.88s
% Computational Cost: add. (2617->96), mult. (5738->149), div. (0->0), fcn. (6110->8), ass. (0->70)
t151 = cos(qJ(5));
t94 = sin(qJ(5));
t95 = sin(qJ(4));
t96 = cos(qJ(4));
t104 = t151 * t95 + t94 * t96;
t164 = t151 * t96 - t94 * t95;
t120 = sin(pkin(9));
t121 = cos(pkin(9));
t150 = sin(qJ(2));
t152 = cos(qJ(2));
t77 = t120 * t150 - t121 * t152;
t48 = t104 * t77;
t50 = t164 * t77;
t26 = -t104 * t50 + t164 * t48;
t179 = m(6) * t26;
t167 = t121 * pkin(2);
t88 = -pkin(3) - t167;
t83 = -t96 * pkin(4) + t88;
t178 = m(6) * t83 - mrSges(6,1) * t164 + mrSges(6,2) * t104;
t168 = t120 * pkin(2);
t87 = pkin(6) + t168;
t153 = pkin(7) + t87;
t71 = t153 * t95;
t72 = t153 * t96;
t109 = -t151 * t71 - t94 * t72;
t58 = t151 * t72 - t94 * t71;
t11 = -t58 * mrSges(6,1) - t109 * mrSges(6,2) + Ifges(6,5) * t164 - Ifges(6,6) * t104;
t177 = t11 * qJD(5);
t137 = t95 * mrSges(5,1);
t106 = t96 * mrSges(5,2) + t137;
t173 = t77 * t106;
t60 = mrSges(6,1) * t104 + mrSges(6,2) * t164;
t172 = t60 * qJD(5);
t78 = -t120 * t152 - t121 * t150;
t47 = t164 * t78;
t49 = t104 * t78;
t110 = t47 * mrSges(6,1) - t49 * mrSges(6,2);
t170 = t110 * qJD(5);
t7 = -Ifges(6,4) * t104 ^ 2 + (Ifges(6,4) * t164 - (-Ifges(6,1) + Ifges(6,2)) * t104) * t164 + t83 * t60;
t92 = t95 ^ 2;
t123 = t96 ^ 2 + t92;
t163 = -t96 * mrSges(5,1) + t95 * mrSges(5,2);
t127 = t48 * mrSges(6,1) / 0.2e1 + t50 * mrSges(6,2) / 0.2e1;
t160 = -m(6) / 0.2e1;
t158 = m(6) * pkin(4);
t155 = pkin(4) * t94;
t154 = pkin(4) * t95;
t142 = t77 * t60;
t62 = t77 * t78;
t141 = t164 * mrSges(6,3);
t140 = t104 * mrSges(6,3);
t117 = t142 / 0.2e1;
t116 = qJD(1) * t179;
t115 = pkin(4) * t151;
t82 = (mrSges(6,1) * t94 + mrSges(6,2) * t151) * pkin(4);
t105 = t82 * qJD(4);
t8 = m(6) * (t47 * t50 + t49 * t48 - t62) + m(5) * (t123 - 0.1e1) * t62;
t101 = qJD(3) * t160 * t26 - t8 * qJD(1);
t97 = t77 * t154 * t160 - t173 / 0.2e1;
t98 = (t151 * t48 - t50 * t94) * t158 / 0.2e1 + t127 + t173 / 0.2e1;
t1 = -t142 / 0.2e1 + t97 + t98;
t3 = t88 * t137 - t92 * Ifges(5,4) + (t88 * mrSges(5,2) + Ifges(5,4) * t96 + (Ifges(5,1) - Ifges(5,2)) * t95) * t96 + t7 + t178 * t154;
t100 = t1 * qJD(1) - t3 * qJD(2);
t4 = t117 - t127;
t99 = -t4 * qJD(1) - t7 * qJD(2);
t79 = t82 * qJD(5);
t6 = qJD(2) * t179 / 0.2e1;
t5 = t117 + t127;
t2 = -t97 + t98 + t117;
t9 = [t8 * qJD(2), (m(6) * (t109 * t48 - t58 * t50) - t50 * t141 - t48 * t140 - t152 * mrSges(3,2) - t150 * mrSges(3,1) + (m(4) * t167 - m(5) * t88 + mrSges(4,1) - t163 - t178) * t78 + (-m(4) * t168 + mrSges(4,2) + (-m(5) * t87 - mrSges(5,3)) * t123) * t77) * qJD(2) + t2 * qJD(4) + t5 * qJD(5) - t101, t6, t2 * qJD(2) + (m(6) * (t47 * t115 + t49 * t155) - t163 * t78 + t110) * qJD(4) + t170, t5 * qJD(2) + qJD(4) * t110 + t170; -t1 * qJD(4) + t4 * qJD(5) + t101, t3 * qJD(4) + t7 * qJD(5), -t116 / 0.2e1, (-t140 * t155 - t115 * t141 + (t109 * t94 - t151 * t58) * t158 + Ifges(5,5) * t96 - Ifges(5,6) * t95 + t163 * t87 + t11) * qJD(4) + t177 - t100, t11 * qJD(4) + t177 - t99; t6, t116 / 0.2e1, 0, -t172 + (-t106 - t60 + (-t104 * t115 + t155 * t164) * m(6)) * qJD(4), -qJD(4) * t60 - t172; t1 * qJD(2), t100, 0, -t79, -t105 - t79; -t4 * qJD(2), t99, 0, t105, 0;];
Cq = t9;
