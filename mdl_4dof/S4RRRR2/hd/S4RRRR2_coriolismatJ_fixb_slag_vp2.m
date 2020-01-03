% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:09
% EndTime: 2019-12-31 17:23:11
% DurationCPUTime: 0.85s
% Computational Cost: add. (2951->122), mult. (6044->164), div. (0->0), fcn. (5422->6), ass. (0->80)
t103 = sin(qJ(4));
t106 = cos(qJ(4));
t104 = sin(qJ(3));
t105 = sin(qJ(2));
t96 = pkin(1) * t105 + pkin(6);
t167 = pkin(7) + t96;
t78 = t167 * t104;
t107 = cos(qJ(3));
t79 = t167 * t107;
t124 = -t103 * t79 - t106 * t78;
t46 = -t103 * t78 + t106 * t79;
t187 = -t46 * mrSges(5,1) - t124 * mrSges(5,2);
t168 = -pkin(7) - pkin(6);
t92 = t168 * t104;
t93 = t168 * t107;
t123 = t103 * t93 + t106 * t92;
t60 = t103 * t92 - t106 * t93;
t186 = -t60 * mrSges(5,1) - t123 * mrSges(5,2);
t185 = m(5) / 0.2e1;
t81 = -t103 * t104 + t106 * t107;
t82 = -t103 * t107 - t104 * t106;
t143 = Ifges(5,5) * t81 + Ifges(5,6) * t82;
t10 = t143 + t187;
t184 = t10 * qJD(4);
t14 = t143 + t186;
t183 = t14 * qJD(4);
t182 = t104 ^ 2 + t107 ^ 2;
t108 = cos(qJ(2));
t163 = t108 * pkin(1);
t63 = t82 * t163;
t64 = t81 * t163;
t144 = t63 * mrSges(5,1) / 0.2e1 - t64 * mrSges(5,2) / 0.2e1;
t175 = t144 + (t103 * t64 + t106 * t63) * pkin(3) * t185;
t149 = t82 * mrSges(5,3);
t174 = (t46 + t60) * t149;
t100 = Ifges(4,4) * t107;
t48 = -mrSges(5,1) * t81 - mrSges(5,2) * t82;
t113 = -Ifges(4,4) * t104 + pkin(3) * t48 + (Ifges(4,1) - Ifges(4,2)) * t107;
t171 = t107 * t100 + t113 * t104;
t170 = t182 * t108;
t89 = mrSges(4,1) * t104 + mrSges(4,2) * t107;
t169 = -mrSges(4,1) * t107 + mrSges(4,2) * t104;
t122 = Ifges(5,4) * t81 ^ 2 + (-Ifges(5,4) * t82 + (Ifges(5,2) - Ifges(5,1)) * t81) * t82;
t166 = m(5) * t105;
t164 = pkin(3) * t104;
t152 = t81 * mrSges(5,3);
t47 = -mrSges(5,1) * t82 + mrSges(5,2) * t81;
t98 = -pkin(3) * t107 - pkin(2);
t86 = t98 - t163;
t146 = t86 * t47;
t145 = t98 * t47;
t142 = pkin(3) * qJD(3);
t110 = t122 + t171;
t77 = t86 * t164;
t97 = -pkin(2) - t163;
t5 = m(5) * t77 + t97 * t89 + t110 + t146;
t137 = t5 * qJD(1);
t7 = t122 + t146;
t136 = t7 * qJD(1);
t111 = (mrSges(4,3) * t182 - mrSges(3,2)) * t108 + (-mrSges(3,1) + t48 + t169) * t105;
t9 = m(5) * (t124 * t63 + t46 * t64) + (t63 * t82 + t64 * t81) * mrSges(5,3) + (t86 * t166 + m(4) * (t105 * t97 + t170 * t96) + t111) * pkin(1);
t135 = t9 * qJD(1);
t131 = t106 * t152;
t129 = t163 / 0.2e1;
t116 = t122 + t146 / 0.2e1 + t145 / 0.2e1 + t174 / 0.2e1;
t87 = t98 * t164;
t109 = (-pkin(2) / 0.2e1 + t97 / 0.2e1) * t89 + (-t60 / 0.2e1 - t46 / 0.2e1) * t149 + (t77 + t87) * t185 + t116;
t2 = t109 + (mrSges(4,1) * t129 + t113) * t104 + (mrSges(4,2) * t129 + t100) * t107 - t175;
t6 = m(5) * t87 - pkin(2) * t89 + t110 + t145;
t121 = t2 * qJD(1) + t6 * qJD(2);
t3 = (-t86 / 0.2e1 - t98 / 0.2e1) * t47 + t144 - t122;
t8 = t122 + t145;
t120 = -qJD(1) * t3 + qJD(2) * t8;
t119 = pkin(3) * t103 * t149 + Ifges(4,5) * t107 - Ifges(4,6) * t104 + t143;
t85 = (mrSges(5,1) * t103 + mrSges(5,2) * t106) * pkin(3);
t114 = t85 * qJD(3);
t80 = t85 * qJD(4);
t4 = t116 + t144 - t174 / 0.2e1;
t1 = t109 - t89 * t163 / 0.2e1 + t171 + t175;
t11 = [qJD(2) * t9 + qJD(3) * t5 + qJD(4) * t7, t1 * qJD(3) + t4 * qJD(4) + t135 + (m(5) * (t123 * t63 + t60 * t64) + t64 * t152 + t63 * t149 + (t98 * t166 + m(4) * (-pkin(2) * t105 + pkin(6) * t170) + t111) * pkin(1)) * qJD(2), t137 + t1 * qJD(2) + (t169 * t96 + t119 + t187) * qJD(3) + t184 + (-t131 + m(5) * (t103 * t124 - t106 * t46)) * t142, t4 * qJD(2) + t10 * qJD(3) + t136 + t184; qJD(3) * t2 - qJD(4) * t3 - t135, qJD(3) * t6 + qJD(4) * t8, (pkin(6) * t169 + t119 + t186) * qJD(3) + t183 + (-t131 + m(5) * (t103 * t123 - t106 * t60)) * t142 + t121, t14 * qJD(3) + t120 + t183; -qJD(2) * t2 - t137, -t121, -t80, -t114 - t80; qJD(2) * t3 - t136, -t120, t114, 0;];
Cq = t11;
