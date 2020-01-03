% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:58
% EndTime: 2019-12-31 19:51:00
% DurationCPUTime: 1.14s
% Computational Cost: add. (3622->177), mult. (7423->216), div. (0->0), fcn. (6999->6), ass. (0->114)
t219 = mrSges(5,3) + mrSges(6,2);
t130 = sin(pkin(8));
t131 = cos(pkin(8));
t132 = sin(qJ(4));
t189 = cos(qJ(4));
t115 = t189 * t130 + t132 * t131;
t218 = -t115 / 0.2e1;
t166 = t130 ^ 2 + t131 ^ 2;
t217 = t166 * mrSges(4,3);
t216 = m(4) * t166 * qJ(3);
t141 = -t132 * t130 + t189 * t131;
t134 = cos(qJ(2));
t186 = t134 * pkin(1);
t92 = t115 * t186;
t93 = t141 * t186;
t214 = t219 * (t92 * t115 + t141 * t93);
t177 = Ifges(6,5) * t141;
t213 = Ifges(6,3) * t218 - t177 / 0.2e1;
t212 = m(6) + m(5);
t168 = t141 * qJ(5);
t188 = pkin(4) * t115;
t72 = -t168 + t188;
t211 = m(6) * t72;
t209 = mrSges(6,1) + mrSges(5,1);
t164 = qJD(1) + qJD(2);
t73 = t115 * mrSges(6,1) - mrSges(6,3) * t141;
t74 = t115 * mrSges(5,1) + mrSges(5,2) * t141;
t18 = 0.2e1 * (t72 / 0.4e1 + t188 / 0.4e1 - t168 / 0.4e1) * m(6) + t74 + t73;
t206 = t164 * t18;
t64 = m(6) * t115;
t205 = t164 * t64;
t195 = m(6) / 0.2e1;
t204 = m(5) / 0.2e1 + t195;
t126 = m(6) * qJ(5) + mrSges(6,3);
t133 = sin(qJ(2));
t187 = t133 * pkin(1);
t156 = qJ(3) + t187;
t203 = t166 * t156;
t202 = -m(6) * pkin(4) - t209;
t201 = mrSges(5,2) - t126;
t200 = (-pkin(4) * t92 + qJ(5) * t93) * t195 + (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t93;
t197 = m(4) / 0.2e1;
t191 = t141 / 0.2e1;
t123 = -t131 * pkin(3) - pkin(2);
t150 = -pkin(4) * t141 - t115 * qJ(5);
t63 = t123 + t150;
t55 = t63 - t186;
t185 = t55 * t73;
t184 = t63 * t73;
t183 = t18 * qJD(4) - t64 * qJD(5);
t181 = t55 + t63;
t180 = t64 * qJD(3);
t179 = t18 * qJD(3);
t178 = Ifges(5,4) * t115;
t116 = t123 - t186;
t176 = t116 * t74;
t175 = t123 * t74;
t75 = -mrSges(6,1) * t141 - t115 * mrSges(6,3);
t30 = t72 * t75;
t106 = Ifges(6,5) * t115;
t78 = -Ifges(6,3) * t141 + t106;
t109 = Ifges(5,4) * t141;
t79 = -Ifges(5,2) * t115 + t109;
t80 = Ifges(5,2) * t141 + t178;
t81 = Ifges(6,1) * t141 + t106;
t82 = Ifges(6,1) * t115 - t177;
t83 = Ifges(5,1) * t141 - t178;
t84 = Ifges(5,1) * t115 + t109;
t139 = t30 + t141 * t213 + t80 * t218 + (t79 + t82 + t84) * t191 + (t78 + t81 + t83) * t115 / 0.2e1;
t3 = t55 * t211 + t139 + t176 + t185;
t174 = t3 * qJD(1);
t136 = (-mrSges(3,2) + t217) * t134 + (-t131 * mrSges(4,1) - mrSges(5,1) * t141 + t130 * mrSges(4,2) + t115 * mrSges(5,2) - mrSges(3,1) + t75) * t133;
t145 = (-pkin(7) - t156) * t130;
t127 = t131 * pkin(7);
t99 = t131 * t156 + t127;
t61 = t132 * t99 - t189 * t145;
t62 = t132 * t145 + t189 * t99;
t160 = t61 * t92 + t62 * t93;
t5 = t136 * pkin(1) + m(6) * (t55 * t187 + t160) + m(5) * (t116 * t187 + t160) + m(4) * ((-pkin(2) - t186) * t187 + t186 * t203) + t214;
t173 = t5 * qJD(1);
t138 = t217 + t219 * (t115 ^ 2 + t141 ^ 2);
t37 = t62 * t141;
t9 = t138 + m(4) * t203 + t212 * (t61 * t115 + t37);
t172 = t9 * qJD(1);
t48 = t115 * t75;
t20 = -t55 * t64 - t48;
t167 = t20 * qJD(1);
t165 = t126 * qJD(4);
t163 = t92 * t195;
t157 = (-pkin(7) - qJ(3)) * t130;
t152 = t181 * t211;
t151 = 0.2e1 * mrSges(6,2) * t191;
t2 = -t30 + (-mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t92 + (-t123 / 0.2e1 - t116 / 0.2e1) * t74 + (-t63 / 0.2e1 - t55 / 0.2e1) * t73 - t152 / 0.2e1 + (-t78 / 0.2e1 + t80 / 0.2e1 - t81 / 0.2e1 - t83 / 0.2e1) * t115 - (t213 + t84 / 0.2e1 + t82 / 0.2e1 + t79 / 0.2e1) * t141 + t200;
t4 = t63 * t211 + t139 + t175 + t184;
t149 = -t2 * qJD(1) + t4 * qJD(2);
t118 = t131 * qJ(3) + t127;
t86 = t189 * t118 + t132 * t157;
t54 = t86 * t141;
t85 = t132 * t118 - t189 * t157;
t10 = t216 + t138 + t212 * (t85 * t115 + t54);
t135 = t138 + t166 * t197 * (0.2e1 * qJ(3) + t187) + t204 * (t37 + t54 + (t61 + t85) * t115);
t6 = -t135 + (m(4) + t212) * t187 / 0.2e1;
t147 = t6 * qJD(1) - t10 * qJD(2);
t142 = -t48 - t181 * t64 / 0.2e1;
t14 = t163 - t142;
t21 = -t63 * t64 - t48;
t146 = t14 * qJD(1) - t21 * qJD(2);
t137 = t150 * mrSges(6,2) + (-Ifges(5,6) + Ifges(6,6)) * t115 + (Ifges(6,4) + Ifges(5,5)) * t141;
t29 = m(6) * t86 + t151;
t23 = m(6) * t62 + t151;
t15 = t163 + t142;
t7 = (t197 + t204) * t187 + t135;
t1 = t139 + t152 / 0.2e1 + t184 / 0.2e1 + t185 / 0.2e1 + t176 / 0.2e1 + t175 / 0.2e1 - t209 * t92 / 0.2e1 + t200;
t8 = [t5 * qJD(2) + t9 * qJD(3) + t3 * qJD(4) + t20 * qJD(5), t7 * qJD(3) + t1 * qJD(4) + t15 * qJD(5) + t173 + (0.2e1 * t204 * (t85 * t92 + t86 * t93) + (t136 + t134 * t216 + (-m(4) * pkin(2) + m(5) * t123 + m(6) * t63) * t133) * pkin(1) + t214) * qJD(2), t7 * qJD(2) + t172, t174 + t1 * qJD(2) + (t201 * t61 + t202 * t62 + t137) * qJD(4) + t23 * qJD(5), t15 * qJD(2) + t23 * qJD(4) + t167; -t6 * qJD(3) - t2 * qJD(4) - t14 * qJD(5) - t173, t10 * qJD(3) + t4 * qJD(4) + t21 * qJD(5), -t147, (t201 * t85 + t202 * t86 + t137) * qJD(4) + t29 * qJD(5) + t149, t29 * qJD(4) - t146; t6 * qJD(2) - t172 + t183, t147 + t183, 0, t206, -t205; t2 * qJD(2) - t174 - t179, -t149 - t179, -t206, t126 * qJD(5), t165; t14 * qJD(2) - t167 + t180, t146 + t180, t205, -t165, 0;];
Cq = t8;
