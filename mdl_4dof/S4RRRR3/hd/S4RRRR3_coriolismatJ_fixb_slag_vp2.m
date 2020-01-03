% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR3
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:18
% EndTime: 2019-12-31 17:24:22
% DurationCPUTime: 1.67s
% Computational Cost: add. (4732->133), mult. (9406->186), div. (0->0), fcn. (9740->6), ass. (0->80)
t116 = sin(qJ(3));
t117 = sin(qJ(2));
t119 = cos(qJ(3));
t120 = cos(qJ(2));
t100 = -t116 * t117 + t119 * t120;
t101 = -t116 * t120 - t119 * t117;
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t135 = t118 * t100 + t101 * t115;
t71 = t100 * t115 - t101 * t118;
t195 = Ifges(5,5) * t135 - Ifges(5,6) * t71;
t173 = -pkin(6) - pkin(5);
t106 = t173 * t117;
t107 = t173 * t120;
t184 = t119 * t106 + t116 * t107;
t192 = t101 * pkin(7) + t184;
t80 = t106 * t116 - t119 * t107;
t56 = pkin(7) * t100 + t80;
t204 = -t115 * t56 + t118 * t192;
t209 = t204 * mrSges(5,2);
t34 = t115 * t192 + t118 * t56;
t212 = t34 * mrSges(5,1);
t215 = -t212 / 0.2e1 - t209 / 0.2e1;
t218 = 0.2e1 * t215 + t195;
t220 = t218 * qJD(4);
t113 = -pkin(2) * t120 - pkin(1);
t172 = Ifges(5,4) * t71;
t176 = -t71 / 0.2e1;
t81 = -t100 * pkin(3) + t113;
t5 = (Ifges(5,1) * t135 - t172) * t71 / 0.2e1 + (Ifges(5,2) * t135 + t172) * t176 + (t81 * mrSges(5,2) + Ifges(5,4) * t135) * t135 + (t81 * mrSges(5,1) + (Ifges(5,1) - Ifges(5,2)) * t135 / 0.2e1) * t71;
t219 = t113 * (-mrSges(4,1) * t101 + mrSges(4,2) * t100) + t5;
t201 = -t80 * mrSges(4,1) - t184 * mrSges(4,2) + Ifges(4,5) * t100 + Ifges(4,6) * t101 + t195;
t205 = -t209 - t212;
t217 = t201 + t205;
t216 = m(5) * t81 - mrSges(5,1) * t135 + mrSges(5,2) * t71;
t211 = t115 * t204 - t118 * t34;
t180 = m(5) * pkin(3);
t191 = -t211 * t180 / 0.2e1 - t215;
t181 = -pkin(3) / 0.2e1;
t140 = t115 * t116;
t94 = (t118 * t119 - t140) * pkin(2);
t152 = t94 * mrSges(5,2);
t139 = t116 * t118;
t93 = (-t115 * t119 - t139) * pkin(2);
t91 = t93 * mrSges(5,1);
t182 = (mrSges(4,1) * t116 + mrSges(4,2) * t119) * pkin(2) - t91 + t152;
t174 = t94 / 0.2e1;
t169 = pkin(3) * t101;
t167 = t117 * pkin(2);
t112 = pkin(2) * t119 + pkin(3);
t89 = -pkin(2) * t140 + t112 * t118;
t154 = t89 * t135;
t90 = pkin(2) * t139 + t112 * t115;
t153 = t90 * t71;
t133 = Ifges(4,4) * t100 + (-Ifges(4,1) + Ifges(4,2)) * t101;
t1 = (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t117) * t117 + m(4) * t113 * t167 + (-mrSges(4,2) * t167 - Ifges(4,4) * t101) * t101 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t120 + (Ifges(3,1) - Ifges(3,2)) * t117) * t120 + (-mrSges(4,1) * t167 + t133) * t100 + t216 * (t167 - t169) + t219;
t149 = t1 * qJD(1);
t147 = t115 * t71;
t145 = t118 * t135;
t2 = -Ifges(4,4) * t101 ^ 2 + t133 * t100 - t216 * t169 + t219;
t144 = t2 * qJD(1);
t143 = t5 * qJD(1);
t18 = -t90 * mrSges(5,1) - t89 * mrSges(5,2);
t142 = qJD(4) * t18;
t138 = t118 * t181;
t134 = t138 - t89 / 0.2e1;
t17 = -m(5) * (t89 * t93 + t90 * t94) + t182;
t122 = m(5) * ((t90 + t93) * t204 + (-t89 + t94) * t34) / 0.2e1 + t215;
t124 = -t153 / 0.2e1 - t154 / 0.2e1 + t135 * t174 + t93 * t176;
t4 = ((t147 / 0.2e1 + t145 / 0.2e1) * pkin(3) + t124) * mrSges(5,3) + t122 + t191;
t131 = t4 * qJD(1) - t17 * qJD(2);
t130 = t18 * qJD(2);
t129 = (t115 * t181 - t90 / 0.2e1) * mrSges(5,1);
t105 = (t115 * mrSges(5,1) + t118 * mrSges(5,2)) * pkin(3);
t14 = -t91 / 0.2e1 + t129 + (t174 + t134) * mrSges(5,2);
t126 = -qJD(2) * t14 + qJD(3) * t105;
t99 = t105 * qJD(4);
t15 = -t152 / 0.2e1 + t91 / 0.2e1 + t134 * mrSges(5,2) + t129;
t3 = t122 - t191 + t201 + (t138 * t135 + t147 * t181 + t124) * mrSges(5,3);
t6 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t5, t3 * qJD(3) + t220 + t149 + (m(5) * (t204 * t90 - t89 * t34) + Ifges(3,5) * t120 - Ifges(3,6) * t117 + (m(4) * (t116 * t184 - t119 * t80) + (-t100 * t119 + t101 * t116) * mrSges(4,3)) * pkin(2) + (-t120 * mrSges(3,1) + t117 * mrSges(3,2)) * pkin(5) + (-t153 - t154) * mrSges(5,3) + t217) * qJD(2), t3 * qJD(2) + t220 + t144 + ((m(5) * t211 + (-t145 - t147) * mrSges(5,3)) * pkin(3) + t217) * qJD(3), t143 + (t195 + t205) * qJD(4) + (qJD(2) + qJD(3)) * t218; qJD(3) * t4 - t149, -qJD(3) * t17 + t142, ((t115 * t94 + t118 * t93) * t180 - t182) * qJD(3) + t15 * qJD(4) + t131, t15 * qJD(3) + t130 + t142; -qJD(2) * t4 - t144, qJD(4) * t14 - t131, -t99, -t126 - t99; -t143, -qJD(3) * t14 - t130, t126, 0;];
Cq = t6;
