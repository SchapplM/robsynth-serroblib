% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:09
% EndTime: 2019-12-31 18:47:11
% DurationCPUTime: 1.07s
% Computational Cost: add. (1904->168), mult. (3500->210), div. (0->0), fcn. (2272->4), ass. (0->105)
t226 = mrSges(6,2) + mrSges(5,3);
t131 = cos(qJ(4));
t127 = t131 ^ 2;
t129 = sin(qJ(4));
t169 = t129 ^ 2 + t127;
t130 = sin(qJ(3));
t121 = t130 * qJ(2);
t132 = cos(qJ(3));
t133 = -pkin(1) - pkin(2);
t97 = t132 * t133 - t121;
t227 = t169 * t97;
t228 = -t97 * mrSges(4,2) + t226 * t227;
t224 = Ifges(5,4) - Ifges(6,5);
t184 = t131 * mrSges(5,2);
t187 = t129 * mrSges(5,1);
t105 = t184 + t187;
t139 = -t224 * t127 + (t224 * t129 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,3)) * t131) * t129;
t222 = -pkin(3) * t105 - t139;
t221 = t169 * t132;
t101 = t131 * mrSges(6,1) + t129 * mrSges(6,3);
t92 = t130 * t101;
t102 = -t131 * mrSges(5,1) + t129 * mrSges(5,2);
t93 = t130 * t102;
t215 = t93 - t92;
t203 = m(6) / 0.2e1;
t152 = -t131 * pkin(4) - t129 * qJ(5);
t99 = -pkin(3) + t152;
t45 = -t97 - t99;
t214 = -t101 + (-t45 + t99) * t203;
t213 = t102 - t101;
t183 = t131 * mrSges(6,3);
t186 = t129 * mrSges(6,1);
t104 = t183 - t186;
t178 = qJ(5) * t131;
t196 = t129 * pkin(4);
t103 = t178 - t196;
t147 = m(6) * t103;
t212 = t104 + t147;
t209 = t92 / 0.2e1 - t93 / 0.2e1;
t208 = m(6) * t45 + t101;
t207 = (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t131;
t167 = qJD(4) * t131;
t206 = -(qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * qJD(4) * t129 - (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t167;
t205 = 2 * qJD(3);
t204 = m(5) / 0.2e1;
t200 = t97 / 0.2e1;
t199 = -t105 / 0.2e1;
t195 = t227 * pkin(7);
t194 = m(6) * qJD(5);
t185 = t129 * t97;
t153 = pkin(3) - t97;
t172 = t130 * t133;
t98 = t132 * qJ(2) + t172;
t96 = -pkin(7) + t98;
t5 = (-m(5) - m(6)) * t96 * t227 + (-m(5) * t153 - mrSges(4,1) + t102 - t208) * t98 + t228;
t181 = t5 * qJD(1);
t170 = -t130 * mrSges(4,1) - t132 * mrSges(4,2);
t140 = t226 * t221 + t170;
t159 = t169 * t96;
t6 = -m(5) * ((pkin(3) + t121) * t130 + (t159 - t172) * t132) - m(6) * (t130 * t45 + t132 * t159) - m(4) * (-t130 * t97 + t98 * t132) - mrSges(3,3) - m(3) * qJ(2) + t140 + t215;
t180 = t6 * qJD(1);
t7 = -t103 * t101 + t153 * t105 - t212 * t45 + t139;
t179 = t7 * qJD(1);
t25 = t208 * t129;
t177 = qJD(1) * t25;
t142 = (-mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1) * t129 + t207;
t148 = t178 / 0.2e1 - t196 / 0.2e1;
t161 = t104 / 0.2e1 + t199;
t12 = ((t103 / 0.2e1 + t148) * m(6) + t142 + t161) * t132;
t176 = t12 * qJD(1);
t173 = t129 * t132;
t171 = t221 * pkin(7);
t120 = m(6) * qJ(5) + mrSges(6,3);
t166 = t120 * qJD(4);
t165 = m(6) * t173;
t164 = t129 * t194;
t163 = t99 / 0.2e1 - t45 / 0.2e1;
t160 = m(6) * t99 - t101;
t35 = (-0.1e1 + t169) * t132 * t130;
t21 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t35;
t144 = (-t98 + t159) * t132;
t134 = (t144 + (t227 + t153) * t130) * t204 + (t144 + (t45 + t227) * t130) * t203 + t209;
t149 = (t130 * t99 + t171) * t203 + (-pkin(3) * t130 + t171) * t204;
t3 = -t134 + t140 + t149 - t209;
t151 = -t3 * qJD(1) + t21 * qJD(2);
t16 = (t101 + (t200 - t163) * m(6)) * t129;
t34 = t160 * t129;
t150 = qJD(1) * t16 + qJD(3) * t34;
t11 = ((-t103 / 0.2e1 + t148) * m(6) + t142 - t161) * t132;
t135 = t214 * t103 + t163 * t104 + t97 * t199 - t222;
t137 = -t183 / 0.2e1 + t187 / 0.2e1 + t186 / 0.2e1 - t147 / 0.2e1 + t184 / 0.2e1;
t2 = t137 * t97 + t135;
t8 = -t160 * t103 - t99 * t104 + t222;
t143 = t2 * qJD(1) - t11 * qJD(2) + t8 * qJD(3);
t141 = (t199 + t212 / 0.2e1) * t132;
t138 = (m(6) * t152 + t213) * qJD(4);
t136 = t137 * t132;
t100 = (m(6) * pkin(7) + mrSges(6,2)) * t131;
t40 = (m(6) * t96 - mrSges(6,2)) * t131;
t17 = t214 * t129 + t185 * t203;
t14 = -t136 - t141;
t13 = -t136 + t141;
t4 = (t102 / 0.2e1 - t101 / 0.2e1) * t130 + t134 + t149;
t1 = t147 * t200 + t135 + t97 * t207 - (mrSges(5,1) + mrSges(6,1)) * t185 / 0.2e1;
t9 = [-qJD(2) * t6 - qJD(3) * t5 - qJD(4) * t7 + qJD(5) * t25, -t180 + 0.2e1 * (t204 + t203) * t35 * qJD(2) + t4 * qJD(3) + t14 * qJD(4), -t181 + t4 * qJD(2) + t1 * qJD(4) + t17 * qJD(5) + ((-pkin(3) * t98 + t195) * t204 + (t98 * t99 + t195) * t203) * t205 + ((-mrSges(4,1) + t213) * t98 + t228) * qJD(3), t14 * qJD(2) + t1 * qJD(3) + t40 * qJD(5) + t96 * t138 - t179 - t206, qJD(3) * t17 + qJD(4) * t40 + t177; -qJD(3) * t3 - qJD(4) * t12 - t132 * t164 + t180, t21 * qJD(3), (t170 + t215) * qJD(3) + t13 * qJD(4) + t149 * t205 + (t226 * qJD(3) * t169 + t164) * t132 + t151, -t176 + t13 * qJD(3) + (t131 * t194 + t138) * t130, (t130 * t167 + (-qJD(1) + qJD(3)) * t173) * m(6); qJD(2) * t3 + qJD(4) * t2 - qJD(5) * t16 + t181, -qJD(4) * t11 - t151, qJD(4) * t8 - qJD(5) * t34, pkin(7) * t138 + t100 * qJD(5) + t143 + t206, qJD(4) * t100 - t150; qJD(2) * t12 - qJD(3) * t2 + t179, qJD(3) * t11 + t176, -t143, t120 * qJD(5), t166; qJD(2) * t165 + qJD(3) * t16 - t177, qJD(1) * t165, t150, -t166, 0;];
Cq = t9;
