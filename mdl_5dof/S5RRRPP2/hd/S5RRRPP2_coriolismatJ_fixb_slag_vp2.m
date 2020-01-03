% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:31
% EndTime: 2019-12-31 20:51:33
% DurationCPUTime: 0.98s
% Computational Cost: add. (1737->193), mult. (3415->216), div. (0->0), fcn. (2166->4), ass. (0->103)
t109 = sin(qJ(3));
t111 = cos(qJ(3));
t182 = t109 ^ 2 + t111 ^ 2;
t112 = cos(qJ(2));
t162 = t112 * pkin(1);
t181 = t182 * t162;
t113 = -pkin(3) - pkin(4);
t149 = t113 * t109;
t150 = t111 * qJ(4);
t180 = t149 + t150;
t179 = Ifges(4,4) - Ifges(6,4) - Ifges(5,5);
t93 = t109 * qJ(5);
t68 = pkin(7) * t109 - t93;
t73 = (pkin(7) - qJ(5)) * t111;
t178 = m(6) * (t109 * t68 + t111 * t73);
t110 = sin(qJ(2));
t163 = t110 * pkin(1);
t89 = pkin(7) + t163;
t47 = t109 * t89 - t93;
t48 = (-qJ(5) + t89) * t111;
t177 = m(6) * (t109 * t47 + t111 * t48);
t176 = mrSges(6,2) + mrSges(5,3);
t69 = -t111 * mrSges(5,1) - t109 * mrSges(5,3);
t175 = -t111 * mrSges(4,1) + t109 * mrSges(4,2) + t69;
t147 = qJD(1) + qJD(2);
t75 = -t109 * mrSges(6,1) + t111 * mrSges(6,2);
t19 = 0.2e1 * (-t150 / 0.4e1 - t149 / 0.4e1 - t180 / 0.4e1) * m(6) - t75;
t174 = t147 * t19;
t94 = t109 * qJ(4);
t173 = -t111 * pkin(3) - t94;
t119 = (Ifges(4,1) + Ifges(5,1) + Ifges(6,1) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3)) * t109 + t179 * t111;
t129 = t179 * t109;
t172 = -t129 * t109 + t119 * t111;
t171 = m(5) / 0.2e1;
t170 = -m(6) / 0.2e1;
t169 = m(6) / 0.2e1;
t90 = -pkin(2) - t162;
t45 = t90 + t173;
t168 = m(5) * t45;
t145 = pkin(2) - t173;
t167 = m(5) * t145;
t72 = pkin(3) * t109 - t150;
t166 = m(5) * t72;
t164 = m(6) * t109;
t70 = mrSges(6,1) * t111 + mrSges(6,2) * t109;
t161 = t180 * t70;
t160 = t45 - t145;
t157 = -t19 * qJD(3) + qJD(4) * t164;
t156 = t182 * mrSges(6,3);
t155 = mrSges(4,2) * t111;
t115 = t72 * t69 + t161 + t172;
t74 = mrSges(5,1) * t109 - mrSges(5,3) * t111;
t136 = t74 + t166;
t105 = t111 * pkin(4);
t38 = t105 - t45;
t26 = t38 * t180;
t76 = mrSges(4,1) * t109 + t155;
t3 = m(6) * t26 + t136 * t45 + t38 * t75 + t90 * t76 + t115;
t154 = t3 * qJD(1);
t120 = -mrSges(3,2) + t182 * (mrSges(5,2) + mrSges(4,3) - mrSges(6,3));
t133 = -t70 - mrSges(3,1) + t175;
t5 = ((m(4) * t90 - m(6) * t38 + t133 + t168) * t110 + (t120 + t177) * t112) * pkin(1) + (m(5) + m(4)) * t181 * t89;
t153 = t5 * qJD(1);
t62 = t109 * t70;
t11 = -t38 * t164 - t62 + (t69 + t168) * t109;
t152 = qJD(1) * t11;
t15 = t156 - t177;
t151 = qJD(1) * t15;
t85 = (m(5) + m(6)) * qJ(4) + t176;
t148 = t85 * qJD(3);
t142 = qJD(5) * t164;
t141 = -t163 / 0.2e1;
t138 = m(5) * t160;
t46 = t105 + t145;
t137 = m(6) * (t38 + t46);
t135 = (-t47 - t68) * t109;
t134 = (-t48 - t73) * t111;
t132 = mrSges(6,3) * t94 + Ifges(5,6) * t109 + (Ifges(5,4) + Ifges(4,5)) * t111;
t128 = (t171 + t169) * t162;
t28 = t46 * t180;
t114 = t161 + (t46 / 0.2e1 + t38 / 0.2e1) * t75 + (-t145 / 0.2e1 + t45 / 0.2e1) * t74 + (t90 / 0.2e1 - pkin(2) / 0.2e1) * t76 + (t160 * t171 + t69) * t72 + (t26 + t28) * t169;
t2 = ((mrSges(4,1) / 0.2e1 + mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1 + pkin(3) * t171 + t113 * t170) * t162 - t129) * t109 + ((mrSges(4,2) / 0.2e1 - mrSges(6,2) / 0.2e1 - mrSges(5,3) / 0.2e1 + 0.2e1 * (-m(5) / 0.4e1 - m(6) / 0.4e1) * qJ(4)) * t162 + t119) * t111 + t114;
t4 = m(6) * t28 - pkin(2) * t76 - t136 * t145 + t46 * t75 + t115;
t127 = t2 * qJD(1) + t4 * qJD(2);
t13 = -t46 * t164 - t62 + (t69 - t167) * t109;
t7 = -t62 + (t69 + t128 + t138 / 0.2e1 - t137 / 0.2e1) * t109;
t124 = -qJD(1) * t7 - qJD(2) * t13;
t20 = t156 - t178;
t9 = (t141 - t134 / 0.2e1 - t135 / 0.2e1) * m(6) - t156;
t123 = -qJD(1) * t9 + qJD(2) * t20;
t25 = (t169 + t170) * t180;
t117 = t25 * qJD(5) + ((-qJ(4) * mrSges(5,2) - Ifges(4,6) - Ifges(6,6)) * t109 + (-pkin(3) * mrSges(5,2) - t113 * mrSges(6,3) - Ifges(6,5)) * t111) * qJD(3);
t116 = qJD(3) * (m(5) * t173 + t175);
t98 = t111 * mrSges(5,2);
t66 = t147 * t164;
t29 = t98 + m(6) * t73 + (m(5) * pkin(7) - mrSges(6,3)) * t111;
t24 = t25 * qJD(3);
t21 = t98 + m(6) * t48 + (m(5) * t89 - mrSges(6,3)) * t111;
t17 = t19 * qJD(5);
t10 = (t134 + t135) * t169 + m(6) * t141 + t156;
t8 = t62 + (-t69 + t128) * t109 + (-t138 + t137) * t109 / 0.2e1;
t1 = t114 + t172 - (t155 + t166 + (mrSges(4,1) + mrSges(5,1) + mrSges(6,1)) * t109) * t162 / 0.2e1 + (m(6) * t180 + t176 * t111) * t162 / 0.2e1;
t6 = [qJD(2) * t5 + qJD(3) * t3 - qJD(4) * t11 + qJD(5) * t15, t1 * qJD(3) + t8 * qJD(4) + t10 * qJD(5) + t153 + (((-m(4) * pkin(2) - m(6) * t46 + t133 - t167) * t110 + (t120 + t178) * t112) * pkin(1) + 0.2e1 * (t171 + m(4) / 0.2e1) * t181 * pkin(7)) * qJD(2), t154 + t1 * qJD(2) + (-t48 * mrSges(6,1) - t47 * mrSges(6,2) + m(6) * (-qJ(4) * t47 + t113 * t48) + t132) * qJD(3) + t21 * qJD(4) + t89 * t116 + t117, qJD(2) * t8 + qJD(3) * t21 - t152, qJD(2) * t10 + t151 + t24; qJD(3) * t2 - qJD(4) * t7 - qJD(5) * t9 - t153, qJD(3) * t4 - qJD(4) * t13 + qJD(5) * t20, (m(6) * (-qJ(4) * t68 + t113 * t73) - t73 * mrSges(6,1) - t68 * mrSges(6,2) + t132) * qJD(3) + t29 * qJD(4) + pkin(7) * t116 + t117 + t127, qJD(3) * t29 + t124, t123 + t24; -qJD(2) * t2 - t154 + t17, -t127 + t17, t85 * qJD(4), t148, t174; qJD(2) * t7 - t142 + t152, -t124 - t142, -t148, 0, -t66; qJD(2) * t9 - t151 + t157, -t123 + t157, -t174, t66, 0;];
Cq = t6;
