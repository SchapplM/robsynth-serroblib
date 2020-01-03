% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:25
% EndTime: 2020-01-03 12:07:27
% DurationCPUTime: 1.23s
% Computational Cost: add. (2637->183), mult. (6035->245), div. (0->0), fcn. (4781->8), ass. (0->120)
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t175 = -Ifges(6,4) * t91 + (-Ifges(6,2) + Ifges(6,1)) * t94;
t176 = t175 * t91;
t78 = -t94 * mrSges(6,1) + t91 * mrSges(6,2);
t127 = -mrSges(5,1) + t78;
t93 = sin(qJ(2));
t95 = cos(qJ(3));
t131 = t93 * t95;
t96 = cos(qJ(2));
t85 = t96 * pkin(1) + pkin(2);
t92 = sin(qJ(3));
t68 = pkin(1) * t131 + t92 * t85;
t89 = sin(pkin(9));
t138 = t89 * t68;
t133 = t92 * t93;
t67 = -pkin(1) * t133 + t95 * t85;
t90 = cos(pkin(9));
t39 = t90 * t67 - t138;
t137 = t89 * t92;
t71 = (t90 * t95 - t137) * pkin(2);
t171 = t39 + t71;
t170 = (-t92 * mrSges(4,1) - t95 * mrSges(4,2)) * pkin(2);
t153 = Ifges(6,4) * t94;
t169 = t153 * t94 + t176;
t116 = -t90 * pkin(3) - pkin(4);
t163 = m(5) * pkin(3);
t168 = m(6) * t116 - t90 * t163;
t87 = t91 ^ 2;
t88 = t94 ^ 2;
t126 = t87 + t88;
t109 = t126 * mrSges(6,3) - mrSges(5,2);
t73 = (t95 * t96 - t133) * pkin(1);
t141 = t73 * mrSges(4,2);
t72 = (-t92 * t96 - t131) * pkin(1);
t142 = t72 * mrSges(4,1);
t43 = t89 * t72 + t90 * t73;
t167 = t109 * t43 - t141 + t142 + (-t93 * mrSges(3,1) - t96 * mrSges(3,2)) * pkin(1);
t83 = t89 * pkin(3) + pkin(8);
t112 = t126 * t83;
t139 = t88 * mrSges(6,3);
t140 = t87 * mrSges(6,3);
t166 = m(6) * t112 + t89 * t163 + t139 + t140;
t165 = m(5) / 0.2e1;
t164 = m(6) / 0.2e1;
t130 = t94 * mrSges(6,2);
t135 = t91 * mrSges(6,1);
t106 = t130 + t135;
t62 = pkin(3) + t67;
t36 = t90 * t62 - t138;
t30 = -pkin(4) - t36;
t20 = t30 * t106;
t162 = -t20 / 0.2e1;
t161 = -t39 / 0.2e1;
t160 = -t43 / 0.2e1;
t155 = t95 * pkin(2);
t84 = pkin(3) + t155;
t65 = -pkin(2) * t137 + t90 * t84;
t56 = -pkin(4) - t65;
t45 = t56 * t106;
t159 = -t45 / 0.2e1;
t69 = t116 * t106;
t158 = -t69 / 0.2e1;
t157 = -t71 / 0.2e1;
t52 = t90 * t68;
t38 = t89 * t67 + t52;
t152 = t38 * mrSges(5,1);
t151 = t38 * t78;
t150 = t39 * mrSges(5,2);
t42 = -t90 * t72 + t89 * t73;
t149 = t42 * mrSges(5,1);
t148 = t42 * t78;
t147 = t67 * mrSges(4,2);
t146 = t68 * mrSges(4,1);
t136 = t90 * t92;
t70 = (t89 * t95 + t136) * pkin(2);
t145 = t70 * mrSges(5,1);
t144 = t70 * t78;
t143 = t71 * mrSges(5,2);
t37 = t89 * t62 + t52;
t108 = -t146 - t147;
t31 = pkin(8) + t37;
t115 = t31 * t126;
t3 = t127 * t38 + t109 * t39 + m(6) * (t39 * t115 + t30 * t38) + m(5) * (-t36 * t38 + t37 * t39) + t108;
t125 = t3 * qJD(1);
t4 = t127 * t42 + m(6) * (t43 * t115 + t30 * t42) + m(5) * (-t36 * t42 + t37 * t43) + m(4) * (t67 * t72 + t68 * t73) + t167;
t124 = t4 * qJD(1);
t99 = -Ifges(6,4) * t88 - t176;
t15 = -t20 + t99;
t123 = t15 * qJD(1);
t122 = -pkin(2) * t92 / 0.2e1;
t121 = -t155 / 0.2e1;
t120 = t140 / 0.2e1;
t119 = t139 / 0.2e1;
t66 = pkin(2) * t136 + t89 * t84;
t57 = pkin(8) + t66;
t114 = t126 * t57;
t113 = t126 * t71;
t111 = Ifges(6,5) * t94 - Ifges(6,6) * t91;
t110 = t69 / 0.2e1 + t169;
t97 = (t43 * t112 + t116 * t42) * t164 - t149 / 0.2e1 + t148 / 0.2e1 + mrSges(5,2) * t160 + t142 / 0.2e1 - t141 / 0.2e1 + (-t90 * t42 + t43 * t89) * t163 / 0.2e1 + t43 * t120 + t43 * t119;
t98 = (-t70 * t36 + t71 * t37 - t65 * t38 + t66 * t39) * t165 + (t31 * t113 + t39 * t114 + t70 * t30 + t56 * t38) * t164;
t2 = -t97 - t145 / 0.2e1 - t143 / 0.2e1 + t144 / 0.2e1 + t151 / 0.2e1 - t147 / 0.2e1 - t146 / 0.2e1 - t152 / 0.2e1 - t150 / 0.2e1 + mrSges(4,1) * t122 + mrSges(4,2) * t121 + t98 + (t120 + t119) * t171;
t6 = t127 * t70 + t170 + t109 * t71 + m(5) * (-t65 * t70 + t66 * t71) + m(6) * (t57 * t113 + t56 * t70);
t105 = t2 * qJD(1) + t6 * qJD(2);
t17 = -t45 + t99;
t8 = t162 + t159 + (mrSges(6,2) * t160 - t153) * t94 + (mrSges(6,1) * t160 - t175) * t91;
t104 = -t8 * qJD(1) - t17 * qJD(2);
t103 = -t135 / 0.2e1 - t130 / 0.2e1;
t10 = t162 + t158 + (mrSges(6,2) * t161 - t153) * t94 + (mrSges(6,1) * t161 - t175) * t91;
t13 = t158 + t159 + (mrSges(6,2) * t157 - t153) * t94 + (mrSges(6,1) * t157 - t175) * t91;
t18 = -t69 + t99;
t100 = t10 * qJD(1) + t13 * qJD(2) + t18 * qJD(3);
t44 = t45 / 0.2e1;
t19 = t20 / 0.2e1;
t14 = t103 * t71 + t110 + t44;
t11 = t103 * t39 + t110 + t19;
t9 = t103 * t43 + t169 + t19 + t44;
t1 = t97 + (t157 + t161) * mrSges(5,2) + (-t67 / 0.2e1 + t121) * mrSges(4,2) + (-t68 / 0.2e1 + t122) * mrSges(4,1) + t98 + t127 * (t70 / 0.2e1 + t38 / 0.2e1) + t171 * mrSges(6,3) * (t87 / 0.2e1 + t88 / 0.2e1);
t5 = [qJD(2) * t4 + qJD(3) * t3 - qJD(5) * t15, t1 * qJD(3) + t9 * qJD(5) + t124 + (t148 - t149 + 0.2e1 * (t43 * t114 + t56 * t42) * t164 + 0.2e1 * (-t65 * t42 + t66 * t43) * t165 + m(4) * (t72 * t95 + t73 * t92) * pkin(2) + t167) * qJD(2), t125 + t1 * qJD(2) + (t166 * t39 + t168 * t38 + t108 - t150 + t151 - t152) * qJD(3) + t11 * qJD(5), 0, -t123 + t9 * qJD(2) + t11 * qJD(3) + (t78 * t31 + t111) * qJD(5); qJD(3) * t2 - qJD(5) * t8 - t124, qJD(3) * t6 - qJD(5) * t17, (t166 * t71 + t168 * t70 - t143 + t144 - t145 + t170) * qJD(3) + t14 * qJD(5) + t105, 0, t14 * qJD(3) + (t78 * t57 + t111) * qJD(5) + t104; -qJD(2) * t2 - qJD(5) * t10 - t125, -qJD(5) * t13 - t105, -t18 * qJD(5), 0, (t78 * t83 + t111) * qJD(5) - t100; 0, 0, 0, 0, -t106 * qJD(5); qJD(2) * t8 + qJD(3) * t10 + t123, qJD(3) * t13 - t104, t100, 0, 0;];
Cq = t5;
