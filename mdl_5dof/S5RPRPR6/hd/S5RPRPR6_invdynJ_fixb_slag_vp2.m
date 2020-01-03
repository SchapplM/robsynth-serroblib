% Calculate vector of inverse dynamics joint torques for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:45
% DurationCPUTime: 1.30s
% Computational Cost: add. (1382->192), mult. (2184->246), div. (0->0), fcn. (1060->12), ass. (0->104)
t89 = sin(qJ(5));
t182 = -t89 / 0.2e1;
t92 = cos(qJ(5));
t163 = t92 / 0.2e1;
t181 = mrSges(4,1) - mrSges(5,2);
t114 = mrSges(6,1) * t89 + mrSges(6,2) * t92;
t180 = -mrSges(5,3) - t114;
t88 = cos(pkin(8));
t75 = pkin(1) * t88 + pkin(2);
t117 = t75 * qJD(1);
t87 = sin(pkin(8));
t162 = pkin(1) * t87;
t90 = sin(qJ(3));
t71 = t90 * t162;
t93 = cos(qJ(3));
t25 = qJD(1) * t71 - t93 * t117;
t102 = qJD(4) + t25;
t179 = -mrSges(6,3) - t181;
t178 = mrSges(4,2) + t180;
t115 = mrSges(6,1) * t92 - mrSges(6,2) * t89;
t153 = Ifges(6,4) * t92;
t154 = Ifges(6,4) * t89;
t40 = t93 * t162 + t75 * t90;
t26 = t40 * qJD(1);
t85 = qJD(1) + qJD(3);
t20 = t85 * qJ(4) + t26;
t177 = ((-Ifges(6,1) * t89 - t153) * t163 + (-Ifges(6,2) * t92 - t154) * t182) * t85 + t20 * t115 + qJD(5) * (-Ifges(6,5) * t89 - Ifges(6,6) * t92) / 0.2e1;
t127 = qJD(3) * t162;
t176 = -qJD(1) * t127 + t75 * qJDD(1);
t135 = pkin(1) * qJDD(1);
t126 = t87 * t135;
t175 = qJD(3) * t117 + t126;
t41 = t114 * t85;
t172 = -mrSges(5,3) * t85 - t41;
t10 = t175 * t93 + t176 * t90;
t84 = qJDD(1) + qJDD(3);
t170 = t84 * qJ(4) + t85 * qJD(4);
t5 = -t10 - t170;
t171 = -t5 * qJ(4) + t102 * t20;
t133 = qJD(5) * t89;
t42 = -t85 * t133 + t84 * t92;
t27 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t42;
t132 = qJD(5) * t92;
t43 = -t85 * t132 - t84 * t89;
t28 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t43;
t169 = t92 * t27 + t89 * t28;
t86 = qJ(1) + pkin(8);
t82 = qJ(3) + t86;
t73 = sin(t82);
t74 = cos(t82);
t168 = -g(1) * t73 + g(2) * t74;
t95 = -pkin(3) - pkin(7);
t167 = t178 * t74 + (-m(6) * t95 - t179) * t73;
t166 = t178 * t73 + t179 * t74;
t145 = t89 * mrSges(6,3);
t50 = -qJD(5) * mrSges(6,2) - t85 * t145;
t142 = t92 * mrSges(6,3);
t51 = qJD(5) * mrSges(6,1) - t85 * t142;
t18 = t95 * t85 + t102;
t14 = -qJD(2) * t89 + t18 * t92;
t11 = -t175 * t90 + t176 * t93;
t101 = qJDD(4) - t11;
t4 = t95 * t84 + t101;
t1 = qJD(5) * t14 + qJDD(2) * t92 + t4 * t89;
t15 = qJD(2) * t92 + t18 * t89;
t109 = -t14 * t89 + t15 * t92;
t2 = -qJD(5) * t15 - qJDD(2) * t89 + t4 * t92;
t98 = t109 * qJD(5) + t1 * t89 + t2 * t92;
t165 = t98 * m(6) + t50 * t132 - t51 * t133 + t169;
t110 = t14 * t92 + t15 * t89;
t164 = -m(5) * t102 - m(6) * t110 - t89 * t50 - t92 * t51 + (m(5) * pkin(3) + t181) * t85;
t91 = sin(qJ(1));
t161 = pkin(1) * t91;
t160 = pkin(3) * t73;
t94 = cos(qJ(1));
t83 = t94 * pkin(1);
t150 = t75 * t93;
t149 = t84 * mrSges(5,2);
t147 = t85 * mrSges(4,2);
t138 = t74 * pkin(3) + t73 * qJ(4);
t81 = cos(t86);
t137 = pkin(2) * t81 + t83;
t130 = m(3) + m(4) + m(5);
t39 = -t71 + t150;
t121 = t137 + t138;
t38 = -pkin(3) - t39;
t80 = sin(t86);
t119 = -pkin(2) * t80 - t161;
t32 = qJD(3) * t150 - t90 * t127;
t29 = -qJD(4) - t32;
t37 = qJ(4) + t40;
t118 = -t20 * t29 - t37 * t5;
t113 = Ifges(6,1) * t92 - t154;
t112 = -Ifges(6,2) * t89 + t153;
t63 = t74 * qJ(4);
t107 = t119 + t63;
t100 = -t20 * t85 + t168;
t34 = Ifges(6,6) * qJD(5) + t112 * t85;
t35 = Ifges(6,5) * qJD(5) + t113 * t85;
t7 = -pkin(3) * t84 + t101;
t97 = t11 * mrSges(4,1) - t10 * mrSges(4,2) + t7 * mrSges(5,2) - t1 * t145 - t2 * t142 + (Ifges(6,4) * t42 + Ifges(6,2) * t43) * t182 + (Ifges(6,1) * t42 + Ifges(6,4) * t43) * t163 + t42 * t113 / 0.2e1 + t43 * t112 / 0.2e1 - t35 * t133 / 0.2e1 - t34 * t132 / 0.2e1 + (Ifges(4,3) + Ifges(5,1)) * t84 + t180 * t5 + (-t15 * t132 + t14 * t133) * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t163 - Ifges(6,6) * t89) * qJDD(5) + t177 * qJD(5);
t69 = t74 * pkin(7);
t17 = -mrSges(6,1) * t43 + mrSges(6,2) * t42;
t3 = [m(6) * t118 + m(5) * (t38 * t7 + t118) + 0.2e1 * t88 * mrSges(3,1) * t135 + t97 + t38 * t149 - t32 * t147 - 0.2e1 * mrSges(3,2) * t126 + m(4) * (t10 * t40 + t11 * t39 + t26 * t32) + t37 * t17 + t172 * t29 + (t39 * mrSges(4,1) - t40 * mrSges(4,2) + t37 * mrSges(5,3)) * t84 + (Ifges(3,3) + Ifges(2,3) + m(3) * (t87 ^ 2 + t88 ^ 2) * pkin(1) ^ 2) * qJDD(1) + t165 * (-pkin(7) + t38) + (m(4) * t25 - t164) * t40 * qJD(3) + (-m(3) * t83 - mrSges(3,1) * t81 + mrSges(3,2) * t80 - m(4) * t137 - m(6) * (t69 + t121) - m(5) * t121 - mrSges(2,1) * t94 + t91 * mrSges(2,2) + t166) * g(2) + (-m(4) * t119 - m(6) * t107 - m(5) * (t107 - t160) + m(3) * t161 + mrSges(3,1) * t80 + mrSges(3,2) * t81 + mrSges(2,1) * t91 + mrSges(2,2) * t94 + t167) * g(1); -t50 * t133 - t51 * t132 + t92 * t28 - t89 * t27 + m(6) * (-t110 * qJD(5) + t1 * t92 - t2 * t89) + t130 * qJDD(2) + (-m(6) - t130) * g(3); -pkin(3) * t149 + qJ(4) * t17 + qJD(4) * t41 + t97 + t170 * mrSges(5,3) + t171 * m(6) + (-t7 * pkin(3) + t171) * m(5) + (-t147 - t172) * t25 + (-m(6) * (t69 + t138) - m(5) * t138 + t166) * g(2) + (-m(6) * t63 - m(5) * (t63 - t160) + t167) * g(1) + t165 * t95 + t164 * t26; t149 + t172 * t85 + (t92 * t50 - t89 * t51) * qJD(5) + (t100 + t98) * m(6) + (t100 + t7) * m(5) + t169; t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t42 + Ifges(6,6) * t43 + Ifges(6,3) * qJDD(5) + g(3) * t114 - t14 * t50 + t15 * t51 + (t89 * t35 / 0.2e1 + t34 * t163 + t109 * mrSges(6,3) - t177) * t85 + t168 * t115;];
tau = t3;
