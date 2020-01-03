% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:51
% EndTime: 2020-01-03 12:01:54
% DurationCPUTime: 1.40s
% Computational Cost: add. (5318->137), mult. (10116->191), div. (0->0), fcn. (9364->8), ass. (0->101)
t125 = cos(qJ(4));
t239 = t125 ^ 2;
t123 = sin(qJ(4));
t201 = sin(qJ(5));
t202 = cos(qJ(5));
t104 = -t201 * t123 + t202 * t125;
t105 = -t202 * t123 - t201 * t125;
t121 = sin(pkin(9));
t124 = sin(qJ(2));
t111 = t121 * t124 * pkin(1);
t122 = cos(pkin(9));
t126 = cos(qJ(2));
t200 = pkin(1) * t126;
t91 = t122 * t200 - t111;
t60 = t105 * t91;
t61 = t104 * t91;
t24 = t104 * t60 - t105 * t61;
t238 = m(6) * t24;
t109 = t123 * mrSges(5,1) + t125 * mrSges(5,2);
t207 = m(6) * pkin(4);
t219 = t60 * mrSges(6,1) / 0.2e1 - t61 * mrSges(6,2) / 0.2e1;
t237 = (t201 * t61 + t202 * t60) * t207 / 0.2e1 - t109 * t91 / 0.2e1 + t219;
t112 = pkin(2) * t121 + pkin(7);
t183 = pkin(8) + t112;
t93 = t183 * t123;
t94 = t183 * t125;
t146 = -t201 * t94 - t202 * t93;
t67 = t201 * t93 - t202 * t94;
t236 = t67 * mrSges(6,1) - t146 * mrSges(6,2);
t116 = pkin(2) + t200;
t165 = t122 * t124;
t149 = pkin(1) * t165 + t121 * t116;
t85 = pkin(7) + t149;
t203 = pkin(8) + t85;
t77 = t203 * t123;
t78 = t203 * t125;
t147 = -t201 * t78 - t202 * t77;
t39 = t201 * t77 - t202 * t78;
t235 = t39 * mrSges(6,1) - t147 * mrSges(6,2);
t181 = Ifges(6,5) * t104 + Ifges(6,6) * t105;
t11 = t181 + t235;
t234 = t11 * qJD(5);
t233 = t239 * Ifges(5,4);
t15 = t181 + t236;
t232 = t15 * qJD(5);
t231 = -Ifges(5,4) * t123 + (Ifges(5,1) - Ifges(5,2)) * t125;
t225 = -Ifges(6,1) + Ifges(6,2);
t74 = -t105 * mrSges(6,1) + t104 * mrSges(6,2);
t224 = t74 * qJD(5);
t177 = Ifges(6,4) * t105;
t98 = Ifges(6,4) * t104;
t215 = (t225 * t105 + t98) * t104;
t75 = -mrSges(6,1) * t104 - mrSges(6,2) * t105;
t221 = -t177 * t105 + t233 + (pkin(4) * t75 + t231) * t123 + t215;
t214 = -mrSges(5,1) * t125 + mrSges(5,2) * t123;
t164 = t123 ^ 2 + t239;
t90 = (t121 * t126 + t165) * pkin(1);
t210 = (t164 * mrSges(5,3) - mrSges(4,2)) * t91 + (-mrSges(4,1) + t75 + t214) * t90 + (-t124 * mrSges(3,1) - t126 * mrSges(3,2)) * pkin(1) + (t104 * t61 + t105 * t60) * mrSges(6,3);
t209 = 0.2e1 * m(6);
t208 = m(6) / 0.2e1;
t199 = pkin(4) * t123;
t198 = pkin(4) * t125;
t151 = t116 * t122 - t111;
t84 = -pkin(3) - t151;
t79 = t84 - t198;
t187 = t79 * t74;
t113 = -pkin(2) * t122 - pkin(3);
t107 = t113 - t198;
t172 = t107 * t74;
t167 = t84 * t109;
t166 = t113 * t109;
t163 = qJD(1) * t238;
t162 = pkin(4) * t202;
t161 = pkin(4) * t201;
t152 = t164 * t91;
t148 = (t79 / 0.2e1 + t107 / 0.2e1) * t74;
t76 = t79 * t199;
t5 = m(6) * t76 + t167 + t187 + t221;
t140 = t5 * qJD(1);
t7 = m(6) * (t147 * t60 - t39 * t61 + t79 * t90) + m(5) * (t85 * t152 + t84 * t90) + m(4) * (t149 * t91 - t151 * t90) + t210;
t139 = -qJD(3) * t238 / 0.2e1 - t7 * qJD(1);
t133 = -Ifges(6,4) * t105 ^ 2 + t215;
t9 = t133 + t187;
t138 = t9 * qJD(1);
t92 = t107 * t199;
t129 = (t76 + t92) * t208;
t1 = t129 + t75 * t199 + t187 / 0.2e1 + t104 * t98 + t172 / 0.2e1 + t167 / 0.2e1 + t166 / 0.2e1 + t233 + t231 * t123 + (t104 * t225 - t177) * t105 - t237;
t8 = m(6) * t92 + t166 + t172 + t221;
t136 = -t1 * qJD(1) - t8 * qJD(2);
t10 = t133 + t172;
t128 = t148 + t133;
t3 = t128 - t219;
t135 = -t3 * qJD(1) - t10 * qJD(2);
t132 = Ifges(5,5) * t125 - Ifges(5,6) * t123 + t181 + (-t104 * t162 + t105 * t161) * mrSges(6,3);
t106 = (mrSges(6,1) * t201 + mrSges(6,2) * t202) * pkin(4);
t131 = t106 * qJD(4);
t103 = t106 * qJD(5);
t6 = t24 * qJD(2) * t209 / 0.4e1;
t4 = t128 + t219;
t2 = t129 + t148 + (t84 / 0.2e1 + t113 / 0.2e1) * t109 + t221 + t237;
t12 = [qJD(2) * t7 + qJD(4) * t5 + qJD(5) * t9, t2 * qJD(4) + t4 * qJD(5) - t139 + (0.2e1 * (t107 * t90 + t146 * t60 - t61 * t67) * t208 + m(5) * (t112 * t152 + t113 * t90) + m(4) * (t121 * t91 - t122 * t90) * pkin(2) + t210) * qJD(2), t6, t2 * qJD(2) + ((t147 * t201 + t202 * t39) * t207 + t132 + t214 * t85 + t235) * qJD(4) + t234 + t140, t4 * qJD(2) + t11 * qJD(4) + t138 + t234; t1 * qJD(4) + t3 * qJD(5) + t139, qJD(4) * t8 + qJD(5) * t10, -t163 / 0.2e1, ((t146 * t201 + t202 * t67) * t207 + t132 + t214 * t112 + t236) * qJD(4) + t232 - t136, t15 * qJD(4) - t135 + t232; t6, t163 / 0.2e1, 0, -t224 + (-t109 - t74 + (t104 * t161 + t105 * t162) * t209 / 0.2e1) * qJD(4), -qJD(4) * t74 - t224; -t1 * qJD(2) - t140, t136, 0, -t103, -t103 - t131; -t3 * qJD(2) - t138, t135, 0, t131, 0;];
Cq = t12;
