% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:56
% EndTime: 2019-12-31 19:29:01
% DurationCPUTime: 2.00s
% Computational Cost: add. (4109->179), mult. (7804->231), div. (0->0), fcn. (8238->6), ass. (0->104)
t115 = sin(qJ(5));
t116 = cos(qJ(5));
t150 = sin(pkin(8));
t151 = cos(pkin(8));
t171 = sin(qJ(2));
t172 = cos(qJ(2));
t96 = t150 * t171 - t151 * t172;
t98 = -t150 * t172 - t151 * t171;
t128 = t115 * t96 - t116 * t98;
t141 = t171 * pkin(6);
t104 = -t171 * qJ(3) - t141;
t143 = t172 * pkin(6);
t105 = t172 * qJ(3) + t143;
t185 = t150 * t104 + t151 * t105;
t186 = t96 * pkin(7) + t185;
t79 = t151 * t104 - t150 * t105;
t48 = -t98 * pkin(7) + t79;
t210 = -t115 * t48 + t116 * t186;
t29 = t115 * t186 + t116 * t48;
t71 = t115 * t98 + t116 * t96;
t221 = -t210 * mrSges(6,1) + t29 * mrSges(6,2) + Ifges(6,5) * t71 - Ifges(6,6) * t128;
t223 = t221 * qJD(5);
t113 = -t172 * pkin(2) - pkin(1);
t120 = t98 * qJ(4) + t113;
t42 = (-pkin(3) - pkin(4)) * t96 - t120;
t222 = -m(6) * t42 + mrSges(6,1) * t71 - mrSges(6,2) * t128;
t217 = t115 * t29 + t116 * t210;
t66 = t96 * pkin(3) + t120;
t213 = m(5) * t66 + mrSges(5,1) * t96 + mrSges(5,3) * t98;
t200 = t71 * mrSges(6,2);
t135 = -t128 * mrSges(6,1) - t200;
t170 = Ifges(6,4) * t128;
t193 = t128 / 0.2e1;
t212 = (Ifges(6,1) * t71 - t170) * t193 - t42 * t135;
t211 = -t98 * mrSges(4,1) - t96 * mrSges(4,2);
t140 = -t200 / 0.2e1;
t68 = Ifges(6,4) * t71;
t208 = Ifges(6,2) * t128 - t68;
t205 = -t71 / 0.2e1;
t204 = t71 / 0.2e1;
t142 = t171 * pkin(2);
t202 = m(4) * t142;
t129 = t115 * mrSges(6,1) + t116 * mrSges(6,2);
t197 = t129 * qJD(5);
t136 = t150 * pkin(2);
t107 = t136 + qJ(4);
t196 = m(5) * t107 + mrSges(5,3);
t194 = -t128 / 0.2e1;
t191 = -Ifges(5,5) + Ifges(4,4);
t190 = mrSges(6,3) * t128;
t184 = 0.2e1 * t98;
t183 = m(5) / 0.2e1;
t182 = m(6) / 0.2e1;
t181 = m(6) / 0.4e1;
t180 = m(4) * pkin(2);
t158 = t96 * mrSges(5,2);
t5 = (-t128 ^ 2 - t71 ^ 2) * mrSges(6,3) + m(6) * (-t128 * t29 - t210 * t71) + (t96 ^ 2 + t98 ^ 2) * (mrSges(4,3) + mrSges(5,2)) + (m(4) + m(5)) * (-t185 * t96 + t79 * t98);
t154 = qJD(1) * t5;
t34 = Ifges(6,2) * t71 + t170;
t35 = Ifges(6,1) * t128 + t68;
t69 = -t98 * pkin(3) + t96 * qJ(4) + t142;
t43 = -t98 * pkin(4) + t69;
t91 = t96 * mrSges(5,3);
t1 = -pkin(1) * (t171 * mrSges(3,1) + t172 * mrSges(3,2)) + t34 * t193 + t35 * t205 + t208 * t204 + t66 * t91 + (-t66 * mrSges(5,1) - mrSges(4,2) * t142 - t191 * t98) * t98 + (mrSges(4,1) * t142 + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t98 + t191 * t96) * t96 + (-Ifges(3,2) + Ifges(3,1)) * t172 * t171 + (-t171 ^ 2 + t172 ^ 2) * Ifges(3,4) + t213 * t69 + (t202 + t211) * t113 - t212 + t222 * t43;
t153 = t1 * qJD(1);
t137 = t151 * pkin(2);
t112 = -t137 - pkin(3);
t106 = -pkin(4) + t112;
t83 = t106 * t116 - t115 * t107;
t84 = t115 * t106 + t107 * t116;
t117 = (-t107 * t96 - t112 * t98) * t183 + (t128 * t83 - t71 * t84) * t182 + (-t150 * t96 + t151 * t98) * t180 / 0.2e1;
t118 = t69 * t183 + t43 * t182 + t202 / 0.2e1;
t9 = -t98 * mrSges(5,1) - t117 + t118 - t135 + t211 + t91;
t152 = t9 * qJD(1);
t13 = (t213 + t222) * t98;
t149 = qJD(1) * t13;
t19 = t84 * mrSges(6,1) + mrSges(6,2) * t83;
t148 = qJD(5) * t19;
t17 = 0.2e1 * t194 * mrSges(6,1) + 0.2e1 * t140;
t147 = t17 * qJD(1);
t122 = m(6) * (-t115 * t71 + t116 * t128);
t21 = -t122 / 0.2e1 + (t181 + t183) * t184;
t146 = t21 * qJD(1);
t145 = t115 * t190;
t139 = t205 + t204;
t138 = t193 + t194;
t4 = t139 * Ifges(6,5) + t138 * Ifges(6,6);
t127 = t4 * qJD(1) + t19 * qJD(2);
t32 = m(6) * (-t115 * t83 + t84 * t116) + t129 + t196;
t119 = t217 * t182;
t123 = m(6) * t217;
t6 = (-t138 * t115 + t139 * t116) * mrSges(6,3) + t119 - t123 / 0.2e1;
t126 = -t6 * qJD(1) + t32 * qJD(2);
t54 = t145 / 0.2e1;
t12 = -t145 / 0.2e1 + t54;
t2 = t34 * t194 + (t35 - t208) * t204 + t212;
t125 = t2 * qJD(1) + t12 * qJD(4);
t124 = -t12 * qJD(1) - qJD(2) * t129;
t20 = t122 / 0.2e1 - m(5) * t98 / 0.2e1 + (t181 + m(5) / 0.4e1) * t184;
t18 = t140 + t200 / 0.2e1;
t11 = t12 * qJD(5);
t10 = t117 + t118;
t7 = t54 - t158 + t123 / 0.2e1 + t119 + 0.2e1 * t185 * t183 + (t115 * t193 + 0.2e1 * t204 * t116) * mrSges(6,3);
t3 = [qJD(2) * t1 + qJD(3) * t5 + qJD(4) * t13 + qJD(5) * t2, t10 * qJD(3) + t7 * qJD(4) - t223 + t153 + (t84 * t190 + t83 * t71 * mrSges(6,3) + mrSges(3,2) * t141 + m(6) * (t210 * t83 + t29 * t84) + Ifges(3,5) * t172 - Ifges(3,6) * t171 - t112 * t158 - mrSges(3,1) * t143 + (mrSges(4,3) * t137 - Ifges(5,4) - Ifges(4,5)) * t96 + (t107 * mrSges(5,2) + mrSges(4,3) * t136 + Ifges(4,6) - Ifges(5,6)) * t98 + (t150 * t180 - mrSges(4,2) + t196) * t79 + (m(5) * t112 - t151 * t180 - mrSges(4,1) - mrSges(5,1)) * t185 + t221) * qJD(2), qJD(2) * t10 + qJD(4) * t20 + qJD(5) * t18 + t154, qJD(2) * t7 + qJD(3) * t20 + t11 + t149, -qJD(2) * t221 + t18 * qJD(3) + t125 + t223; -qJD(3) * t9 - qJD(4) * t6 + qJD(5) * t4 - t153, qJD(4) * t32 + t148, -t152, t126, t127 - t148; qJD(2) * t9 + qJD(4) * t21 + qJD(5) * t17 - t154, t152, 0, t146, t147; qJD(2) * t6 - qJD(3) * t21 + t11 - t149, -t126 + t197, -t146, 0, -t124 - t197; -qJD(2) * t4 - qJD(3) * t17 - t125, -qJD(4) * t129 - t127, -t147, t124, 0;];
Cq = t3;
