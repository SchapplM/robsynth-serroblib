% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 04:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:41:15
% EndTime: 2019-05-06 04:41:50
% DurationCPUTime: 19.30s
% Computational Cost: add. (313828->221), mult. (973532->319), div. (0->0), fcn. (841013->16), ass. (0->132)
t120 = cos(qJ(3));
t107 = sin(pkin(7));
t109 = cos(pkin(13));
t111 = cos(pkin(6));
t108 = sin(pkin(6));
t110 = cos(pkin(7));
t153 = t108 * t110;
t131 = t107 * t111 + t109 * t153;
t129 = t131 * t120;
t115 = sin(qJ(3));
t106 = sin(pkin(13));
t156 = t106 * t108;
t147 = t115 * t156;
t152 = t110 * t115;
t155 = t107 * t115;
t125 = t111 * t155 + (t106 * t120 + t109 * t152) * t108;
t86 = t125 * qJD(1);
t74 = -t86 * qJD(3) + (t129 - t147) * qJDD(1);
t163 = pkin(9) * t106;
t135 = -pkin(2) * t109 - t107 * t163;
t151 = qJD(1) * t108;
t158 = pkin(9) * qJDD(1);
t132 = qJD(1) * t135 * t151 + t110 * t158;
t146 = qJD(2) * t151;
t154 = t108 * t109;
t122 = qJD(1) ^ 2;
t116 = sin(qJ(1));
t121 = cos(qJ(1));
t145 = t116 * g(1) - t121 * g(2);
t157 = qJ(2) * t108;
t96 = qJDD(1) * pkin(1) + t122 * t157 + t145;
t160 = t111 * t96;
t136 = -g(3) * t154 - 0.2e1 * t106 * t146 + t109 * t160;
t128 = t131 * qJD(1);
t93 = pkin(9) * t128;
t141 = -t121 * g(1) - t116 * g(2);
t97 = -t122 * pkin(1) + qJDD(1) * t157 + t141;
t57 = (pkin(2) * qJDD(1) + qJD(1) * t93) * t111 + (-t132 * t108 - t97) * t106 + t136;
t148 = t106 * t160 + (0.2e1 * t146 + t97) * t109;
t98 = (pkin(2) * t111 - t153 * t163) * qJD(1);
t58 = (-qJD(1) * t98 + t107 * t158) * t111 + (-g(3) * t106 + t132 * t109) * t108 + t148;
t142 = -t111 * g(3) + qJDD(2);
t66 = (-t96 + t135 * qJDD(1) + (t106 * t98 - t109 * t93) * qJD(1)) * t108 + t142;
t164 = -t115 * t58 + (t107 * t66 + t110 * t57) * t120;
t113 = sin(qJ(5));
t118 = cos(qJ(5));
t114 = sin(qJ(4));
t119 = cos(qJ(4));
t149 = t120 * t58 + t57 * t152 + t66 * t155;
t102 = qJD(1) * t147;
t85 = t120 * t128 - t102;
t73 = -t85 * pkin(3) - t86 * pkin(10);
t130 = -t107 * t154 + t110 * t111;
t94 = t130 * qJD(1) + qJD(3);
t90 = t94 ^ 2;
t91 = t130 * qJDD(1) + qJDD(3);
t33 = -t90 * pkin(3) + t91 * pkin(10) + t85 * t73 + t149;
t144 = -t107 * t57 + t110 * t66;
t75 = t85 * qJD(3) + t125 * qJDD(1);
t36 = (-t85 * t94 - t75) * pkin(10) + (t86 * t94 - t74) * pkin(3) + t144;
t143 = -t114 * t33 + t119 * t36;
t77 = -t114 * t86 + t119 * t94;
t53 = t77 * qJD(4) + t114 * t91 + t119 * t75;
t71 = qJDD(4) - t74;
t78 = t114 * t94 + t119 * t86;
t84 = -qJD(1) * t129 + qJD(4) + t102;
t24 = (t77 * t84 - t53) * pkin(11) + (t77 * t78 + t71) * pkin(4) + t143;
t161 = t114 * t36 + t119 * t33;
t52 = -t78 * qJD(4) - t114 * t75 + t119 * t91;
t69 = t84 * pkin(4) - t78 * pkin(11);
t76 = t77 ^ 2;
t26 = -t76 * pkin(4) + t52 * pkin(11) - t84 * t69 + t161;
t162 = t113 * t24 + t118 * t26;
t112 = sin(qJ(6));
t117 = cos(qJ(6));
t60 = -t113 * t78 + t118 * t77;
t61 = t113 * t77 + t118 * t78;
t46 = -t60 * pkin(5) - t61 * pkin(12);
t70 = qJDD(5) + t71;
t82 = qJD(5) + t84;
t81 = t82 ^ 2;
t21 = -t81 * pkin(5) + t70 * pkin(12) + t60 * t46 + t162;
t32 = -t91 * pkin(3) - t90 * pkin(10) + t86 * t73 - t164;
t124 = -t52 * pkin(4) - t76 * pkin(11) + t78 * t69 + t32;
t39 = -t61 * qJD(5) - t113 * t53 + t118 * t52;
t40 = t60 * qJD(5) + t113 * t52 + t118 * t53;
t22 = (-t60 * t82 - t40) * pkin(12) + (t61 * t82 - t39) * pkin(5) + t124;
t47 = -t112 * t61 + t117 * t82;
t30 = t47 * qJD(6) + t112 * t70 + t117 * t40;
t38 = qJDD(6) - t39;
t48 = t112 * t82 + t117 * t61;
t41 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t59 = qJD(6) - t60;
t42 = -t59 * mrSges(7,2) + t47 * mrSges(7,3);
t18 = m(7) * (-t112 * t21 + t117 * t22) - t30 * mrSges(7,3) + t38 * mrSges(7,1) - t48 * t41 + t59 * t42;
t29 = -t48 * qJD(6) - t112 * t40 + t117 * t70;
t43 = t59 * mrSges(7,1) - t48 * mrSges(7,3);
t19 = m(7) * (t112 * t22 + t117 * t21) + t29 * mrSges(7,3) - t38 * mrSges(7,2) + t47 * t41 - t59 * t43;
t45 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t50 = t82 * mrSges(6,1) - t61 * mrSges(6,3);
t14 = m(6) * t162 - t70 * mrSges(6,2) + t39 * mrSges(6,3) - t112 * t18 + t117 * t19 + t60 * t45 - t82 * t50;
t137 = -t113 * t26 + t118 * t24;
t126 = m(7) * (-t70 * pkin(5) - t81 * pkin(12) + t61 * t46 - t137) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t47 * t42 + t48 * t43;
t49 = -t82 * mrSges(6,2) + t60 * mrSges(6,3);
t15 = m(6) * t137 + t70 * mrSges(6,1) - t40 * mrSges(6,3) - t61 * t45 + t82 * t49 - t126;
t62 = -t77 * mrSges(5,1) + t78 * mrSges(5,2);
t67 = -t84 * mrSges(5,2) + t77 * mrSges(5,3);
t11 = m(5) * t143 + t71 * mrSges(5,1) - t53 * mrSges(5,3) + t113 * t14 + t118 * t15 - t78 * t62 + t84 * t67;
t68 = t84 * mrSges(5,1) - t78 * mrSges(5,3);
t12 = m(5) * t161 - t71 * mrSges(5,2) + t52 * mrSges(5,3) - t113 * t15 + t118 * t14 + t77 * t62 - t84 * t68;
t79 = -t94 * mrSges(4,2) + t85 * mrSges(4,3);
t80 = t94 * mrSges(4,1) - t86 * mrSges(4,3);
t10 = m(4) * t144 - t74 * mrSges(4,1) + t75 * mrSges(4,2) + t119 * t11 + t114 * t12 - t85 * t79 + t86 * t80;
t133 = -t111 * mrSges(3,2) + mrSges(3,3) * t154;
t100 = t133 * qJD(1);
t134 = t111 * mrSges(3,1) - mrSges(3,3) * t156;
t127 = -m(6) * t124 + t39 * mrSges(6,1) - t40 * mrSges(6,2) - t112 * t19 - t117 * t18 + t60 * t49 - t61 * t50;
t123 = m(5) * t32 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t77 * t67 + t78 * t68 - t127;
t72 = -t85 * mrSges(4,1) + t86 * mrSges(4,2);
t13 = m(4) * t164 + t91 * mrSges(4,1) - t75 * mrSges(4,3) - t86 * t72 + t94 * t79 - t123;
t9 = m(4) * t149 - t91 * mrSges(4,2) + t74 * mrSges(4,3) - t114 * t11 + t119 * t12 + t85 * t72 - t94 * t80;
t140 = t115 * t9 + t120 * t13;
t139 = -mrSges(3,1) * t109 + mrSges(3,2) * t106;
t95 = t139 * t151;
t4 = m(3) * (-t106 * t97 + t136) - t107 * t10 + t140 * t110 + t134 * qJDD(1) + (t111 * t100 - t95 * t156) * qJD(1);
t99 = t134 * qJD(1);
t6 = m(3) * t142 + t110 * t10 + t140 * t107 + (-m(3) * t96 + t139 * qJDD(1) + (-t100 * t109 + t106 * t99) * qJD(1)) * t108;
t8 = m(3) * (-g(3) * t156 + t148) + t120 * t9 - t115 * t13 + t133 * qJDD(1) + (-t111 * t99 + t95 * t154) * qJD(1);
t150 = t111 * t6 + t4 * t154 + t8 * t156;
t2 = m(2) * t141 - t122 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t109 * t8;
t1 = m(2) * t145 + qJDD(1) * mrSges(2,1) - t122 * mrSges(2,2) - t108 * t6 + (t106 * t8 + t109 * t4) * t111;
t3 = [-m(1) * g(1) - t116 * t1 + t121 * t2, t2, t8, t9, t12, t14, t19; -m(1) * g(2) + t121 * t1 + t116 * t2, t1, t4, t13, t11, t15, t18; (-m(1) - m(2)) * g(3) + t150, -m(2) * g(3) + t150, t6, t10, t123, -t127, t126;];
f_new  = t3;
