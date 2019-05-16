% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 16:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 15:41:23
% EndTime: 2019-05-08 15:42:38
% DurationCPUTime: 24.20s
% Computational Cost: add. (446156->229), mult. (1097703->325), div. (0->0), fcn. (933942->16), ass. (0->131)
t107 = sin(pkin(7));
t114 = sin(qJ(3));
t120 = cos(qJ(3));
t109 = cos(pkin(7));
t108 = sin(pkin(6));
t115 = sin(qJ(2));
t121 = cos(qJ(2));
t143 = qJD(1) * qJD(2);
t102 = (qJDD(1) * t115 + t121 * t143) * t108;
t110 = cos(pkin(6));
t104 = t110 * qJDD(1) + qJDD(2);
t105 = t110 * qJD(1) + qJD(2);
t123 = qJD(1) ^ 2;
t116 = sin(qJ(1));
t122 = cos(qJ(1));
t132 = -t122 * g(1) - t116 * g(2);
t161 = pkin(9) * t108;
t100 = -t123 * pkin(1) + qJDD(1) * t161 + t132;
t137 = t116 * g(1) - t122 * g(2);
t99 = qJDD(1) * pkin(1) + t123 * t161 + t137;
t153 = t110 * t99;
t133 = -t115 * t100 + t121 * t153;
t146 = qJD(1) * t115;
t159 = pkin(10) * t109;
t145 = qJD(1) * t121;
t139 = t108 * t145;
t152 = t105 * t107;
t92 = (t109 * t139 + t152) * pkin(10);
t147 = qJD(1) * t108;
t160 = pkin(10) * t107;
t96 = (-pkin(2) * t121 - t115 * t160) * t147;
t57 = -t102 * t159 + t104 * pkin(2) + t105 * t92 + (-g(3) * t121 - t96 * t146) * t108 + t133;
t154 = t109 * t57;
t103 = (qJDD(1) * t121 - t115 * t143) * t108;
t128 = t103 * t109 + t104 * t107;
t155 = t121 * t100 + t115 * t153;
t140 = t108 * t146;
t95 = t105 * pkin(2) - t140 * t159;
t58 = -t105 * t95 + (-g(3) * t115 + t96 * t145) * t108 + t128 * pkin(10) + t155;
t158 = t110 * g(3);
t63 = -t102 * t160 - t103 * pkin(2) - t158 + (-t99 + (t115 * t95 - t121 * t92) * qJD(1)) * t108;
t162 = -t114 * t58 + (t107 * t63 + t154) * t120;
t148 = t109 * t121;
t151 = t107 * t114;
t87 = t105 * t151 + (t114 * t148 + t115 * t120) * t147;
t72 = -t87 * qJD(3) - t114 * t102 + t128 * t120;
t86 = (-t114 * t115 + t120 * t148) * t147 + t120 * t152;
t112 = sin(qJ(5));
t118 = cos(qJ(5));
t113 = sin(qJ(4));
t119 = cos(qJ(4));
t142 = t114 * t154 + t120 * t58 + t63 * t151;
t75 = -t86 * pkin(3) - t87 * pkin(11);
t88 = -t107 * t103 + t109 * t104 + qJDD(3);
t93 = t109 * t105 - t107 * t139 + qJD(3);
t91 = t93 ^ 2;
t34 = -t91 * pkin(3) + t88 * pkin(11) + t86 * t75 + t142;
t136 = -t107 * t57 + t109 * t63;
t73 = t86 * qJD(3) + t120 * t102 + t128 * t114;
t36 = (-t86 * t93 - t73) * pkin(11) + (t87 * t93 - t72) * pkin(3) + t136;
t156 = t113 * t36 + t119 * t34;
t79 = -t113 * t87 + t119 * t93;
t80 = t113 * t93 + t119 * t87;
t60 = -t79 * pkin(4) - t80 * pkin(12);
t71 = qJDD(4) - t72;
t85 = qJD(4) - t86;
t84 = t85 ^ 2;
t24 = -t84 * pkin(4) + t71 * pkin(12) + t79 * t60 + t156;
t33 = -t88 * pkin(3) - t91 * pkin(11) + t87 * t75 - t162;
t48 = -t80 * qJD(4) - t113 * t73 + t119 * t88;
t49 = t79 * qJD(4) + t113 * t88 + t119 * t73;
t27 = (-t79 * t85 - t49) * pkin(12) + (t80 * t85 - t48) * pkin(4) + t33;
t157 = t112 * t27 + t118 * t24;
t150 = t108 * t115;
t149 = t108 * t121;
t111 = sin(qJ(6));
t117 = cos(qJ(6));
t135 = -t112 * t24 + t118 * t27;
t66 = -t112 * t80 + t118 * t85;
t40 = t66 * qJD(5) + t112 * t71 + t118 * t49;
t47 = qJDD(5) - t48;
t67 = t112 * t85 + t118 * t80;
t78 = qJD(5) - t79;
t18 = (t66 * t78 - t40) * pkin(13) + (t66 * t67 + t47) * pkin(5) + t135;
t39 = -t67 * qJD(5) - t112 * t49 + t118 * t71;
t53 = t78 * pkin(5) - t67 * pkin(13);
t65 = t66 ^ 2;
t19 = -t65 * pkin(5) + t39 * pkin(13) - t78 * t53 + t157;
t43 = -t111 * t67 + t117 * t66;
t30 = t43 * qJD(6) + t111 * t39 + t117 * t40;
t44 = t111 * t66 + t117 * t67;
t38 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t76 = qJD(6) + t78;
t41 = -t76 * mrSges(7,2) + t43 * mrSges(7,3);
t46 = qJDD(6) + t47;
t16 = m(7) * (-t111 * t19 + t117 * t18) - t30 * mrSges(7,3) + t46 * mrSges(7,1) - t44 * t38 + t76 * t41;
t29 = -t44 * qJD(6) - t111 * t40 + t117 * t39;
t42 = t76 * mrSges(7,1) - t44 * mrSges(7,3);
t17 = m(7) * (t111 * t18 + t117 * t19) + t29 * mrSges(7,3) - t46 * mrSges(7,2) + t43 * t38 - t76 * t42;
t45 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t51 = -t78 * mrSges(6,2) + t66 * mrSges(6,3);
t13 = m(6) * t135 + t47 * mrSges(6,1) - t40 * mrSges(6,3) + t111 * t17 + t117 * t16 - t67 * t45 + t78 * t51;
t52 = t78 * mrSges(6,1) - t67 * mrSges(6,3);
t14 = m(6) * t157 - t47 * mrSges(6,2) + t39 * mrSges(6,3) - t111 * t16 + t117 * t17 + t66 * t45 - t78 * t52;
t59 = -t79 * mrSges(5,1) + t80 * mrSges(5,2);
t69 = t85 * mrSges(5,1) - t80 * mrSges(5,3);
t12 = m(5) * t156 - t71 * mrSges(5,2) + t48 * mrSges(5,3) - t112 * t13 + t118 * t14 + t79 * t59 - t85 * t69;
t134 = -t113 * t34 + t119 * t36;
t23 = -t71 * pkin(4) - t84 * pkin(12) + t80 * t60 - t134;
t126 = t29 * mrSges(7,1) + t43 * t41 - m(7) * (-t39 * pkin(5) - t65 * pkin(13) + t67 * t53 + t23) - t30 * mrSges(7,2) - t44 * t42;
t124 = m(6) * t23 - t39 * mrSges(6,1) + t40 * mrSges(6,2) - t66 * t51 + t67 * t52 - t126;
t68 = -t85 * mrSges(5,2) + t79 * mrSges(5,3);
t15 = m(5) * t134 + t71 * mrSges(5,1) - t49 * mrSges(5,3) - t80 * t59 + t85 * t68 - t124;
t81 = -t93 * mrSges(4,2) + t86 * mrSges(4,3);
t82 = t93 * mrSges(4,1) - t87 * mrSges(4,3);
t10 = m(4) * t136 - t72 * mrSges(4,1) + t73 * mrSges(4,2) + t113 * t12 + t119 * t15 - t86 * t81 + t87 * t82;
t125 = m(5) * t33 - t48 * mrSges(5,1) + t49 * mrSges(5,2) + t112 * t14 + t118 * t13 - t79 * t68 + t80 * t69;
t74 = -t86 * mrSges(4,1) + t87 * mrSges(4,2);
t11 = m(4) * t162 + t88 * mrSges(4,1) - t73 * mrSges(4,3) - t87 * t74 + t93 * t81 - t125;
t9 = m(4) * t142 - t88 * mrSges(4,2) + t72 * mrSges(4,3) - t113 * t15 + t119 * t12 + t86 * t74 - t93 * t82;
t131 = t120 * t11 + t114 * t9;
t138 = (-mrSges(3,1) * t121 + mrSges(3,2) * t115) * t147 ^ 2;
t98 = -t105 * mrSges(3,2) + mrSges(3,3) * t139;
t4 = m(3) * (-g(3) * t149 + t133) - t102 * mrSges(3,3) + t104 * mrSges(3,1) - t115 * t138 + t105 * t98 - t107 * t10 + t131 * t109;
t97 = t105 * mrSges(3,1) - mrSges(3,3) * t140;
t6 = m(3) * (-t108 * t99 - t158) + t102 * mrSges(3,2) - t103 * mrSges(3,1) + t109 * t10 + t131 * t107 + (t115 * t97 - t121 * t98) * t147;
t8 = m(3) * (-g(3) * t150 + t155) + t103 * mrSges(3,3) - t104 * mrSges(3,2) + t121 * t138 - t105 * t97 + t120 * t9 - t114 * t11;
t144 = t110 * t6 + t4 * t149 + t8 * t150;
t2 = m(2) * t132 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t115 * t4 + t121 * t8;
t1 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t108 * t6 + (t115 * t8 + t121 * t4) * t110;
t3 = [-m(1) * g(1) - t116 * t1 + t122 * t2, t2, t8, t9, t12, t14, t17; -m(1) * g(2) + t122 * t1 + t116 * t2, t1, t4, t11, t15, t13, t16; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t6, t10, t125, t124, -t126;];
f_new  = t3;
