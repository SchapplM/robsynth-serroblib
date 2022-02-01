% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:51
% EndTime: 2022-01-20 10:05:55
% DurationCPUTime: 0.92s
% Computational Cost: add. (1170->197), mult. (1884->277), div. (0->0), fcn. (1162->14), ass. (0->136)
t130 = pkin(1) * qJD(2);
t117 = qJD(1) * t130;
t88 = sin(qJ(2));
t125 = qJDD(1) * t88;
t91 = cos(qJ(2));
t163 = pkin(1) * t125 + t91 * t117;
t152 = t91 * pkin(1);
t68 = qJDD(1) * t152;
t76 = qJDD(1) + qJDD(2);
t33 = t76 * pkin(2) - t88 * t117 + t68;
t84 = sin(pkin(8));
t86 = cos(pkin(8));
t17 = t163 * t86 + t84 * t33;
t79 = qJD(1) + qJD(2);
t12 = t76 * qJ(4) + t79 * qJD(4) + t17;
t83 = sin(pkin(9));
t85 = cos(pkin(9));
t8 = -t85 * qJDD(3) + t83 * t12;
t155 = t8 * t83;
t9 = t83 * qJDD(3) + t85 * t12;
t162 = t9 * t85 + t155;
t82 = qJ(1) + qJ(2);
t72 = sin(t82);
t73 = cos(t82);
t161 = g(1) * t72 - g(2) * t73;
t131 = pkin(1) * qJD(1);
t121 = t91 * t131;
t122 = t88 * t131;
t51 = t84 * t122;
t38 = t86 * t121 - t51;
t127 = qJD(4) - t38;
t142 = t85 * t79;
t48 = -qJD(5) + t142;
t126 = qJD(5) + t48;
t42 = t79 * pkin(2) + t121;
t52 = t86 * t122;
t26 = t84 * t42 + t52;
t22 = t79 * qJ(4) + t26;
t18 = -t85 * qJD(3) + t83 * t22;
t151 = t18 * t83;
t19 = t83 * qJD(3) + t85 * t22;
t160 = -t126 * t19 - t79 * t151;
t159 = pkin(1) * t88;
t158 = pkin(2) * t72;
t71 = pkin(8) + t82;
t60 = sin(t71);
t157 = g(1) * t60;
t154 = t86 * pkin(2);
t89 = sin(qJ(1));
t153 = t89 * pkin(1);
t57 = t84 * t159;
t103 = t86 * t91 * t130 - qJD(2) * t57;
t32 = qJD(4) + t103;
t150 = t32 * t79;
t87 = sin(qJ(5));
t149 = t32 * t87;
t143 = t85 * t76;
t46 = -qJDD(5) + t143;
t148 = t46 * t85;
t147 = t76 * t87;
t90 = cos(qJ(5));
t146 = t76 * t90;
t75 = t79 ^ 2;
t77 = t83 ^ 2;
t145 = t77 * t75;
t144 = t83 * t87;
t141 = t85 * t87;
t140 = t85 * t90;
t139 = t86 * t88;
t138 = t87 * t46;
t137 = t87 * t90;
t136 = t90 * t46;
t67 = pkin(2) + t152;
t135 = pkin(1) * t139 + t84 * t67;
t134 = g(1) * t73 + g(2) * t72;
t133 = t85 ^ 2 + t77;
t81 = t90 ^ 2;
t132 = t87 ^ 2 - t81;
t129 = qJD(5) * t87;
t128 = qJD(5) * t90;
t61 = cos(t71);
t66 = pkin(2) * t73;
t123 = t61 * pkin(3) + t60 * qJ(4) + t66;
t120 = t79 * t128;
t119 = t48 * t129;
t16 = -t163 * t84 + t86 * t33;
t100 = qJDD(4) - t16;
t13 = -t76 * pkin(3) + t100;
t116 = -g(2) * t61 - t13;
t115 = t133 * t76;
t25 = t86 * t42 - t51;
t114 = t46 - t143;
t113 = t46 + t143;
t112 = t86 * t67 - t57;
t111 = t79 * t126;
t110 = t127 * t90;
t109 = qJD(1) * (-qJD(2) + t79);
t108 = qJD(2) * (-qJD(1) - t79);
t106 = t68 + t161;
t105 = qJD(4) - t25;
t104 = t19 * t85 + t151;
t102 = -t85 * pkin(4) - t83 * pkin(7) - pkin(3);
t101 = -t60 * pkin(3) + t61 * qJ(4) - t158;
t39 = t102 - t154;
t59 = t84 * pkin(2) + qJ(4);
t99 = -t59 * t141 + t90 * t39;
t14 = t102 * t79 + t105;
t27 = t60 * t141 + t61 * t90;
t29 = -t61 * t141 + t60 * t90;
t7 = t102 * t76 + t100;
t98 = -g(1) * t27 - g(2) * t29 + (t87 * t7 + t90 * t9 + (t90 * t14 - t87 * t19) * qJD(5)) * t85 + t90 * t155;
t28 = -t60 * t140 + t61 * t87;
t30 = t61 * t140 + t60 * t87;
t97 = -g(1) * t28 - g(2) * t30 + t128 * t151 + t8 * t144;
t96 = g(3) * t83 - t126 * t14 - t9;
t95 = -g(1) * t61 - g(2) * t60 + t162;
t94 = t127 * t87 + t59 * t128;
t93 = -t48 ^ 2 - t145;
t92 = cos(qJ(1));
t74 = t92 * pkin(1);
t62 = -pkin(3) - t154;
t43 = t85 * t157;
t41 = t83 * t129 * t142;
t37 = (t84 * t91 + t139) * t130;
t36 = t84 * t121 + t52;
t35 = -pkin(3) - t112;
t34 = qJ(4) + t135;
t24 = (-0.2e1 * t87 * t120 + t76 * t81) * t77;
t23 = t102 - t112;
t21 = -t79 * pkin(3) + t105;
t20 = 0.2e1 * (t132 * t79 * qJD(5) - t76 * t137) * t77;
t11 = (t113 * t87 + (t48 + t142) * t128) * t83;
t10 = t41 + (-t113 * t90 + t119) * t83;
t5 = t90 * t7;
t2 = -t87 * t9 + t5 + (-t87 * t14 - t90 * t19) * qJD(5);
t1 = [qJDD(1), g(1) * t89 - g(2) * t92, g(1) * t92 + g(2) * t89, t76, (t88 * t108 + t76 * t91) * pkin(1) + t106, ((-qJDD(1) - t76) * t88 + t91 * t108) * pkin(1) + t134, t17 * t135 + t26 * t103 + t16 * t112 - t25 * t37 - g(1) * (-t153 - t158) - g(2) * (t66 + t74), t43 + (-t35 * t76 - t37 * t79 + t116) * t85, t34 * t115 + t133 * t150 + t95, t13 * t35 + t21 * t37 - g(1) * (t101 - t153) - g(2) * (t74 + t123) + t162 * t34 + t104 * t32, t24, t20, t10, t11, t148, -(-t23 * t129 + t90 * t37) * t48 - t23 * t136 + (-(-t34 * t128 - t149) * t48 + t34 * t138 - t2) * t85 + (t79 * t149 + (t120 + t147) * t34) * t77 + t97, (t32 * t140 + t87 * t37) * t48 + (t34 * t140 + t87 * t23) * t46 + (t34 * t76 + t150) * t90 * t77 + (t90 * t23 * t48 + (-t151 + (-t48 * t85 - t77 * t79) * t34) * t87) * qJD(5) + t98; 0, 0, 0, t76, t109 * t159 + t106, (t91 * t109 - t125) * pkin(1) + t134, t25 * t36 - t26 * t38 + (t16 * t86 + t17 * t84 + t161) * pkin(2), t43 + (t36 * t79 - t62 * t76 + t116) * t85, t127 * t79 * t133 + t59 * t115 + t95, t13 * t62 - t21 * t36 - g(1) * t101 - g(2) * t123 + (t127 * t19 + t9 * t59) * t85 + (t127 * t18 + t8 * t59) * t83, t24, t20, t10, t11, t148, -t99 * t46 - t2 * t85 + (t39 * t129 + t90 * t36 + t85 * t94) * t48 + (t59 * t147 + t79 * t94) * t77 + t97, (t59 * t140 + t87 * t39) * t46 + (t85 * t110 - t87 * t36) * t48 + (-t18 * t144 + t48 * t99) * qJD(5) + (t59 * t146 + (-t59 * t129 + t110) * t79) * t77 + t98; 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, -t8 * t85 + t9 * t83 - g(3), 0, 0, 0, 0, 0, (t114 * t87 + (t48 - t142) * t128) * t83, t41 + (t114 * t90 - t119) * t83; 0, 0, 0, 0, 0, 0, 0, -t143, -t133 * t75, -t104 * t79 - t116 - t157, 0, 0, 0, 0, 0, t87 * t93 - t136, t90 * t93 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137 * t145, -t132 * t145, (-t87 * t111 + t146) * t83, (-t90 * t111 - t147) * t83, -t46, -g(1) * t29 + g(2) * t27 + t160 * t90 + t96 * t87 + t5, g(1) * t30 - g(2) * t28 + t96 * t90 + (-t160 - t7) * t87;];
tau_reg = t1;
