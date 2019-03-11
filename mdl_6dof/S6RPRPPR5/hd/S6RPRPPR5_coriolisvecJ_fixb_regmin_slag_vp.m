% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:38
% EndTime: 2019-03-09 02:51:47
% DurationCPUTime: 3.21s
% Computational Cost: add. (3326->323), mult. (8713->428), div. (0->0), fcn. (6638->8), ass. (0->164)
t157 = sin(pkin(9));
t159 = cos(pkin(9));
t162 = sin(qJ(3));
t236 = cos(qJ(3));
t132 = t236 * t157 + t162 * t159;
t242 = t132 * qJD(1);
t110 = qJD(6) + t242;
t198 = t236 * t159;
t190 = qJD(1) * t198;
t144 = qJD(3) * t190;
t215 = t162 * t157;
t197 = qJD(1) * t215;
t103 = qJD(3) * t197 - t144;
t156 = sin(pkin(10));
t158 = cos(pkin(10));
t161 = sin(qJ(6));
t163 = cos(qJ(6));
t129 = t163 * t156 + t161 * t158;
t225 = t110 * t129;
t243 = -t161 * t156 + t163 * t158;
t247 = -t243 * t103 - t225 * t110;
t130 = -t198 + t215;
t75 = t243 * t130;
t202 = qJD(6) * t163;
t203 = qJD(6) * t161;
t224 = -t156 * t203 + t158 * t202 + t243 * t242;
t191 = t129 * t103 - t110 * t224;
t117 = -t190 + t197;
t92 = t156 * qJD(3) - t158 * t117;
t94 = t158 * qJD(3) + t156 * t117;
t182 = t161 * t92 - t163 * t94;
t246 = t110 * t182;
t245 = t242 * qJD(3);
t232 = pkin(7) + qJ(2);
t138 = t232 * t157;
t133 = qJD(1) * t138;
t139 = t232 * t159;
t134 = qJD(1) * t139;
t83 = t236 * t133 + t162 * t134;
t244 = qJD(4) + t83;
t91 = t162 * t138 - t236 * t139;
t241 = -qJD(6) + t110;
t124 = t132 * qJD(3);
t104 = qJD(1) * t124;
t223 = qJ(4) * t104;
t233 = pkin(3) + qJ(5);
t155 = qJD(3) * qJ(4);
t84 = -t162 * t133 + t236 * t134;
t65 = -t117 * pkin(4) + t84;
t56 = qJD(5) + t155 + t65;
t240 = -t103 * t233 + (qJD(5) - t56) * t242 + t223;
t237 = t242 ^ 2;
t239 = t156 * t103 - t158 * t237;
t195 = qJD(3) * t236;
t196 = qJD(2) * t236;
t69 = (qJD(2) * t157 + qJD(3) * t139) * t162 + t138 * t195 - t159 * t196;
t23 = -t182 * qJD(6) - t104 * t243;
t238 = t117 ^ 2;
t235 = t104 * pkin(3);
t234 = t158 * pkin(8);
t231 = -pkin(8) - t233;
t222 = t103 * qJ(4);
t178 = -qJD(4) * t242 + t222;
t28 = t117 * qJD(5) + t233 * t104 + t178;
t189 = qJD(1) * t196;
t201 = qJD(1) * qJD(2);
t194 = t162 * t201;
t204 = qJD(3) * t162;
t61 = -t133 * t204 + t134 * t195 + t157 * t189 + t159 * t194;
t38 = -t103 * pkin(4) - qJD(3) * qJD(5) + t61;
t8 = t156 * t38 + t158 * t28;
t123 = t157 * t204 - t159 * t195;
t177 = t123 * qJ(4) - t132 * qJD(4);
t32 = t130 * qJD(5) + t233 * t124 + t177;
t70 = t132 * qJD(2) - t91 * qJD(3);
t45 = -t123 * pkin(4) + t70;
t13 = t156 * t45 + t158 * t32;
t151 = -t159 * pkin(2) - pkin(1);
t137 = t151 * qJD(1) + qJD(2);
t166 = -qJ(4) * t242 + t137;
t43 = t233 * t117 + t166;
t208 = pkin(4) * t242 + t244;
t50 = -t233 * qJD(3) + t208;
t18 = t156 * t50 + t158 * t43;
t221 = t117 * qJ(4);
t60 = t233 * t242 + t221;
t21 = t156 * t65 + t158 * t60;
t174 = -t132 * qJ(4) + t151;
t63 = t233 * t130 + t174;
t90 = t236 * t138 + t162 * t139;
t73 = t132 * pkin(4) + t90;
t26 = t156 * t73 + t158 * t63;
t51 = t161 * t94 + t163 * t92;
t229 = t117 * t51;
t228 = t51 * t110;
t227 = t182 * t117;
t66 = t117 * pkin(3) + t166;
t226 = t66 * t242;
t220 = t242 * t117;
t218 = t158 * t104;
t212 = t69 * qJD(3);
t211 = t70 * qJD(3);
t200 = -pkin(5) * t158 - pkin(4);
t209 = -t200 * t242 + t244;
t205 = t157 ^ 2 + t159 ^ 2;
t36 = t158 * t38;
t4 = -t103 * pkin(5) + t36 + (-pkin(8) * t104 - t28) * t156;
t5 = pkin(8) * t218 + t8;
t199 = -t161 * t5 + t163 * t4;
t17 = -t156 * t43 + t158 * t50;
t193 = t205 * qJD(1) ^ 2;
t192 = -t133 * t195 - t134 * t204 - t157 * t194 + t159 * t189;
t7 = -t156 * t28 + t36;
t187 = t8 * t156 + t7 * t158;
t186 = t161 * t4 + t163 * t5;
t10 = pkin(5) * t242 - t94 * pkin(8) + t17;
t11 = -t92 * pkin(8) + t18;
t1 = t163 * t10 - t161 * t11;
t2 = t161 * t10 + t163 * t11;
t68 = t158 * t73;
t15 = t132 * pkin(5) + t68 + (-pkin(8) * t130 - t63) * t156;
t19 = t130 * t234 + t26;
t185 = t163 * t15 - t161 * t19;
t184 = t161 * t15 + t163 * t19;
t183 = -t156 * t18 - t158 * t17;
t154 = qJD(3) * qJD(4);
t55 = -t154 - t192;
t180 = -t158 * t103 - t156 * t237;
t176 = 0.2e1 * t205 * t201;
t22 = t104 * t129 - t92 * t202 - t94 * t203;
t135 = t231 * t156;
t59 = t158 * t65;
t173 = qJD(5) * t158 + qJD(6) * t135 - t117 * pkin(5) + t59 + (-pkin(8) * t242 - t60) * t156;
t136 = t231 * t158;
t172 = qJD(5) * t156 - qJD(6) * t136 + t234 * t242 + t21;
t171 = t84 * qJD(3) - t61;
t170 = -t83 * qJD(3) - t192;
t76 = t129 * t130;
t37 = -t104 * pkin(4) - t55;
t74 = -t130 * pkin(4) - t91;
t168 = t104 * t74 + t124 * t56 + t130 * t37;
t149 = t156 * pkin(5) + qJ(4);
t105 = qJD(3) * t117;
t85 = t103 * t132;
t80 = t130 * pkin(3) + t174;
t79 = pkin(3) * t242 + t221;
t78 = -t155 - t84;
t77 = -qJD(3) * pkin(3) + t244;
t62 = t124 * pkin(3) + t177;
t49 = t178 + t235;
t48 = t200 * t130 - t91;
t44 = -t124 * pkin(4) - t69;
t42 = t158 * t45;
t34 = t92 * pkin(5) + t56;
t33 = t200 * t124 - t69;
t31 = qJD(6) * t76 - t243 * t124;
t30 = qJD(6) * t75 + t129 * t124;
t27 = t200 * t104 - t55;
t25 = -t156 * t63 + t68;
t20 = -t156 * t60 + t59;
t12 = -t156 * t32 + t42;
t9 = t124 * t234 + t13;
t6 = -t123 * pkin(5) + t42 + (-pkin(8) * t124 - t32) * t156;
t3 = [0, 0, 0, 0, 0, t176, qJ(2) * t176, -t123 * t242 - t85, t103 * t130 - t132 * t104 + t123 * t117 - t124 * t242, -t123 * qJD(3), -t124 * qJD(3), 0, t151 * t104 + t137 * t124 - t211, -t151 * t103 - t137 * t123 + t212, -t90 * t103 + t91 * t104 + t69 * t117 - t77 * t123 + t78 * t124 + t55 * t130 + t61 * t132 + t242 * t70, -t80 * t104 - t62 * t117 - t66 * t124 - t49 * t130 + t211, t80 * t103 + t66 * t123 - t49 * t132 - t242 * t62 - t212, t49 * t80 + t55 * t91 + t61 * t90 + t66 * t62 + t78 * t69 + t77 * t70, -t25 * t103 + t12 * t242 - t17 * t123 + t7 * t132 - t158 * t168 + t44 * t92, t26 * t103 + t18 * t123 - t13 * t242 - t8 * t132 + t156 * t168 + t44 * t94, -t12 * t94 - t13 * t92 + (t104 * t26 + t124 * t18 + t130 * t8) * t158 + (-t104 * t25 - t124 * t17 - t130 * t7) * t156, t17 * t12 + t18 * t13 + t7 * t25 + t8 * t26 + t37 * t74 + t56 * t44, -t182 * t30 + t22 * t76, t182 * t31 + t22 * t75 - t76 * t23 - t30 * t51, -t76 * t103 + t30 * t110 + t123 * t182 + t22 * t132, -t75 * t103 - t31 * t110 + t51 * t123 - t23 * t132, -t110 * t123 - t85 (-t161 * t9 + t163 * t6) * t110 - t185 * t103 + t199 * t132 - t1 * t123 + t33 * t51 + t48 * t23 - t27 * t75 + t34 * t31 + (-t110 * t184 - t132 * t2) * qJD(6) -(t161 * t6 + t163 * t9) * t110 + t184 * t103 - t186 * t132 + t2 * t123 - t33 * t182 + t48 * t22 + t27 * t76 + t34 * t30 + (-t1 * t132 - t110 * t185) * qJD(6); 0, 0, 0, 0, 0, -t193, -qJ(2) * t193, 0, 0, 0, 0, 0, 0.2e1 * t245, t144 + (-t117 - t197) * qJD(3), -t237 - t238, -0.2e1 * t245, t103 + t105, t235 + t222 - t78 * t117 + (-qJD(4) - t77) * t242, t117 * t92 + t239, t117 * t94 - t180 (t156 * t92 + t158 * t94) * t242 + (t156 ^ 2 + t158 ^ 2) * t104, t56 * t117 - t7 * t156 + t8 * t158 + t183 * t242, 0, 0, 0, 0, 0, t191 + t229, -t227 - t247; 0, 0, 0, 0, 0, 0, 0, t220, t237 - t238, t144 + (t117 - t197) * qJD(3), 0, 0, -t137 * t242 + t171, t137 * t117 + t170, pkin(3) * t103 - t223 + (-t78 - t84) * t242 + (t77 - t244) * t117, t79 * t117 - t171 + t226, -t66 * t117 + t242 * t79 + 0.2e1 * t154 - t170, -t61 * pkin(3) - t55 * qJ(4) - t244 * t78 - t66 * t79 - t77 * t84, t17 * t117 + t37 * t156 - t240 * t158 - t20 * t242 + t208 * t92, -t18 * t117 + t240 * t156 + t37 * t158 + t208 * t94 + t21 * t242, t20 * t94 + t21 * t92 + (qJD(5) * t94 - t18 * t242 - t7) * t158 + (qJD(5) * t92 + t17 * t242 - t8) * t156, t37 * qJ(4) + t183 * qJD(5) - t17 * t20 - t18 * t21 - t187 * t233 + t208 * t56, t182 * t225 + t22 * t243, -t22 * t129 + t182 * t224 + t225 * t51 - t23 * t243, -t227 + t247, t191 - t229, t110 * t117 -(-t161 * t135 + t163 * t136) * t103 + t149 * t23 + t27 * t129 + t1 * t117 + t209 * t51 + t224 * t34 + (t161 * t172 - t163 * t173) * t110 (t163 * t135 + t161 * t136) * t103 + t149 * t22 + t27 * t243 - t2 * t117 - t209 * t182 - t225 * t34 + (t161 * t173 + t163 * t172) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103 + t105, -t220, -qJD(3) ^ 2 - t237, t78 * qJD(3) + t226 + t61, -qJD(3) * t92 + t180, -qJD(3) * t94 + t239 (t156 * t94 - t158 * t92) * t242, -t56 * qJD(3) + (-t156 * t17 + t158 * t18) * t242 + t187, 0, 0, 0, 0, 0, -qJD(3) * t51 + t247, qJD(3) * t182 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242 * t94 - t218, t156 * t104 - t242 * t92, -t92 ^ 2 - t94 ^ 2, t17 * t94 + t18 * t92 + t37, 0, 0, 0, 0, 0, t23 - t246, t22 - t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182 * t51, t182 ^ 2 - t51 ^ 2, t22 + t228, -t23 - t246, -t103, t182 * t34 + t241 * t2 + t199, t241 * t1 + t34 * t51 - t186;];
tauc_reg  = t3;
