% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:37
% EndTime: 2021-01-15 15:32:42
% DurationCPUTime: 1.30s
% Computational Cost: add. (1436->156), mult. (3466->238), div. (0->0), fcn. (3690->6), ass. (0->139)
t132 = sin(qJ(3));
t131 = sin(pkin(8));
t134 = cos(qJ(3));
t204 = t131 * t134;
t209 = cos(pkin(8));
t105 = t209 * t132 + t204;
t243 = t105 / 0.2e1;
t135 = cos(qJ(2));
t142 = -t131 * t132 + t209 * t134;
t242 = t142 * t135;
t133 = sin(qJ(2));
t219 = t133 / 0.2e1;
t82 = t105 * t133;
t84 = t142 * t133;
t224 = t84 * t142 / 0.2e1 + t82 * t243;
t145 = t219 + t224;
t83 = t135 * t105;
t144 = -t135 * t133 + t242 * t84 + t82 * t83;
t228 = t144 * qJD(1);
t241 = t145 * qJD(4) + t228;
t139 = t219 - t224;
t240 = -qJD(4) * t139 - t228;
t102 = t105 ^ 2;
t235 = t142 ^ 2 + t102;
t239 = qJD(4) * t235;
t238 = t235 * qJD(2);
t237 = t131 / 0.2e1;
t236 = qJ(4) + pkin(6);
t118 = t236 * t134;
t108 = t209 * t118;
t117 = t236 * t132;
t206 = t131 * t117;
t234 = t108 - t206;
t72 = t209 * t117 + t131 * t118;
t153 = t72 * t105 + t142 * t234;
t231 = qJD(4) * t153;
t229 = t139 * qJD(2);
t227 = t144 * qJD(2);
t226 = t145 * qJD(2);
t223 = -qJD(1) * t139 + qJD(2) * t153;
t156 = t108 / 0.2e1;
t122 = t131 * pkin(3) + qJ(5);
t221 = -t122 / 0.2e1;
t124 = -t209 * pkin(3) - pkin(4);
t220 = -t124 / 0.2e1;
t218 = t105 * pkin(4);
t217 = t132 * pkin(3);
t164 = t83 / 0.2e1;
t160 = t135 * t209;
t158 = -t160 / 0.2e1;
t172 = t135 * t204;
t186 = -t172 / 0.2e1 + t132 * t158;
t47 = t164 + t186;
t195 = t47 * qJD(1);
t216 = -t105 * qJD(4) + t195;
t215 = -qJD(3) * t234 - t195;
t214 = qJD(3) * pkin(3);
t208 = t142 * qJ(5);
t126 = -t134 * pkin(3) - pkin(2);
t54 = -pkin(4) * t142 - t105 * qJ(5) + t126;
t55 = -t208 + t217 + t218;
t13 = t54 * t105 - t142 * t55;
t207 = t13 * qJD(2);
t202 = t135 * t132;
t127 = t217 / 0.2e1;
t19 = t127 + (pkin(4) / 0.2e1 + t220) * t105 - (qJ(5) / 0.2e1 + t122 / 0.2e1) * t142;
t199 = t19 * qJD(2);
t171 = -t209 / 0.2e1;
t138 = t105 * t171 + t142 * t237;
t37 = (-t132 / 0.2e1 + t138) * pkin(3);
t198 = t37 * qJD(2);
t44 = t126 * t105 - t142 * t217;
t197 = t44 * qJD(2);
t165 = -t83 / 0.2e1;
t157 = t160 / 0.2e1;
t187 = t172 / 0.2e1 + t132 * t157;
t46 = t165 + t187;
t196 = t46 * qJD(1);
t192 = t72 * qJD(3);
t191 = t82 * qJD(3);
t190 = t84 * qJD(3);
t162 = t202 / 0.2e1;
t189 = t131 * t162 + t134 * t158;
t163 = -t202 / 0.2e1;
t188 = t131 * t163 + t134 * t157;
t185 = qJD(2) * t134;
t184 = t102 * qJD(2);
t183 = t142 * qJD(2);
t93 = t142 * qJD(3);
t182 = t142 * qJD(4);
t181 = t105 * qJD(2);
t180 = t105 * qJD(5);
t121 = -t132 ^ 2 + t134 ^ 2;
t179 = t121 * qJD(2);
t178 = t132 * qJD(3);
t177 = t133 * qJD(2);
t176 = t134 * qJD(3);
t175 = t135 * qJD(2);
t174 = pkin(2) * t132 * qJD(2);
t173 = pkin(2) * t185;
t170 = t105 * t177;
t169 = t142 * t181;
t168 = t132 * t185;
t167 = t242 / 0.2e1;
t166 = -t242 / 0.2e1;
t161 = t234 * t242 + t83 * t72;
t136 = -t135 * t55 / 0.2e1;
t143 = t83 * t220 + t221 * t242;
t2 = t136 + t143;
t5 = t54 * t55;
t154 = t2 * qJD(1) + t5 * qJD(2);
t152 = (t83 * t105 + t142 * t242) * qJD(2);
t12 = t126 * t217;
t141 = t83 * t171 + t237 * t242;
t3 = (t162 + t141) * pkin(3);
t151 = -t3 * qJD(1) + t12 * qJD(2);
t14 = -t55 * t105 - t142 * t54;
t48 = t166 + t188;
t150 = -t48 * qJD(1) + t14 * qJD(2);
t45 = t105 * t217 + t126 * t142;
t49 = t167 + t189;
t147 = -t49 * qJD(1) + t45 * qJD(2);
t69 = t156 - t108 / 0.2e1;
t146 = t69 * qJD(2) + t122 * qJD(3);
t95 = t105 * qJD(3);
t53 = t164 + t187;
t52 = t165 + t186;
t51 = t166 + t189;
t50 = t167 + t188;
t42 = t47 * qJD(3);
t40 = t47 * qJD(2);
t38 = 0.2e1 * t156 - t206;
t36 = t138 * pkin(3) + t127;
t28 = t52 * qJD(3) - t142 * t177;
t27 = t52 * qJD(2) - t190;
t20 = -t142 * t221 + t124 * t243 + t127 - t208 / 0.2e1 + t218 / 0.2e1;
t4 = (t141 + t163) * pkin(3);
t1 = t136 - t143;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, 0, 0, 0, t227; 0, 0, -t177, -t175, 0, 0, 0, 0, 0, -t134 * t177 - t135 * t178, t132 * t177 - t135 * t176, t28, t51 * qJD(3) + t170, t152, (t133 * t126 + t161) * qJD(2) + t4 * qJD(3) + t241, t28, t152, t50 * qJD(3) - t170, (t133 * t54 + t161) * qJD(2) + t1 * qJD(3) + t53 * qJD(5) + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132 * t175 - t133 * t176, t133 * t178 - t134 * t175, t27, t51 * qJD(2) + t191, 0, t4 * qJD(2) + (-t131 * t82 - t209 * t84) * t214, t27, 0, t50 * qJD(2) - t191, t1 * qJD(2) + (-t82 * t122 + t84 * t124) * qJD(3) + t84 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, 0, 0, 0, t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * qJD(2) + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t49 * qJD(3), 0, -t3 * qJD(3) + t240, -t42, 0, -t48 * qJD(3), t2 * qJD(3) - t46 * qJD(5) + t240; 0, 0, 0, 0, t132 * t176, t121 * qJD(3), 0, 0, 0, -pkin(2) * t178, -pkin(2) * t176, t44 * qJD(3), t45 * qJD(3), t239, t12 * qJD(3) + t231, t13 * qJD(3) + t142 * t180, t239, t14 * qJD(3) + t102 * qJD(5), t5 * qJD(3) - t54 * t180 + t231; 0, 0, 0, 0, t168, t179, t176, -t178, 0, -pkin(6) * t176 - t174, pkin(6) * t178 - t173, t197 + t215, t147 + t192, (-t105 * t131 - t142 * t209) * t214, (-t131 * t72 - t209 * t234) * t214 + t36 * qJD(4) + t151, t207 + t215, (-t122 * t105 + t124 * t142) * qJD(3) + qJD(5) * t142, t150 - t192, (-t122 * t72 + t124 * t234) * qJD(3) + t20 * qJD(4) + t38 * qJD(5) + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t36 * qJD(3) + t223, 0, t238, 0, t20 * qJD(3) + t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t93, t184, t38 * qJD(3) - t181 * t54 - t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t49 * qJD(2), 0, t3 * qJD(2), t40, 0, t48 * qJD(2), -t2 * qJD(2); 0, 0, 0, 0, -t168, -t179, 0, 0, 0, t174, t173, -t197 + t216, -t147 - t182, 0, t37 * qJD(4) - t151, -t207 + t216, 0, -t150 + t182, -t19 * qJD(4) + t69 * qJD(5) - t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t122 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t183, 0, t198, -t181, 0, t183, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, 0, 0, 0, t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t93, -t238, -t37 * qJD(3) - t223, t95, -t238, -t93, t19 * qJD(3) - t180 - t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, t183, 0, -t198, t181, 0, -t183, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, 0, -t184, t196 - t69 * qJD(3) + (qJD(2) * t54 + qJD(4)) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
