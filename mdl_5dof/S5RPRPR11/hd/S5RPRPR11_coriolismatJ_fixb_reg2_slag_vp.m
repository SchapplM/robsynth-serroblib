% Calculate inertial parameters regressor of coriolis matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR11_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:28:00
% EndTime: 2019-12-31 18:28:04
% DurationCPUTime: 2.56s
% Computational Cost: add. (3571->182), mult. (6770->218), div. (0->0), fcn. (7616->6), ass. (0->147)
t187 = qJD(3) - qJD(5);
t149 = sin(qJ(5));
t148 = cos(pkin(8));
t239 = pkin(6) + qJ(2);
t137 = t239 * t148;
t150 = sin(qJ(3));
t147 = sin(pkin(8));
t173 = t239 * t147;
t245 = cos(qJ(3));
t104 = t150 * t137 + t245 * t173;
t175 = t245 * t147;
t227 = t150 * t148;
t131 = t175 + t227;
t164 = -t131 * pkin(7) + t104;
t244 = cos(qJ(5));
t105 = t245 * t137 - t150 * t173;
t129 = t150 * t147 - t245 * t148;
t75 = t129 * pkin(7) + t105;
t15 = -t149 * t75 + t244 * t164;
t276 = t187 * t15;
t170 = t149 * t129 + t244 * t131;
t240 = t170 ^ 2;
t255 = -t244 * t129 + t149 * t131;
t241 = t255 ^ 2;
t34 = t241 - t240;
t275 = t34 * qJD(1);
t269 = t149 * t164 + t244 * t75;
t274 = t187 * t269;
t264 = t170 * qJD(1);
t273 = t255 * t264;
t125 = t129 ^ 2;
t250 = t131 ^ 2;
t76 = t250 - t125;
t272 = t76 * qJD(1);
t271 = t76 * qJD(3);
t203 = t255 * qJD(5);
t268 = t255 * qJD(3) - t203;
t266 = qJD(2) * t170;
t265 = qJD(2) * t255;
t254 = t125 + t250;
t263 = t254 * qJD(1);
t262 = t254 * qJD(2);
t261 = t255 * qJD(1);
t200 = t170 * qJD(5);
t259 = qJD(3) * t170 - t200;
t186 = -t244 / 0.2e1;
t161 = t104 * t131 - t105 * t129;
t253 = t161 * qJD(1);
t252 = t161 * qJD(2);
t249 = t170 / 0.2e1;
t248 = -pkin(3) - pkin(4);
t247 = -t149 / 0.2e1;
t246 = t149 / 0.2e1;
t243 = pkin(3) * t131;
t242 = t129 * pkin(3);
t143 = -t148 * pkin(2) - pkin(1);
t230 = t131 * qJ(4);
t160 = -t143 + t230;
t55 = t248 * t129 + t160;
t232 = t129 * qJ(4);
t74 = t248 * t131 - t232;
t3 = t55 * t74;
t237 = t3 * qJD(1);
t9 = t15 * t170 + t255 * t269;
t236 = t9 * qJD(1);
t180 = t269 * t244;
t168 = -t180 / 0.2e1;
t38 = t180 / 0.2e1;
t11 = t38 + t168;
t234 = t11 * qJD(1);
t233 = t11 * qJD(4);
t86 = -t160 + t242;
t96 = t232 + t243;
t18 = t86 * t96;
t226 = t18 * qJD(1);
t21 = -t170 * t55 + t255 * t74;
t225 = t21 * qJD(1);
t22 = t170 * t74 + t255 * t55;
t224 = t22 * qJD(1);
t166 = t255 * t186;
t179 = t244 * t255;
t167 = t179 / 0.2e1;
t23 = t166 + t167 + (t249 - t170 / 0.2e1) * t149;
t223 = t23 * qJD(1);
t113 = t243 / 0.2e1;
t156 = t113 + t232 / 0.2e1 + t131 * pkin(4) / 0.2e1;
t133 = t149 * qJ(4) - t244 * t248;
t134 = t244 * qJ(4) + t149 * t248;
t159 = t133 * t249 - t255 * t134 / 0.2e1;
t26 = t156 + t159;
t222 = t26 * qJD(1);
t27 = -t240 - t241;
t221 = t27 * qJD(1);
t40 = t96 * t129 + t86 * t131;
t218 = t40 * qJD(1);
t41 = t86 * t129 - t96 * t131;
t217 = t41 * qJD(1);
t153 = -t227 / 0.2e1 - t175 / 0.2e1;
t155 = t170 * t186 + t247 * t255;
t48 = t153 + t155;
t216 = t48 * qJD(1);
t81 = -t179 / 0.2e1;
t52 = t81 + t167;
t213 = t52 * qJD(1);
t212 = t52 * qJD(4);
t209 = t96 * qJD(1);
t139 = t147 ^ 2 + t148 ^ 2;
t199 = qJD(3) * qJ(4);
t196 = t104 * qJD(3);
t98 = t105 * qJD(3);
t195 = t250 * qJD(1);
t194 = t129 * qJD(1);
t118 = t129 * qJD(3);
t193 = t131 * qJD(1);
t120 = t131 * qJD(3);
t192 = t131 * qJD(4);
t135 = t139 * qJ(2);
t191 = t135 * qJD(1);
t190 = t139 * qJD(1);
t189 = t149 * qJD(3);
t188 = t149 * qJD(4);
t185 = t55 * t261;
t184 = t55 * t264;
t181 = t170 * t261;
t178 = t255 * t193;
t177 = t170 * t193;
t176 = t55 * t193;
t103 = t129 * t193;
t102 = t129 * t120;
t174 = t143 * t193;
t172 = t244 * qJD(3);
t171 = t244 * qJD(4);
t169 = t269 * t186;
t6 = -t169 + t168;
t99 = t133 * t149 + t134 * t244;
t163 = -t6 * qJD(1) + t99 * qJD(3);
t138 = -qJD(5) * t244 + t172;
t136 = t187 * t149;
t122 = t131 * qJD(2);
t95 = t113 - t243 / 0.2e1;
t51 = t52 * qJD(5);
t47 = t153 - t155;
t25 = t156 - t159;
t24 = t166 + t81 + (t246 - t247) * t170;
t10 = t11 * qJD(5);
t7 = -0.2e1 * t15 * t246 - t169 + t38;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139 * qJD(2), t135 * qJD(2), -t102, -t271, 0, t102, 0, 0, t143 * t120, -t143 * t118, t262, t252, -t102, 0, t271, 0, 0, t102, t40 * qJD(3) - t129 * t192, t262, t41 * qJD(3) + qJD(4) * t250, t18 * qJD(3) - t192 * t86 + t252, t268 * t170, -t187 * t34, 0, -t259 * t255, 0, 0, t21 * qJD(3) + t192 * t255 + t200 * t55, t22 * qJD(3) + t170 * t192 - t203 * t55, t27 * qJD(2), t9 * qJD(2) + t3 * qJD(3) + t192 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, t191, 0, 0, 0, 0, 0, 0, 0, 0, t263, t253, 0, 0, 0, 0, 0, 0, 0, t263, 0, t95 * qJD(3) + t253, 0, 0, 0, 0, 0, 0, 0, 0, t221, t25 * qJD(3) + t47 * qJD(4) + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t272, -t118, t103, -t120, 0, -t98 + t174, -t143 * t194 + t196, 0, 0, -t103, -t118, t272, 0, t120, t103, -t98 + t218, (-t230 + t242) * qJD(3) - t129 * qJD(4), -t196 + t217, t226 + t95 * qJD(2) + (-t105 * pkin(3) - t104 * qJ(4)) * qJD(3) + t105 * qJD(4), t181, -t275, -t268, -t273, -t259, 0, t225 - t274, t224 - t276, (t133 * t255 + t134 * t170) * qJD(3) + t24 * qJD(4), t237 + t25 * qJD(2) + (-t133 * t269 - t134 * t15) * qJD(3) + t7 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t118, t195, -t193 * t86 + t98, 0, 0, 0, 0, 0, 0, t178, t177, t24 * qJD(3) + t51, t47 * qJD(2) + t7 * qJD(3) + t10 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, t275, t268, t273, t259, 0, t184 + t274, -t185 + t276, t212, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, -t191, 0, 0, 0, 0, 0, 0, t120, -t118, -t263, -t253, 0, 0, 0, 0, 0, 0, t120, -t263, t118, qJD(3) * t96 - t192 - t253, 0, 0, 0, 0, 0, 0, t259, -t268, -t221, t26 * qJD(3) + t48 * qJD(4) - t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t194, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, t194, t209, 0, 0, 0, 0, 0, 0, t264, -t261, 0, t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t264, t261, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t272, 0, -t103, 0, 0, -t122 - t174, (qJD(1) * t143 + qJD(2)) * t129, 0, 0, t103, 0, -t272, 0, 0, -t103, -t122 - t218, 0, -t129 * qJD(2) - t217, -qJD(2) * t96 - t226, -t181, t275, 0, t273, 0, 0, -t225 - t266, -t224 + t265, -t23 * qJD(4), -t26 * qJD(2) - t6 * qJD(4) - t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, t194, 0, 0, 0, 0, 0, 0, 0, 0, -t193, 0, -t194, -t209, 0, 0, 0, 0, 0, 0, -t264, t261, 0, -t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), 0, 0, 0, 0, 0, 0, t134 * qJD(5) + t188, -t133 * qJD(5) + t171, 0, t99 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t199, 0, 0, 0, 0, 0, 0, t189, t172, -t223, t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134 * t187, -t133 * t187, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, -t195, (qJD(1) * t86 + qJD(2)) * t131, 0, 0, 0, 0, 0, 0, -t178, -t177, t23 * qJD(3) + t51, -t48 * qJD(2) + t6 * qJD(3) + t10 - t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t199, 0, 0, 0, 0, 0, 0, -t136, -t138, t223, -t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t138, t213, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, -t275, 0, -t273, 0, 0, -t184 + t266, t185 - t265, -t212, -t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, -t261, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * qJD(3) - t188, t133 * qJD(3) - t171, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t172, -t213, -t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
