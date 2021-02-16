% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:47
% EndTime: 2021-01-15 16:33:53
% DurationCPUTime: 1.77s
% Computational Cost: add. (1989->177), mult. (4551->258), div. (0->0), fcn. (4758->6), ass. (0->159)
t191 = qJD(3) + qJD(4);
t236 = cos(qJ(4));
t190 = t236 * pkin(3);
t150 = t190 + pkin(4);
t238 = -t150 / 0.2e1;
t169 = t190 / 0.2e1 + t238;
t158 = cos(qJ(2));
t157 = cos(qJ(3));
t184 = t236 * t157;
t154 = sin(qJ(4));
t155 = sin(qJ(3));
t219 = t154 * t155;
t162 = -t184 / 0.2e1 + t219 / 0.2e1;
t166 = t184 - t219;
t242 = t166 / 0.2e1;
t75 = (t242 + t162) * t158;
t204 = t75 * qJD(1);
t185 = t236 * t155;
t218 = t154 * t157;
t135 = t185 + t218;
t246 = pkin(6) + pkin(7);
t143 = t246 * t155;
t144 = t246 * t157;
t253 = t236 * t143 + t154 * t144;
t60 = -t135 * qJ(5) - t253;
t259 = -t191 * t60 - t204;
t222 = t166 * qJ(5);
t138 = t236 * t144;
t221 = t154 * t143;
t252 = -t138 + t221;
t58 = t252 - t222;
t258 = -t58 / 0.2e1;
t156 = sin(qJ(2));
t112 = t135 * t156;
t227 = t112 * t58;
t256 = t227 / 0.2e1;
t255 = t191 * t253;
t113 = t158 * t135;
t254 = -t113 / 0.2e1;
t89 = t191 * t135;
t151 = -t157 * pkin(3) - pkin(2);
t200 = t166 * qJD(2);
t251 = -t151 * t200 + t204;
t231 = t155 * pkin(3);
t73 = t135 * t231 + t151 * t166;
t250 = -t73 * qJD(2) + t204;
t198 = t135 * qJD(2);
t163 = -t218 / 0.2e1 - t185 / 0.2e1;
t240 = t135 / 0.2e1;
t74 = (t240 + t163) * t158;
t205 = t74 * qJD(1);
t249 = -t151 * t198 + t205;
t72 = t151 * t135 - t166 * t231;
t248 = -t72 * qJD(2) + t205;
t247 = t135 ^ 2;
t245 = t58 * pkin(4);
t241 = -t166 / 0.2e1;
t239 = -t138 / 0.2e1;
t237 = t156 / 0.2e1;
t235 = pkin(3) * t154;
t114 = t166 * t156;
t234 = t114 * pkin(4);
t233 = t166 * pkin(4);
t232 = t135 * pkin(4);
t224 = -qJD(5) * t166 + t204;
t125 = t135 * qJD(5);
t223 = t205 - t125;
t107 = t151 - t233;
t81 = t107 * t135;
t108 = t231 + t232;
t217 = t158 * t108;
t115 = t166 * t158;
t21 = t112 * t113 + t114 * t115 - t158 * t156;
t216 = t21 * qJD(1);
t33 = -t108 * t166 + t81;
t215 = t33 * qJD(2);
t80 = t107 * t166;
t34 = t108 * t135 + t80;
t214 = t34 * qJD(2);
t35 = t166 * t232 - t81;
t213 = t35 * qJD(2);
t36 = -t247 * pkin(4) - t80;
t212 = t36 * qJD(2);
t180 = t154 * t242;
t44 = (t238 - pkin(4) / 0.2e1) * t135 + (t180 - t155 / 0.2e1) * pkin(3);
t211 = t44 * qJD(2);
t165 = pkin(4) / 0.2e1 + t169;
t48 = t165 * t166;
t210 = t48 * qJD(2);
t209 = t58 * qJD(4);
t132 = t166 ^ 2;
t63 = t132 - t247;
t208 = t63 * qJD(2);
t94 = t132 + t247;
t203 = t94 * qJD(2);
t97 = t239 + t138 / 0.2e1;
t202 = t97 * qJD(2);
t201 = qJD(2) * t157;
t199 = t166 * qJD(4);
t197 = t135 * qJD(4);
t148 = -t155 ^ 2 + t157 ^ 2;
t196 = t148 * qJD(2);
t195 = t155 * qJD(3);
t194 = t156 * qJD(2);
t193 = t157 * qJD(3);
t192 = t158 * qJD(2);
t189 = pkin(2) * t155 * qJD(2);
t188 = pkin(2) * t201;
t187 = qJD(4) * t235;
t186 = pkin(4) * t198;
t181 = t155 * t201;
t178 = t236 * qJD(3);
t177 = t236 * qJD(4);
t174 = pkin(3) * t177;
t164 = t150 * t254 + t115 * t235 / 0.2e1;
t1 = t217 / 0.2e1 - (t58 / 0.2e1 + t258) * t112 + t164;
t9 = t107 * t108;
t172 = -t1 * qJD(1) + t9 * qJD(2);
t10 = pkin(4) * t81;
t161 = -t227 / 0.2e1 + t256;
t4 = (t254 + t113 / 0.2e1) * pkin(4) + t161;
t171 = t4 * qJD(1) + t10 * qJD(2);
t18 = -t60 * t135 - t166 * t58;
t168 = -t112 * t240 + t114 * t241;
t31 = t237 + t168;
t170 = -t31 * qJD(1) + t18 * qJD(2);
t116 = (t190 - t150) * t235;
t30 = t165 * t114;
t159 = t169 * t58;
t6 = -t245 / 0.2e1 - t159;
t160 = -t30 * qJD(1) - t6 * qJD(2) - t116 * qJD(3);
t65 = 0.2e1 * t239 + t221;
t95 = t166 * t198;
t93 = t97 * qJD(3);
t92 = t97 * qJD(4);
t88 = t191 * t166;
t79 = qJD(3) * t235 - t202;
t78 = pkin(3) * t178;
t77 = t163 * t158 + t254;
t76 = (t162 + t241) * t158;
t68 = t75 * qJD(2);
t66 = t74 * qJD(2);
t57 = -t191 * t235 + t202;
t56 = (-t178 - t177) * pkin(3);
t47 = -t233 / 0.2e1 + t169 * t166;
t46 = t65 - t222;
t43 = pkin(3) * t180 + t135 * t238 + t231 / 0.2e1 + t232 / 0.2e1;
t32 = t237 - t168;
t29 = -t234 / 0.2e1 + t169 * t114;
t25 = t191 * t74;
t24 = t191 * t75;
t23 = t77 * qJD(2) - t191 * t114;
t22 = t76 * qJD(2) + t191 * t112;
t20 = t135 * t194 + t191 * t76;
t19 = -t166 * t194 + t191 * t77;
t5 = t245 / 0.2e1 - t159;
t3 = 0.2e1 * t254 * pkin(4) + t161;
t2 = t112 * t258 - t217 / 0.2e1 + t164 + t256;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * qJD(2); 0, 0, -t194, -t192, 0, 0, 0, 0, 0, -t157 * t194 - t158 * t195, t155 * t194 - t158 * t193, 0, 0, 0, 0, 0, t19, t20, t19, t20, (t113 * t135 + t115 * t166) * qJD(2), t216 + (t156 * t107 - t113 * t60 - t115 * t58) * qJD(2) + t2 * qJD(3) + t3 * qJD(4) + t32 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155 * t192 - t156 * t193, t156 * t195 - t157 * t192, 0, 0, 0, 0, 0, t23, t22, t23, t22, 0, t2 * qJD(2) + (-t112 * t235 - t114 * t150) * qJD(3) + t29 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t22, t23, t22, 0, t3 * qJD(2) + t29 * qJD(3) - qJD(4) * t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t24, -t25, -t24, 0, -t1 * qJD(3) + t4 * qJD(4) - t31 * qJD(5) - t216; 0, 0, 0, 0, t155 * t193, t148 * qJD(3), 0, 0, 0, -pkin(2) * t195, -pkin(2) * t193, t166 * t89, t191 * t63, 0, 0, 0, t72 * qJD(3) + t151 * t197, t73 * qJD(3) + t151 * t199, t33 * qJD(3) - t35 * qJD(4), t34 * qJD(3) - t36 * qJD(4), t94 * qJD(5), qJD(3) * t9 + qJD(4) * t10 + qJD(5) * t18; 0, 0, 0, 0, t181, t196, t193, -t195, 0, -pkin(6) * t193 - t189, pkin(6) * t195 - t188, t95, t208, t88, -t89, 0, qJD(3) * t252 + t65 * qJD(4) - t248, -t250 + t255, t58 * qJD(3) + t46 * qJD(4) - t205 + t215, t214 + t259, (-t135 * t235 - t150 * t166) * qJD(3) + t47 * qJD(4), (t58 * t150 + t235 * t60) * qJD(3) + t5 * qJD(4) + t43 * qJD(5) + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t208, t88, -t89, 0, t65 * qJD(3) + qJD(4) * t252 - t249, -t251 + t255, t46 * qJD(3) - t205 + t209 - t213, -t212 + t259, -pkin(4) * t199 + t47 * qJD(3), pkin(4) * t209 + t5 * qJD(3) + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t43 * qJD(3) + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t68, t66, t68, 0, t1 * qJD(2) + t30 * qJD(4); 0, 0, 0, 0, -t181, -t196, 0, 0, 0, t189, t188, -t95, -t208, 0, 0, 0, t92 + t248, t250, t92 - t215 + t223, -t214 + t224, t48 * qJD(4), t6 * qJD(4) + t44 * qJD(5) - t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, -t174, -t187, -t174, 0, t116 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t56, t57, t56, t210, -pkin(4) * t187 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, -t200, 0, t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t68, t66, t68, 0, -t4 * qJD(2) - t30 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t208, 0, 0, 0, -t93 + t249, t251, -t93 + t213 + t223, t212 + t224, -t48 * qJD(3), -pkin(4) * t125 - t6 * qJD(3) - t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t78, t79, t78, -t210, t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, -t200, 0, -t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88, -t203, pkin(4) * t197 - t44 * qJD(3) - t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t200, 0, -t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t200, 0, t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
