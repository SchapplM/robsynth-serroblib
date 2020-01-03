% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:14
% EndTime: 2019-12-31 18:49:19
% DurationCPUTime: 2.10s
% Computational Cost: add. (3577->176), mult. (6897->224), div. (0->0), fcn. (7890->6), ass. (0->142)
t133 = qJD(3) + qJD(4);
t135 = sin(pkin(8));
t136 = cos(pkin(8));
t230 = sin(qJ(3));
t232 = cos(qJ(3));
t120 = -t230 * t135 + t232 * t136;
t137 = sin(qJ(4));
t154 = t232 * t135 + t230 * t136;
t231 = cos(qJ(4));
t107 = t231 * t120 - t137 * t154;
t207 = t107 ^ 2;
t204 = t137 * t120;
t240 = t231 * t154;
t242 = t240 + t204;
t246 = t242 ^ 2;
t252 = t207 - t246;
t254 = t252 * qJD(1);
t250 = t107 * qJ(5);
t253 = -t250 / 0.2e1;
t245 = t242 * pkin(4);
t54 = -t250 + t245;
t236 = t107 / 0.2e1;
t127 = -t136 * pkin(2) - pkin(1);
t110 = -t120 * pkin(3) + t127;
t161 = -pkin(4) * t107 - qJ(5) * t242;
t45 = t110 + t161;
t251 = t45 * t107;
t249 = t107 * qJD(1);
t248 = t107 * qJD(2);
t185 = t107 * qJD(4);
t247 = t107 * qJD(3) + t185;
t218 = t45 * t242;
t244 = t242 * qJD(1);
t226 = pkin(6) + qJ(2);
t162 = t230 * t226;
t163 = t232 * t226;
t73 = (-t230 * pkin(7) - t162) * t136 + (-t232 * pkin(7) - t163) * t135;
t142 = t120 * t226;
t74 = t120 * pkin(7) + t142;
t48 = -t137 * t74 + t231 * t73;
t153 = t133 * t48;
t143 = t240 / 0.2e1;
t139 = t137 * t73;
t173 = t231 * t74;
t239 = t173 + t139;
t243 = t239 * pkin(4);
t138 = t139 / 0.2e1;
t184 = t107 * qJD(5);
t149 = t154 * qJD(1);
t148 = t154 * qJD(3);
t102 = t143 - t240 / 0.2e1;
t26 = t173 + 0.2e1 * t138;
t176 = t102 * qJD(2) - t26 * qJD(3) - qJD(4) * t239;
t238 = -pkin(4) / 0.2e1;
t237 = -qJ(5) / 0.2e1;
t227 = t137 * pkin(3);
t126 = qJ(5) + t227;
t235 = -t126 / 0.2e1;
t175 = t231 * pkin(3);
t128 = -t175 - pkin(4);
t234 = t128 / 0.2e1;
t233 = t137 / 0.2e1;
t27 = t138 - t139 / 0.2e1;
t63 = 0.2e1 * t143 + t204;
t225 = -t63 * qJD(2) + t27 * qJD(3);
t224 = -qJD(2) * t242 - t27 * qJD(4);
t223 = -qJD(3) * t239 - t26 * qJD(4);
t222 = pkin(3) * qJD(4);
t152 = t154 * pkin(3);
t53 = t152 + t54;
t3 = t45 * t53;
t220 = t3 * qJD(1);
t4 = t45 * t54;
t219 = t4 * qJD(1);
t7 = t107 * t239 - t242 * t48;
t216 = t7 * qJD(1);
t8 = -t107 * t53 + t218;
t215 = t8 * qJD(1);
t9 = -t242 * t53 - t251;
t214 = t9 * qJD(1);
t10 = -t107 * t54 + t218;
t211 = t10 * qJD(1);
t11 = -t242 * t54 - t251;
t206 = t11 * qJD(1);
t140 = (t231 * t236 + t233 * t242) * pkin(3) + t242 * t235 + t128 * t236;
t155 = t107 * t238 + t237 * t242;
t13 = t140 - t155;
t205 = t13 * qJD(1);
t202 = t27 * qJD(1);
t164 = t126 * t236;
t90 = t250 / 0.2e1;
t28 = t164 - t152 / 0.2e1 + t90 + (t234 + t238) * t242;
t201 = t28 * qJD(1);
t31 = 0.2e1 * t253 + t245;
t200 = t31 * qJD(1);
t32 = t207 + t246;
t199 = t32 * qJD(1);
t34 = t107 * t152 - t110 * t242;
t197 = t34 * qJD(1);
t35 = -t110 * t107 - t152 * t242;
t196 = t35 * qJD(1);
t193 = t63 * qJD(1);
t77 = t120 ^ 2 - t154 ^ 2;
t192 = t77 * qJD(1);
t123 = t135 ^ 2 + t136 ^ 2;
t190 = qJD(1) * t110;
t188 = t102 * qJD(1);
t86 = t102 * qJD(4);
t187 = t246 * qJD(1);
t181 = t242 * qJD(4);
t180 = t120 * qJD(1);
t119 = t120 * qJD(3);
t121 = t123 * qJ(2);
t179 = t121 * qJD(1);
t178 = t123 * qJD(1);
t130 = qJD(4) * t175;
t177 = t130 + qJD(5);
t174 = t137 * t222;
t172 = t45 * t244;
t169 = t107 * t244;
t168 = t242 * t249;
t167 = t107 * t190;
t166 = t242 * t190;
t111 = (t231 * t126 + t128 * t137) * pkin(3);
t141 = (t239 * t231 / 0.2e1 - t48 * t233) * pkin(3) - t48 * t235 + t239 * t234;
t156 = t243 / 0.2e1 + t48 * t237;
t2 = t141 + t156;
t160 = t2 * qJD(1) + t111 * qJD(3);
t51 = qJD(3) * t242 + t63 * qJD(4);
t146 = t127 * t149;
t145 = t120 * t149;
t134 = qJ(5) * qJD(5);
t129 = qJD(3) * t175;
t124 = t133 * qJ(5);
t122 = t126 * qJD(5);
t30 = t253 + t90;
t29 = t152 / 0.2e1 + t253 + t245 / 0.2e1 + t164 + t242 * t234;
t15 = qJD(3) * t227 + t202;
t14 = -t133 * t227 - t202;
t12 = t140 + t155;
t1 = t141 - t156;
t5 = [0, 0, 0, 0, 0, t123 * qJD(2), t121 * qJD(2), t120 * t148, t77 * qJD(3), 0, 0, 0, t127 * t148, t127 * t119, t247 * t242, t133 * t252, 0, 0, 0, -t34 * qJD(3) + t110 * t181, -t35 * qJD(3) + t110 * t185, t8 * qJD(3) + t10 * qJD(4) + t184 * t242, t32 * qJD(2), t9 * qJD(3) + t11 * qJD(4) + qJD(5) * t246, t7 * qJD(2) + t3 * qJD(3) + t4 * qJD(4) - qJD(5) * t218; 0, 0, 0, 0, 0, t178, t179, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, t86, t199, 0, t29 * qJD(3) + t30 * qJD(4) - t102 * qJD(5) + t216; 0, 0, 0, 0, 0, 0, 0, t145, t192, t119, -t148, 0, -t142 * qJD(3) + t146, t127 * t180 + (t135 * t163 + t136 * t162) * qJD(3), t168, t254, t247, -t51, 0, -t197 + t223, -t153 - t196, t215 + t223, (t128 * t107 - t126 * t242) * qJD(3) + t12 * qJD(4) + t184, t153 + t214, t220 + t29 * qJD(2) + (t48 * t126 + t128 * t239) * qJD(3) + t1 * qJD(4) + t26 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t254, t247, -t63 * qJD(3) - t181, 0, t166 + t176, -t153 + t167, t176 + t211, t12 * qJD(3) + t161 * qJD(4) + t184, t153 + t206, t219 + t30 * qJD(2) + t1 * qJD(3) + (t48 * qJ(5) - t243) * qJD(4) + t239 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t247, t187, -t172 - t176; 0, 0, 0, 0, 0, -t178, -t179, 0, 0, 0, 0, 0, t148, t119, 0, 0, 0, 0, 0, t51, t247, t51, -t199, -t247, -t28 * qJD(3) + t31 * qJD(4) - t63 * qJD(5) - t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t180, 0, 0, 0, 0, 0, t244, t249, t244, 0, -t249, -t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t249, t193, 0, -t249, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193; 0, 0, 0, 0, 0, 0, 0, -t145, -t192, 0, 0, 0, -t154 * qJD(2) - t146, (-qJD(1) * t127 - qJD(2)) * t120, -t168, -t254, 0, -t86, 0, t197 + t224, -t248 + t196, -t215 + t224, t13 * qJD(4), t248 - t214, t28 * qJD(2) + t2 * qJD(4) + t27 * qJD(5) - t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t180, 0, 0, 0, 0, 0, -t244, -t249, -t244, 0, t249, t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, -t130, -t174, 0, t177, t111 * qJD(4) + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, 0, t14, -t129 - t130, t14, t205, t177 + t129, (-pkin(4) * t137 + t231 * qJ(5)) * t222 + t122 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t133 * t126 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, -t254, 0, t102 * qJD(3), 0, -t166 + t225, -t248 - t167, -t211 + t225, -t13 * qJD(3), t248 - t206, -t31 * qJD(2) - t2 * qJD(3) - t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, -t249, -t193, 0, t249, -t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, t15, t129, t15, -t205, qJD(5) - t129, t134 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, 0, -t187, t172 - t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -qJ(5) * qJD(4) - t126 * qJD(3) - t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
