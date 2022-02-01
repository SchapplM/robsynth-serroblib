% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:08
% DurationCPUTime: 0.36s
% Computational Cost: add. (560->80), mult. (1164->119), div. (0->0), fcn. (730->8), ass. (0->63)
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t79 = t58 ^ 2 + t60 ^ 2;
t57 = qJD(1) + qJD(3);
t49 = cos(pkin(8)) * pkin(1) + pkin(2);
t46 = t49 * qJD(1);
t63 = sin(qJ(3));
t87 = pkin(1) * sin(pkin(8));
t74 = qJD(3) * t87;
t72 = qJD(1) * t74;
t65 = cos(qJ(3));
t78 = qJD(3) * t65;
t70 = t46 * t78 - t63 * t72;
t15 = t57 * qJD(4) + t70;
t94 = t79 * t15;
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t38 = t64 * t58 + t62 * t60;
t29 = t38 * t57;
t93 = t57 * t79;
t82 = t63 * t46;
t20 = qJD(3) * t82 + t65 * t72;
t66 = t63 * t49 + t65 * t87;
t30 = t66 * qJD(3);
t92 = -t30 * t57 - t20;
t75 = qJD(1) * t87;
t25 = t65 * t75 + t82;
t91 = t25 * t57 - t20;
t24 = t65 * t46 - t63 * t75;
t90 = t24 - qJD(4);
t86 = t60 * pkin(4);
t50 = -pkin(3) - t86;
t11 = t50 * t57 - t90;
t81 = t64 * t60;
t83 = t62 * t58;
t37 = -t81 + t83;
t34 = t37 * qJD(5);
t89 = -t11 * t34 + t20 * t38;
t35 = t38 * qJD(5);
t88 = t11 * t35 + t20 * t37;
t77 = t57 * t83;
t76 = t57 * t81;
t71 = -t65 * t49 + t63 * t87 - pkin(3);
t68 = t79 * (t57 * qJ(4) + t25);
t67 = t49 * t78 - t63 * t74;
t53 = t60 * pkin(7);
t42 = t60 * qJ(4) + t53;
t41 = (-pkin(7) - qJ(4)) * t58;
t39 = qJD(5) * t76;
t33 = t35 * qJD(5);
t32 = t34 * qJD(5);
t31 = qJ(4) + t66;
t27 = -t76 + t77;
t26 = qJD(4) + t67;
t23 = t71 - t86;
t22 = t57 * t35;
t21 = -qJD(5) * t77 + t39;
t19 = t60 * t31 + t53;
t18 = (-pkin(7) - t31) * t58;
t16 = -t57 * pkin(3) - t90;
t2 = t21 * t38 - t29 * t34;
t1 = -t21 * t37 - t38 * t22 + t34 * t27 - t29 * t35;
t3 = [0, 0, 0, 0, 0, t92, -t67 * t57 - t70, t92 * t60, t26 * t93 + t94, t16 * t30 + t20 * t71 + t68 * t26 + t31 * t94, t2, t1, -t32, -t33, 0, t23 * t22 + t30 * t27 + ((-t18 * t62 - t19 * t64) * qJD(5) - t38 * t26) * qJD(5) + t88, t23 * t21 + t30 * t29 + ((-t18 * t64 + t19 * t62) * qJD(5) + t37 * t26) * qJD(5) + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32; 0, 0, 0, 0, 0, t91, t24 * t57 - t70, t91 * t60, -t90 * t93 + t94, -t20 * pkin(3) + qJ(4) * t94 - t16 * t25 - t90 * t68, t2, t1, -t32, -t33, 0, t50 * t22 - t25 * t27 + ((-t41 * t62 - t42 * t64) * qJD(5) + t90 * t38) * qJD(5) + t88, t50 * t21 - t25 * t29 + ((-t41 * t64 + t42 * t62) * qJD(5) - t90 * t37) * qJD(5) + t89; 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t57 ^ 2, -t68 * t57 + t20, 0, 0, 0, 0, 0, 0.2e1 * t29 * qJD(5), t39 + (-t27 - t77) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t27, -t27 ^ 2 + t29 ^ 2, t39 + (t27 - t77) * qJD(5), 0, 0, -t11 * t29 - t38 * t15, t11 * t27 + t37 * t15;];
tauc_reg = t3;
