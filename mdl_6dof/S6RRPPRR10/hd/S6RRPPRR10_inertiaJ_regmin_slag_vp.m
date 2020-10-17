% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:49:39
% EndTime: 2019-05-06 11:49:42
% DurationCPUTime: 0.74s
% Computational Cost: add. (704->100), mult. (1350->179), div. (0->0), fcn. (1538->8), ass. (0->74)
t66 = sin(qJ(2));
t61 = sin(pkin(10));
t62 = cos(pkin(10));
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t39 = t68 * t61 + t65 * t62;
t40 = -t65 * t61 + t68 * t62;
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t73 = t67 * t39 + t64 * t40;
t96 = t73 * t66;
t75 = -t64 * t39 + t67 * t40;
t95 = t75 * t66;
t51 = t61 * pkin(4) + qJ(3);
t28 = t39 * pkin(5) + t51;
t94 = 0.2e1 * t28;
t93 = 0.2e1 * t51;
t92 = -0.2e1 * t66;
t91 = 0.2e1 * t66;
t69 = cos(qJ(2));
t90 = 0.2e1 * t69;
t89 = 0.2e1 * qJ(3);
t88 = t64 * pkin(5);
t87 = t66 * pkin(5);
t86 = t67 * pkin(5);
t29 = t40 * t69;
t63 = -pkin(2) - qJ(4);
t76 = -t66 * qJ(3) - pkin(1);
t37 = t63 * t69 + t76;
t53 = t66 * pkin(7);
t46 = t66 * pkin(3) + t53;
t42 = t62 * t46;
t16 = t66 * pkin(4) + t42 + (pkin(8) * t69 - t37) * t61;
t24 = t62 * t37 + t61 * t46;
t81 = t62 * t69;
t17 = -pkin(8) * t81 + t24;
t9 = t65 * t16 + t68 * t17;
t5 = -t29 * pkin(9) + t9;
t85 = t67 * t5;
t84 = -pkin(8) + t63;
t83 = t39 * t66;
t82 = t61 * t69;
t54 = t69 * pkin(7);
t47 = t69 * pkin(3) + t54;
t49 = t61 ^ 2 + t62 ^ 2;
t59 = t66 ^ 2;
t78 = t69 ^ 2 + t59;
t77 = t69 * qJ(3);
t35 = pkin(4) * t81 + t47;
t30 = t39 * t69;
t8 = t68 * t16 - t65 * t17;
t4 = t30 * pkin(9) + t8 + t87;
t1 = t67 * t4 - t64 * t5;
t43 = t84 * t61;
t44 = t84 * t62;
t25 = -t65 * t43 + t68 * t44;
t74 = -t66 * pkin(2) + t77;
t23 = -t61 * t37 + t42;
t10 = t23 * t62 + t24 * t61;
t26 = t68 * t43 + t65 * t44;
t72 = t63 * t66 + t77;
t70 = qJ(3) ^ 2;
t45 = -t69 * pkin(2) + t76;
t38 = t49 * t63;
t34 = t40 * t66;
t18 = t29 * pkin(5) + t35;
t15 = -t64 * t29 - t67 * t30;
t14 = t67 * t29 - t64 * t30;
t13 = -t39 * pkin(9) + t26;
t12 = -t40 * pkin(9) + t25;
t7 = t64 * t12 + t67 * t13;
t6 = t67 * t12 - t64 * t13;
t2 = t64 * t4 + t85;
t3 = [1, 0, 0, t59, t66 * t90, 0, 0, 0, pkin(1) * t90, pkin(1) * t92, 0.2e1 * t78 * pkin(7), t45 * t90, t45 * t92, t78 * pkin(7) ^ 2 + t45 ^ 2, 0.2e1 * t23 * t66 + 0.2e1 * t47 * t81, -0.2e1 * t24 * t66 - 0.2e1 * t47 * t82 (t23 * t61 - t24 * t62) * t90, t23 ^ 2 + t24 ^ 2 + t47 ^ 2, t30 ^ 2, 0.2e1 * t30 * t29, -t30 * t91, -t29 * t91, t59, 0.2e1 * t35 * t29 + 0.2e1 * t8 * t66, -0.2e1 * t35 * t30 - 0.2e1 * t9 * t66, t15 ^ 2, -0.2e1 * t15 * t14, t15 * t91, t14 * t92, t59, 0.2e1 * t1 * t66 + 0.2e1 * t18 * t14, 0.2e1 * t18 * t15 - 0.2e1 * t2 * t66; 0, 0, 0, 0, 0, t66, t69, 0, -t53, -t54, t74, t53, t54, t74 * pkin(7), t47 * t61 + t72 * t62, t47 * t62 - t72 * t61, -t10, t47 * qJ(3) + t10 * t63, -t30 * t40, -t40 * t29 + t30 * t39, t34, -t83, 0, t25 * t66 + t51 * t29 + t35 * t39, -t26 * t66 - t51 * t30 + t35 * t40, t15 * t75, -t14 * t75 - t15 * t73, t95, -t96, 0, t28 * t14 + t18 * t73 + t6 * t66, t28 * t15 + t18 * t75 - t7 * t66; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t89, pkin(2) ^ 2 + t70, t61 * t89, t62 * t89, -0.2e1 * t38, t49 * t63 ^ 2 + t70, t40 ^ 2, -0.2e1 * t40 * t39, 0, 0, 0, t39 * t93, t40 * t93, t75 ^ 2, -0.2e1 * t75 * t73, 0, 0, 0, t73 * t94, t75 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, t53, t62 * t66, -t61 * t66, 0, t10, 0, 0, 0, 0, 0, t34, -t83, 0, 0, 0, 0, 0, t95, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, -t49, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t82, 0, t47, 0, 0, 0, 0, 0, t29, -t30, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t62, 0, qJ(3), 0, 0, 0, 0, 0, t39, t40, 0, 0, 0, 0, 0, t73, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, t66, t8, -t9, 0, 0, t15, -t14, t66, t66 * t86 + t1, -t85 + (-t4 - t87) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t39, 0, t25, -t26, 0, 0, t75, -t73, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t39, 0, 0, 0, 0, 0, t75, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t86, -0.2e1 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t66, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t73, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t86, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
