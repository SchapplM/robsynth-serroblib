% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:48:04
% EndTime: 2019-05-06 01:48:07
% DurationCPUTime: 0.79s
% Computational Cost: add. (741->105), mult. (1254->160), div. (0->0), fcn. (1419->6), ass. (0->76)
t48 = sin(qJ(4));
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t75 = cos(qJ(4));
t27 = t48 * t51 + t75 * t49;
t24 = t27 ^ 2;
t26 = t48 * t49 - t75 * t51;
t25 = t26 ^ 2;
t93 = t24 + t25;
t80 = t48 * pkin(3);
t39 = pkin(9) + t80;
t47 = sin(qJ(5));
t45 = t47 ^ 2;
t50 = cos(qJ(5));
t46 = t50 ^ 2;
t64 = t45 + t46;
t92 = t64 * t39;
t91 = t93 * t50;
t89 = t93 * t47;
t88 = -0.2e1 * t26;
t87 = 0.2e1 * t27;
t86 = -0.2e1 * t47;
t85 = -0.2e1 * t50;
t84 = 2 * qJ(2);
t83 = pkin(9) * t27;
t82 = t27 * pkin(5);
t81 = t47 * pkin(9);
t79 = t50 * pkin(9);
t52 = -pkin(1) - pkin(7);
t29 = (-pkin(8) + t52) * t49;
t41 = t51 * t52;
t30 = -t51 * pkin(8) + t41;
t14 = t48 * t29 - t75 * t30;
t32 = -pkin(5) * t47 + t50 * qJ(6);
t6 = t32 * t26 + t14;
t78 = t6 * t47;
t77 = t6 * t50;
t61 = t75 * pkin(3);
t40 = -t61 - pkin(4);
t76 = pkin(4) - t40;
t15 = t75 * t29 + t48 * t30;
t37 = t49 * pkin(3) + qJ(2);
t9 = t27 * pkin(4) + t26 * pkin(9) + t37;
t5 = t50 * t15 + t47 * t9;
t74 = t14 * t50;
t56 = t50 * pkin(5) + t47 * qJ(6);
t31 = -pkin(4) - t56;
t22 = -t61 + t31;
t73 = t26 * t22;
t72 = t26 * t31;
t21 = t26 * t47;
t20 = t26 * t50;
t71 = t27 * t39;
t70 = t47 * t39;
t69 = t47 * t50;
t19 = t50 * t27;
t68 = t50 * t39;
t67 = -t22 - t31;
t65 = pkin(9) * t64;
t63 = t27 * qJ(6);
t62 = t26 * t87;
t60 = t47 * t15 - t50 * t9;
t11 = t64 * t27;
t58 = pkin(4) * t26 - t83;
t2 = t63 + t5;
t3 = t60 - t82;
t1 = t2 * t50 + t3 * t47;
t57 = t72 + t83;
t55 = t71 + t73;
t54 = -t26 * t40 - t71;
t35 = 0.2e1 * t69;
t18 = t47 * t27;
t17 = t26 * t69;
t13 = t14 * t47;
t10 = (t45 - t46) * t26;
t4 = [1, 0, 0, -2 * pkin(1), t84, pkin(1) ^ 2 + qJ(2) ^ 2, t51 ^ 2, -0.2e1 * t51 * t49, 0, 0, 0, t49 * t84, t51 * t84, t25, t62, 0, 0, 0, t37 * t87, t37 * t88, t46 * t25, -0.2e1 * t25 * t69, t19 * t88, t47 * t62, t24, -0.2e1 * t14 * t21 - 0.2e1 * t27 * t60, -0.2e1 * t14 * t20 - 0.2e1 * t5 * t27, -0.2e1 * t6 * t21 - 0.2e1 * t3 * t27, 0.2e1 * (t2 * t47 - t3 * t50) * t26, 0.2e1 * t2 * t27 + 0.2e1 * t6 * t20, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t91, -t89, 0, t91, t1 * t27 + t6 * t26; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t24 + t25; 0, 0, 0, 0, 0, 0, 0, 0, t51, -t49, 0, t41, -t49 * t52, 0, 0, -t26, -t27, 0, -t14, -t15, -t17, t10, t18, t19, 0, t54 * t47 - t74, t54 * t50 + t13, -t55 * t47 - t77, t1, t55 * t50 - t78, t1 * t39 + t6 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t49, 0, 0, 0, 0, 0, -t26, -t27, 0, 0, 0, 0, 0, -t20, t21, -t20, t11, -t21, t27 * t92 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t61, -0.2e1 * t80, t45, t35, 0, 0, 0, t40 * t85, 0.2e1 * t40 * t47, t22 * t85, 0.2e1 * t92, t22 * t86, t64 * t39 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, -t14, -t15, -t17, t10, t18, t19, 0, t58 * t47 - t74, t58 * t50 + t13, -t57 * t47 - t77, t1, t57 * t50 - t78, t1 * pkin(9) + t6 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, 0, 0, 0, 0, -t20, t21, -t20, t11, -t21, pkin(9) * t11 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t61, -t80, t45, t35, 0, 0, 0, t76 * t50, -t76 * t47, t67 * t50, t65 + t92, t67 * t47, pkin(9) * t92 + t22 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, t35, 0, 0, 0, 0.2e1 * pkin(4) * t50, pkin(4) * t86, t31 * t85, 0.2e1 * t65, t31 * t86, t64 * pkin(9) ^ 2 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t21, t27, -t60, -t5, -t60 + 0.2e1 * t82, t56 * t26, 0.2e1 * t63 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, -t18, 0, t19, t32 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t50, 0, -t70, -t68, -t70, t32, t68, t32 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t50, 0, -t81, -t79, -t81, t32, t79, t32 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t20, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
