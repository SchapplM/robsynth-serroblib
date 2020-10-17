% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:13:29
% EndTime: 2019-05-05 19:13:31
% DurationCPUTime: 0.53s
% Computational Cost: add. (590->70), mult. (1038->121), div. (0->0), fcn. (1309->8), ass. (0->58)
t46 = sin(pkin(10));
t47 = cos(pkin(10));
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t31 = -t46 * t50 + t47 * t52;
t32 = -t46 * t52 - t47 * t50;
t49 = sin(qJ(5));
t68 = cos(qJ(5));
t17 = t68 * t31 + t49 * t32;
t14 = t17 ^ 2;
t48 = sin(qJ(6));
t77 = t48 * t17;
t51 = cos(qJ(6));
t76 = t51 * t17;
t55 = t49 * t31 - t68 * t32;
t75 = t55 ^ 2;
t74 = (t31 * t47 - t32 * t46) * pkin(3);
t11 = t48 * t55;
t12 = t51 * t55;
t41 = t50 * pkin(3) + qJ(2);
t21 = -t32 * pkin(4) + t41;
t73 = 0.2e1 * t21;
t72 = 0.2e1 * qJ(2);
t71 = pkin(3) * t46;
t53 = -pkin(1) - pkin(7);
t34 = (-qJ(4) + t53) * t50;
t42 = t52 * t53;
t35 = -t52 * qJ(4) + t42;
t19 = -t46 * t34 + t47 * t35;
t56 = -t31 * pkin(8) + t19;
t20 = t47 * t34 + t46 * t35;
t9 = t32 * pkin(8) + t20;
t4 = t49 * t9 - t68 * t56;
t70 = t4 * t51;
t40 = t47 * pkin(3) + pkin(4);
t24 = t68 * t40 - t49 * t71;
t22 = -pkin(5) - t24;
t69 = pkin(5) - t22;
t65 = t48 * t51;
t63 = -0.2e1 * t17 * t55;
t62 = t31 ^ 2 + t32 ^ 2;
t61 = -pkin(5) * t17 - pkin(9) * t55;
t60 = -t75 - t14;
t25 = -t49 * t40 - t68 * t71;
t23 = pkin(9) - t25;
t59 = t17 * t22 - t23 * t55;
t58 = t19 * t31 - t20 * t32;
t45 = t51 ^ 2;
t44 = t48 ^ 2;
t36 = 0.2e1 * t65;
t10 = t48 * t76;
t7 = (-t44 + t45) * t17;
t6 = pkin(5) * t55 - t17 * pkin(9) + t21;
t5 = t49 * t56 + t68 * t9;
t3 = t4 * t48;
t2 = t48 * t6 + t51 * t5;
t1 = -t48 * t5 + t51 * t6;
t8 = [1, 0, 0, -2 * pkin(1), t72 (pkin(1) ^ 2) + qJ(2) ^ 2, t52 ^ 2, -0.2e1 * t52 * t50, 0, 0, 0, t50 * t72, t52 * t72, -0.2e1 * t58, t19 ^ 2 + t20 ^ 2 + t41 ^ 2, t14, t63, 0, 0, 0, t55 * t73, t17 * t73, t45 * t14, -0.2e1 * t14 * t65, 0.2e1 * t55 * t76, t48 * t63, t75, 0.2e1 * t1 * t55 + 0.2e1 * t4 * t77, -0.2e1 * t2 * t55 + 0.2e1 * t4 * t76; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t62, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t48, t60 * t51; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t52, -t50, 0, t42, -t50 * t53, -t74 (t19 * t47 + t20 * t46) * pkin(3), 0, 0, t17, -t55, 0, -t4, -t5, t10, t7, t11, t12, 0, t59 * t48 - t70, t59 * t51 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t50, 0, t74, 0, 0, 0, 0, 0, t17, -t55, 0, 0, 0, 0, 0, t76, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t46 ^ 2 + t47 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t24, 0.2e1 * t25, t44, t36, 0, 0, 0, -0.2e1 * t22 * t51, 0.2e1 * t22 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, t55, t17, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t55, 0, -t4, -t5, t10, t7, t11, t12, 0, t61 * t48 - t70, t61 * t51 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t55, 0, 0, 0, 0, 0, t76, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, t25, t44, t36, 0, 0, 0, t69 * t51, -t69 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t44, t36, 0, 0, 0, 0.2e1 * pkin(5) * t51, -0.2e1 * pkin(5) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t77, t55, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t48 * t23, -t51 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t48 * pkin(9), -t51 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
