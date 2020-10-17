% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:08:50
% EndTime: 2019-05-05 17:08:52
% DurationCPUTime: 0.53s
% Computational Cost: add. (638->88), mult. (1120->150), div. (0->0), fcn. (1331->8), ass. (0->58)
t46 = sin(pkin(9));
t48 = cos(pkin(9));
t51 = cos(qJ(3));
t71 = sin(qJ(3));
t26 = -t46 * t71 + t48 * t51;
t27 = -t46 * t51 - t48 * t71;
t78 = (t26 * t48 - t27 * t46) * pkin(3);
t25 = t27 ^ 2;
t76 = t26 ^ 2;
t63 = t25 + t76;
t52 = -pkin(1) - pkin(7);
t61 = t71 * t52;
t32 = -t71 * qJ(4) + t61;
t59 = (-qJ(4) + t52) * t51;
t18 = t46 * t32 - t48 * t59;
t77 = t18 ^ 2;
t45 = sin(pkin(10));
t47 = cos(pkin(10));
t49 = sin(qJ(6));
t50 = cos(qJ(6));
t28 = t49 * t45 - t50 * t47;
t69 = t26 * t28;
t75 = 0.2e1 * t69;
t40 = -t48 * pkin(3) - pkin(4);
t33 = -t47 * pkin(5) + t40;
t74 = 0.2e1 * t33;
t73 = 2 * qJ(2);
t35 = t46 * pkin(3) + qJ(5);
t72 = pkin(8) + t35;
t70 = t18 * t26;
t31 = t50 * t45 + t49 * t47;
t8 = t26 * t31;
t68 = t26 * t40;
t67 = t28 * t27;
t9 = t31 * t27;
t66 = t45 * t26;
t65 = t47 * t26;
t41 = t71 * pkin(3) + qJ(2);
t16 = -t27 * pkin(4) - t26 * qJ(5) + t41;
t20 = t48 * t32 + t46 * t59;
t6 = t45 * t16 + t47 * t20;
t64 = t45 ^ 2 + t47 ^ 2;
t60 = t64 * t27;
t5 = t47 * t16 - t45 * t20;
t58 = t6 * t45 + t5 * t47;
t57 = t5 * t45 - t6 * t47;
t56 = t20 * t27 + t70;
t55 = t27 * t35 + t68;
t23 = t72 * t47;
t22 = t72 * t45;
t15 = -t49 * t22 + t50 * t23;
t14 = -t50 * t22 - t49 * t23;
t7 = pkin(5) * t66 + t18;
t4 = -pkin(8) * t66 + t6;
t3 = -t27 * pkin(5) - pkin(8) * t65 + t5;
t2 = t49 * t3 + t50 * t4;
t1 = t50 * t3 - t49 * t4;
t10 = [1, 0, 0, -2 * pkin(1), t73, pkin(1) ^ 2 + qJ(2) ^ 2, t51 ^ 2, -0.2e1 * t51 * t71, 0, 0, 0, t71 * t73, t51 * t73, 0.2e1 * t56, t20 ^ 2 + t41 ^ 2 + t77, 0.2e1 * t18 * t66 - 0.2e1 * t5 * t27, 0.2e1 * t18 * t65 + 0.2e1 * t6 * t27, -0.2e1 * t58 * t26, t5 ^ 2 + t6 ^ 2 + t77, t69 ^ 2, t8 * t75, t27 * t75, 0.2e1 * t8 * t27, t25, -0.2e1 * t1 * t27 + 0.2e1 * t7 * t8, 0.2e1 * t2 * t27 - 0.2e1 * t69 * t7; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t63, -t56, -t63 * t45, -t63 * t47, 0, t57 * t27 - t70, 0, 0, 0, 0, 0, -t26 * t8 - t9 * t27, t26 * t69 + t27 * t67; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, t64 * t25 + t76, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t51, -t71, 0, t51 * t52, -t61, -t78 (-t18 * t48 + t20 * t46) * pkin(3), -t18 * t47 + t55 * t45, t18 * t45 + t55 * t47, -t57, t18 * t40 - t57 * t35, -t69 * t31, t28 * t69 - t31 * t8, -t9, t67, 0, -t14 * t27 + t7 * t28 + t33 * t8, t15 * t27 + t7 * t31 - t33 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t71, 0, t78, t65, -t66, -t60, -t35 * t60 - t68, 0, 0, 0, 0, 0, -t69, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t46 ^ 2 + t48 ^ 2) * pkin(3) ^ 2, -0.2e1 * t40 * t47, 0.2e1 * t40 * t45, 0.2e1 * t64 * t35, t64 * t35 ^ 2 + t40 ^ 2, t31 ^ 2, -0.2e1 * t31 * t28, 0, 0, 0, t28 * t74, t31 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t47 * t27, t45 * t27, -t64 * t26, t58, 0, 0, 0, 0, 0, t67, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, t18, 0, 0, 0, 0, 0, t8, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t45, 0, t40, 0, 0, 0, 0, 0, t28, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t8, -t27, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t28, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
