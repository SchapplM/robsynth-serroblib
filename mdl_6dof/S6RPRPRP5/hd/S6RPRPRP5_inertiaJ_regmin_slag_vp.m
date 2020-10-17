% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:48:52
% EndTime: 2019-05-05 17:48:54
% DurationCPUTime: 0.70s
% Computational Cost: add. (1161->103), mult. (2239->189), div. (0->0), fcn. (2660->8), ass. (0->60)
t48 = sin(pkin(10));
t50 = cos(pkin(10));
t52 = sin(qJ(5));
t74 = cos(qJ(5));
t81 = -t52 * t48 + t74 * t50;
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t53 = sin(qJ(3));
t75 = cos(qJ(3));
t34 = t53 * t49 - t75 * t51;
t80 = -0.2e1 * t34;
t35 = t74 * t48 + t52 * t50;
t79 = -0.2e1 * t35;
t42 = -t50 * pkin(4) - pkin(3);
t78 = 0.2e1 * t42;
t43 = -t51 * pkin(2) - pkin(1);
t77 = 0.2e1 * t43;
t36 = t75 * t49 + t53 * t51;
t20 = t34 * pkin(3) - t36 * qJ(4) + t43;
t67 = pkin(7) + qJ(2);
t37 = t67 * t49;
t39 = t67 * t51;
t25 = -t53 * t37 + t75 * t39;
t12 = t48 * t20 + t50 * t25;
t70 = t48 * t36;
t10 = -pkin(8) * t70 + t12;
t11 = t50 * t20 - t48 * t25;
t69 = t50 * t36;
t8 = t34 * pkin(4) - pkin(8) * t69 + t11;
t4 = t74 * t10 + t52 * t8;
t76 = t34 * pkin(5);
t66 = pkin(8) + qJ(4);
t38 = t66 * t50;
t60 = t66 * t48;
t22 = t52 * t38 + t74 * t60;
t73 = t22 * t34;
t24 = t74 * t38 - t52 * t60;
t72 = t24 * t34;
t16 = t81 * t36;
t71 = t81 * t16;
t26 = t81 * t34;
t27 = t35 * t34;
t65 = t48 ^ 2 + t50 ^ 2;
t64 = t49 ^ 2 + t51 ^ 2;
t63 = t34 * qJ(6);
t61 = t52 * t10 - t74 * t8;
t23 = t75 * t37 + t53 * t39;
t59 = -pkin(3) * t36 - qJ(4) * t34;
t58 = pkin(5) * t81 + t35 * qJ(6);
t57 = t11 * t50 + t12 * t48;
t56 = -t11 * t48 + t12 * t50;
t14 = pkin(4) * t70 + t23;
t30 = t35 ^ 2;
t19 = t42 - t58;
t15 = t35 * t36;
t13 = t35 * t15;
t5 = t15 * pkin(5) - t16 * qJ(6) + t14;
t2 = t61 - t76;
t1 = t63 + t4;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t51, -0.2e1 * pkin(1) * t49, 0.2e1 * t64 * qJ(2), t64 * qJ(2) ^ 2 + pkin(1) ^ 2, t36 ^ 2, t36 * t80, 0, 0, 0, t34 * t77, t36 * t77, 0.2e1 * t11 * t34 + 0.2e1 * t23 * t70, -0.2e1 * t12 * t34 + 0.2e1 * t23 * t69, -0.2e1 * t57 * t36, t11 ^ 2 + t12 ^ 2 + t23 ^ 2, t16 ^ 2, -0.2e1 * t16 * t15, 0.2e1 * t16 * t34, t15 * t80, t34 ^ 2, 0.2e1 * t14 * t15 - 0.2e1 * t34 * t61, 0.2e1 * t14 * t16 - 0.2e1 * t4 * t34, 0.2e1 * t5 * t15 - 0.2e1 * t2 * t34, -0.2e1 * t1 * t15 + 0.2e1 * t2 * t16, 0.2e1 * t1 * t34 - 0.2e1 * t5 * t16, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t51, t49, 0, -pkin(1), 0, 0, 0, 0, 0, t34, t36, t50 * t34, -t48 * t34, -t65 * t36, t57, 0, 0, 0, 0, 0, t26, -t27, t26, -t13 - t71, t27, t1 * t35 - t2 * t81; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 ^ 2 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, 0, -t23, -t25, -t23 * t50 + t59 * t48, t23 * t48 + t59 * t50, t56, -t23 * pkin(3) + t56 * qJ(4), t16 * t35, -t13 + t71, t27, t26, 0, -t14 * t81 + t42 * t15 - t73, t14 * t35 + t42 * t16 - t72, t19 * t15 - t5 * t81 - t73, t1 * t81 - t24 * t15 + t22 * t16 + t2 * t35, -t19 * t16 - t5 * t35 + t72, t1 * t24 + t5 * t19 + t2 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t81 + t35 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t50, -0.2e1 * pkin(3) * t48, 0.2e1 * t65 * qJ(4), t65 * qJ(4) ^ 2 + pkin(3) ^ 2, t30, -t81 * t79, 0, 0, 0, -t81 * t78, t35 * t78, -0.2e1 * t19 * t81, 0.2e1 * t22 * t35 + 0.2e1 * t24 * t81, t19 * t79, t19 ^ 2 + t22 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t69, 0, t23, 0, 0, 0, 0, 0, t15, t16, t15, 0, -t16, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t48, 0, -pkin(3), 0, 0, 0, 0, 0, -t81, t35, -t81, 0, -t35, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t34, -t61, -t4, -t61 + 0.2e1 * t76, -pkin(5) * t16 - t15 * qJ(6), 0.2e1 * t63 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t35, t81, 0, t35, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t81, 0, -t22, -t24, -t22, -pkin(5) * t35 + qJ(6) * t81, t24, -t22 * pkin(5) + t24 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
