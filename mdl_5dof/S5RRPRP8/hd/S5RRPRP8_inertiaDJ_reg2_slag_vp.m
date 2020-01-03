% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:27
% EndTime: 2019-12-31 20:04:30
% DurationCPUTime: 0.73s
% Computational Cost: add. (695->104), mult. (1562->175), div. (0->0), fcn. (1239->4), ass. (0->62)
t68 = cos(qJ(2));
t60 = t68 * qJD(2);
t66 = sin(qJ(4));
t67 = sin(qJ(2));
t83 = cos(qJ(4));
t72 = t67 * t83;
t89 = qJD(2) * t72 - t66 * t60;
t82 = t67 * t66;
t84 = pkin(6) - pkin(7);
t88 = t84 * t82;
t87 = -t68 * pkin(2) - t67 * qJ(3);
t86 = 2 * qJD(3);
t85 = -pkin(2) - pkin(3);
t47 = t84 * t68;
t21 = t83 * t47 + t88;
t81 = qJ(3) * t60 + t67 * qJD(3);
t61 = qJD(4) * t66;
t80 = t67 * qJD(2);
t79 = -0.2e1 * pkin(1) * qJD(2);
t45 = -pkin(1) + t87;
t62 = qJD(4) * t83;
t17 = -t68 * t61 + t67 * t62 - t89;
t57 = t83 * t85;
t27 = qJ(3) * t61 - t83 * qJD(3) - qJD(4) * t57;
t44 = t83 * qJ(3) + t66 * t85;
t28 = t66 * qJD(3) + t44 * qJD(4);
t37 = t68 * t83 + t82;
t38 = -t68 * t66 + t72;
t78 = -t44 * t17 + t27 * t37 + t28 * t38;
t77 = -t27 * t66 - t28 * t83 + t44 * t62;
t76 = pkin(6) * t80;
t75 = pkin(6) * t60;
t73 = t67 * t60;
t71 = t84 * t83;
t40 = t67 * t71;
t20 = -t66 * t47 + t40;
t35 = t68 * pkin(3) - t45;
t43 = -t66 * qJ(3) + t57;
t5 = -qJD(4) * t40 + t47 * t61 + t89 * t84;
t22 = t85 * t80 + t81;
t69 = t87 * qJD(2) + t68 * qJD(3);
t6 = t47 * t62 - t71 * t60 + (-qJD(2) + qJD(4)) * t88;
t49 = -0.2e1 * t73;
t48 = 0.2e1 * t73;
t46 = (-t67 ^ 2 + t68 ^ 2) * qJD(2);
t42 = -pkin(4) + t43;
t29 = pkin(2) * t80 - t81;
t26 = 0.2e1 * t28;
t25 = 0.2e1 * t27;
t19 = t37 * pkin(4) + t35;
t18 = t37 * qJD(2) - t67 * t61 - t68 * t62;
t16 = t44 * t27;
t12 = -t37 * qJ(5) + t21;
t11 = -t38 * qJ(5) + t20;
t9 = 0.2e1 * t38 * t18;
t8 = 0.2e1 * t37 * t17;
t7 = t17 * pkin(4) + t22;
t4 = -0.2e1 * t38 * t17 - 0.2e1 * t18 * t37;
t3 = -t83 * t18 - t66 * t17 + (-t83 * t37 + t38 * t66) * qJD(4);
t2 = t18 * qJ(5) + t38 * qJD(5) + t6;
t1 = t17 * qJ(5) + t37 * qJD(5) + t5;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0.2e1 * t46, 0, t49, 0, 0, t67 * t79, t68 * t79, 0, 0, t48, 0, -0.2e1 * t46, 0, 0, t49, -0.2e1 * t29 * t68 + 0.2e1 * t45 * t80, 0, -0.2e1 * t29 * t67 - 0.2e1 * t45 * t60, 0.2e1 * t45 * t29, t9, t4, 0, t8, 0, 0, 0.2e1 * t35 * t17 + 0.2e1 * t22 * t37, 0.2e1 * t35 * t18 + 0.2e1 * t22 * t38, -0.2e1 * t21 * t17 - 0.2e1 * t20 * t18 + 0.2e1 * t5 * t37 + 0.2e1 * t6 * t38, -0.2e1 * t20 * t6 - 0.2e1 * t21 * t5 + 0.2e1 * t35 * t22, t9, t4, 0, t8, 0, 0, 0.2e1 * t19 * t17 + 0.2e1 * t7 * t37, 0.2e1 * t19 * t18 + 0.2e1 * t7 * t38, 0.2e1 * t1 * t37 - 0.2e1 * t11 * t18 - 0.2e1 * t12 * t17 + 0.2e1 * t2 * t38, -0.2e1 * t12 * t1 - 0.2e1 * t11 * t2 + 0.2e1 * t19 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, -t80, 0, -t75, t76, 0, 0, 0, t60, 0, 0, t80, 0, -t75, t69, -t76, t69 * pkin(6), 0, 0, -t18, 0, t17, 0, t6, -t5, -t43 * t18 + t78, -t20 * t28 - t21 * t27 - t6 * t43 - t5 * t44, 0, 0, -t18, 0, t17, 0, t2, -t1, -t42 * t18 + t78, -t1 * t44 - t11 * t28 - t12 * t27 - t2 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, qJ(3) * t86, 0, 0, 0, 0, 0, 0, t26, -t25, 0, -0.2e1 * t43 * t28 - 0.2e1 * t16, 0, 0, 0, 0, 0, 0, t26, -t25, 0, -0.2e1 * t42 * t28 - 0.2e1 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t6 * t83 - t5 * t66 + (-t20 * t66 + t83 * t21) * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2 * t83 - t1 * t66 + (-t11 * t66 + t83 * t12) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t62, 0, -t43 * t61 + t77, 0, 0, 0, 0, 0, 0, t61, t62, 0, -t42 * t61 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t17, 0, -t6, t5, 0, 0, 0, 0, t18, 0, -t17, 0, -t2, t1, -t18 * pkin(4), -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t27, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t27, 0, -t28 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t62, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t62, 0, -pkin(4) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
