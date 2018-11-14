% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:32
% EndTime: 2018-11-14 13:52:32
% DurationCPUTime: 0.21s
% Computational Cost: add. (320->84), mult. (351->85), div. (0->0), fcn. (146->6), ass. (0->60)
t37 = qJ(1) + qJ(2);
t32 = sin(t37);
t33 = cos(t37);
t75 = -g(1) * t33 - g(2) * t32;
t65 = pkin(1) * qJDD(1);
t36 = qJD(1) + qJD(2);
t74 = t36 ^ 2;
t42 = -pkin(2) - pkin(3);
t38 = sin(qJ(2));
t66 = pkin(1) * qJD(1);
t60 = t38 * t66;
t10 = t36 * qJ(3) + t60;
t40 = cos(qJ(2));
t63 = qJD(2) * t40;
t12 = pkin(1) * t63 + qJD(3);
t18 = t38 * pkin(1) + qJ(3);
t35 = qJDD(1) + qJDD(2);
t28 = t35 * qJ(3);
t29 = t36 * qJD(3);
t57 = qJD(2) * t66;
t69 = t38 * t65 + t40 * t57;
t4 = t28 + t29 + t69;
t73 = t10 * t12 + t4 * t18;
t72 = t4 * qJ(3) + t10 * qJD(3);
t39 = sin(qJ(1));
t71 = t39 * pkin(1);
t70 = t10 * t40;
t68 = t33 * pkin(2) + t32 * qJ(3);
t67 = -t38 * t57 + t40 * t65;
t64 = qJD(2) * t38;
t41 = cos(qJ(1));
t62 = t41 * pkin(1) + t68;
t61 = pkin(1) * t64;
t59 = t40 * t66;
t58 = t36 * t64;
t25 = -t40 * pkin(1) - pkin(2);
t17 = t33 * qJ(3);
t56 = -t32 * pkin(2) + t17;
t55 = g(1) * t32 - g(2) * t33 + t67;
t54 = t69 + t75;
t11 = t36 * t60;
t31 = t35 * pkin(2);
t5 = -t31 + qJDD(3) - t67;
t53 = t42 * t32 + t17;
t52 = g(1) * t39 - g(2) * t41;
t51 = qJD(3) - t59;
t50 = -qJDD(3) + t55;
t49 = t4 + t75;
t48 = t31 + t50;
t47 = t36 * t59 - t54;
t46 = -t10 * t36 - t48;
t45 = -pkin(1) * t58 + t48;
t44 = t12 * t36 + t18 * t35 + t49;
t30 = t35 * pkin(3);
t19 = t33 * pkin(3);
t15 = -pkin(3) + t25;
t9 = -t36 * pkin(2) + t51;
t7 = t42 * t36 + t51;
t3 = -t30 + t5;
t1 = [0, 0, 0, 0, 0, qJDD(1), t52, g(1) * t41 + g(2) * t39, 0, 0, 0, 0, 0, 0, 0, t35 (t35 * t40 - t58) * pkin(1) + t55 (-t35 * t38 - t36 * t63) * pkin(1) - t54, 0 (t52 + (t38 ^ 2 + t40 ^ 2) * t65) * pkin(1), 0, 0, 0, t35, 0, 0, -t25 * t35 + t45, 0, t44, t5 * t25 + t9 * t61 - g(1) * (t56 - t71) - g(2) * t62 + t73, 0, 0, 0, 0, 0, t35, -t15 * t35 + t30 + t45, t44, 0, t3 * t15 + t7 * t61 - g(1) * (t53 - t71) - g(2) * (t19 + t62) + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t11 + t55, t47, 0, 0, 0, 0, 0, t35, 0, 0, 0.2e1 * t31 + t11 + t50, 0, 0.2e1 * t28 + 0.2e1 * t29 - t47, -t5 * pkin(2) - g(1) * t56 - g(2) * t68 + (-t38 * t9 - t70) * t66 + t72, 0, 0, 0, 0, 0, t35, -t42 * t35 + t11 + t30 + t48, t51 * t36 + t28 + t49, 0, t3 * t42 - g(1) * t53 - g(2) * (t19 + t68) + (-t38 * t7 - t70) * t66 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, -t74, t46, 0, 0, 0, 0, 0, 0, -t35, -t74, 0, -t30 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + g(3);];
tau_reg  = t1;
