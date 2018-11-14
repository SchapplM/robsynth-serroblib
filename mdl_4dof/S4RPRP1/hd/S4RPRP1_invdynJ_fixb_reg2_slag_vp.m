% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRP1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:41
% EndTime: 2018-11-14 13:48:41
% DurationCPUTime: 0.22s
% Computational Cost: add. (344->78), mult. (563->91), div. (0->0), fcn. (304->10), ass. (0->57)
t48 = cos(pkin(6));
t36 = t48 * pkin(1) + pkin(2);
t47 = sin(pkin(6));
t76 = pkin(1) * t47;
t66 = qJD(3) * t76;
t79 = -qJD(1) * t66 + t36 * qJDD(1);
t45 = qJ(1) + pkin(6);
t41 = qJ(3) + t45;
t34 = sin(t41);
t35 = cos(t41);
t78 = g(1) * t34 - g(2) * t35;
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t13 = t49 * t36 + t51 * t76;
t43 = qJDD(1) + qJDD(3);
t37 = t43 * qJ(4);
t44 = qJD(1) + qJD(3);
t38 = t44 * qJD(4);
t77 = t37 + t38;
t69 = pkin(1) * qJDD(1);
t74 = t43 * pkin(3);
t22 = t36 * qJD(1);
t73 = t49 * t22;
t72 = t35 * pkin(3) + t34 * qJ(4);
t40 = cos(t45);
t52 = cos(qJ(1));
t71 = t52 * pkin(1) + pkin(2) * t40;
t67 = qJD(1) * t76;
t5 = t51 * t22 - t49 * t67;
t70 = qJD(4) - t5;
t68 = qJD(3) * t51;
t65 = t47 * t69;
t64 = -t34 * pkin(3) + t35 * qJ(4);
t63 = t22 * t68 + t79 * t49 + t51 * t65;
t62 = qJD(3) * t73 + t49 * t65 - t79 * t51;
t39 = sin(t45);
t50 = sin(qJ(1));
t60 = -t50 * pkin(1) - pkin(2) * t39;
t59 = g(1) * t50 - g(2) * t52;
t8 = t36 * t68 - t49 * t66;
t58 = g(1) * t35 + g(2) * t34 - t63;
t12 = t51 * t36 - t49 * t76;
t6 = t51 * t67 + t73;
t57 = -t62 + t78;
t2 = qJDD(4) + t62 - t74;
t56 = t5 * t44 + t58;
t55 = t6 * t44 + t57;
t9 = t13 * qJD(3);
t54 = -t9 * t44 + t57;
t46 = qJDD(2) - g(3);
t11 = -pkin(3) - t12;
t10 = qJ(4) + t13;
t7 = qJD(4) + t8;
t4 = t44 * qJ(4) + t6;
t3 = -t44 * pkin(3) + t70;
t1 = t63 + t77;
t14 = [0, 0, 0, 0, 0, qJDD(1), t59, g(1) * t52 + g(2) * t50, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t39 - g(2) * t40 + 0.2e1 * t48 * t69, g(1) * t40 + g(2) * t39 - 0.2e1 * t65, 0 (t59 + (t47 ^ 2 + t48 ^ 2) * t69) * pkin(1), 0, 0, 0, 0, 0, t43, t12 * t43 + t54, -t13 * t43 - t8 * t44 + t58, 0, -g(1) * t60 - g(2) * t71 - t62 * t12 + t63 * t13 - t5 * t9 + t6 * t8, 0, 0, 0, t43, 0, 0, -qJDD(4) + (pkin(3) - t11) * t43 + t54, 0, t10 * t43 + t7 * t44 - t58 + t77, t1 * t10 + t4 * t7 + t2 * t11 + t3 * t9 - g(1) * (t60 + t64) - g(2) * (t71 + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t55, t56, 0, 0, 0, 0, 0, t43, 0, 0, -qJDD(4) + t55 + 0.2e1 * t74, 0, 0.2e1 * t37 + 0.2e1 * t38 - t56, -t2 * pkin(3) - g(1) * t64 - g(2) * t72 + t1 * qJ(4) - t3 * t6 + t70 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, -t44 ^ 2, -t4 * t44 + t2 - t78;];
tau_reg  = t14;
