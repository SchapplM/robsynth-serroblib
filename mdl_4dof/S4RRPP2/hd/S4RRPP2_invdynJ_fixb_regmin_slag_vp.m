% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x12]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:31
% EndTime: 2018-11-14 13:52:31
% DurationCPUTime: 0.19s
% Computational Cost: add. (313->82), mult. (339->81), div. (0->0), fcn. (142->6), ass. (0->59)
t37 = qJ(1) + qJ(2);
t32 = sin(t37);
t33 = cos(t37);
t73 = -g(1) * t33 - g(2) * t32;
t36 = qJD(1) + qJD(2);
t72 = t36 ^ 2;
t42 = -pkin(2) - pkin(3);
t38 = sin(qJ(2));
t64 = pkin(1) * qJD(1);
t58 = t38 * t64;
t10 = t36 * qJ(3) + t58;
t40 = cos(qJ(2));
t61 = qJD(2) * t40;
t12 = pkin(1) * t61 + qJD(3);
t18 = t38 * pkin(1) + qJ(3);
t35 = qJDD(1) + qJDD(2);
t28 = t35 * qJ(3);
t29 = t36 * qJD(3);
t55 = qJD(2) * t64;
t63 = pkin(1) * qJDD(1);
t67 = t38 * t63 + t40 * t55;
t4 = t28 + t29 + t67;
t71 = t10 * t12 + t4 * t18;
t70 = t4 * qJ(3) + t10 * qJD(3);
t39 = sin(qJ(1));
t69 = t39 * pkin(1);
t68 = t10 * t40;
t66 = t33 * pkin(2) + t32 * qJ(3);
t65 = -t38 * t55 + t40 * t63;
t62 = qJD(2) * t38;
t41 = cos(qJ(1));
t60 = t41 * pkin(1) + t66;
t59 = pkin(1) * t62;
t57 = t40 * t64;
t56 = t36 * t62;
t25 = -t40 * pkin(1) - pkin(2);
t17 = t33 * qJ(3);
t54 = -t32 * pkin(2) + t17;
t53 = g(1) * t32 - g(2) * t33 + t65;
t52 = t67 + t73;
t11 = t36 * t58;
t31 = t35 * pkin(2);
t5 = -t31 + qJDD(3) - t65;
t51 = t42 * t32 + t17;
t50 = qJD(3) - t57;
t49 = -qJDD(3) + t53;
t48 = t4 + t73;
t47 = t31 + t49;
t46 = t36 * t57 - t52;
t45 = -t10 * t36 - t47;
t44 = -pkin(1) * t56 + t47;
t43 = t12 * t36 + t18 * t35 + t48;
t30 = t35 * pkin(3);
t19 = t33 * pkin(3);
t15 = -pkin(3) + t25;
t9 = -t36 * pkin(2) + t50;
t7 = t42 * t36 + t50;
t3 = -t30 + t5;
t1 = [qJDD(1), g(1) * t39 - g(2) * t41, g(1) * t41 + g(2) * t39, t35 (t35 * t40 - t56) * pkin(1) + t53 (-t35 * t38 - t36 * t61) * pkin(1) - t52, -t25 * t35 + t44, t43, t5 * t25 + t9 * t59 - g(1) * (t54 - t69) - g(2) * t60 + t71, -t15 * t35 + t30 + t44, t43, t3 * t15 + t7 * t59 - g(1) * (t51 - t69) - g(2) * (t19 + t60) + t71; 0, 0, 0, t35, t11 + t53, t46, 0.2e1 * t31 + t11 + t49, 0.2e1 * t28 + 0.2e1 * t29 - t46, -t5 * pkin(2) - g(1) * t54 - g(2) * t66 + (-t38 * t9 - t68) * t64 + t70, -t42 * t35 + t11 + t30 + t47, t50 * t36 + t28 + t48, t3 * t42 - g(1) * t51 - g(2) * (t19 + t66) + (-t38 * t7 - t68) * t64 + t70; 0, 0, 0, 0, 0, 0, -t35, -t72, t45, -t35, -t72, -t30 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + g(3);];
tau_reg  = t1;
