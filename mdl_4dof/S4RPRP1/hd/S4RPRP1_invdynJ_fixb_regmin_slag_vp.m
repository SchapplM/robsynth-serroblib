% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:43
% EndTime: 2019-03-08 18:29:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (293->69), mult. (447->78), div. (0->0), fcn. (238->10), ass. (0->52)
t41 = cos(pkin(6));
t32 = t41 * pkin(1) + pkin(2);
t40 = sin(pkin(6));
t69 = pkin(1) * t40;
t60 = qJD(3) * t69;
t72 = -qJD(1) * t60 + t32 * qJDD(1);
t38 = qJ(1) + pkin(6);
t35 = qJ(3) + t38;
t30 = sin(t35);
t31 = cos(t35);
t71 = g(1) * t30 - g(2) * t31;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t64 = t42 * t32 + t44 * t69;
t36 = qJDD(1) + qJDD(3);
t33 = t36 * qJ(4);
t37 = qJD(1) + qJD(3);
t34 = t37 * qJD(4);
t70 = t33 + t34;
t67 = t36 * pkin(3);
t19 = t32 * qJD(1);
t66 = t42 * t19;
t65 = t31 * pkin(3) + t30 * qJ(4);
t61 = qJD(1) * t69;
t5 = t44 * t19 - t42 * t61;
t63 = qJD(4) - t5;
t62 = qJD(3) * t44;
t59 = qJDD(1) * t69;
t58 = -t30 * pkin(3) + t31 * qJ(4);
t57 = t19 * t62 + t72 * t42 + t44 * t59;
t56 = -qJD(3) * t66 - t42 * t59 + t72 * t44;
t43 = sin(qJ(1));
t45 = cos(qJ(1));
t54 = g(1) * t43 - g(2) * t45;
t53 = t32 * t62 - t42 * t60;
t52 = g(1) * t31 + g(2) * t30 - t57;
t51 = t44 * t32 - t42 * t69;
t6 = t44 * t61 + t66;
t50 = t56 + t71;
t2 = qJDD(4) - t56 - t67;
t49 = t5 * t37 + t52;
t48 = t6 * t37 + t50;
t8 = t64 * qJD(3);
t47 = -t8 * t37 + t50;
t39 = qJDD(2) - g(3);
t10 = -pkin(3) - t51;
t9 = qJ(4) + t64;
t7 = qJD(4) + t53;
t4 = t37 * qJ(4) + t6;
t3 = -t37 * pkin(3) + t63;
t1 = t57 + t70;
t11 = [qJDD(1), t54, g(1) * t45 + g(2) * t43 (t54 + (t40 ^ 2 + t41 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t36, t36 * t51 + t47, -t64 * t36 - t53 * t37 + t52, -qJDD(4) + (pkin(3) - t10) * t36 + t47, t9 * t36 + t7 * t37 - t52 + t70, t1 * t9 + t4 * t7 + t2 * t10 + t3 * t8 - g(1) * (-pkin(2) * sin(t38) - t43 * pkin(1) + t58) - g(2) * (pkin(2) * cos(t38) + t45 * pkin(1) + t65); 0, 0, 0, t39, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, t36, t48, t49, -qJDD(4) + t48 + 0.2e1 * t67, 0.2e1 * t33 + 0.2e1 * t34 - t49, -t2 * pkin(3) - g(1) * t58 - g(2) * t65 + t1 * qJ(4) - t3 * t6 + t63 * t4; 0, 0, 0, 0, 0, 0, 0, -t36, -t37 ^ 2, -t4 * t37 + t2 - t71;];
tau_reg  = t11;
