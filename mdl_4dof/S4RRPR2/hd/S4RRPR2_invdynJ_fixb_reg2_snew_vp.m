% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR2
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:35
% DurationCPUTime: 0.21s
% Computational Cost: add. (918->54), mult. (1066->64), div. (0->0), fcn. (514->6), ass. (0->43)
t68 = pkin(2) + pkin(3);
t46 = (qJD(1) + qJD(2));
t44 = t46 ^ 2;
t45 = qJDD(1) + qJDD(2);
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t67 = pkin(1) * (t51 * t44 + t48 * t45);
t40 = t45 * pkin(2);
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t64 = t49 * g(1) - t52 * g(2);
t30 = qJDD(1) * pkin(1) + t64;
t59 = t52 * g(1) + t49 * g(2);
t31 = -qJD(1) ^ 2 * pkin(1) - t59;
t16 = t51 * t30 - t48 * t31;
t58 = -qJDD(3) + t16;
t15 = -t44 * qJ(3) - t40 - t58;
t12 = -t45 * pkin(3) + t15;
t47 = sin(qJ(4));
t50 = cos(qJ(4));
t37 = t45 * qJ(3);
t17 = t48 * t30 + t51 * t31;
t65 = (2 * qJD(3) * t46) + t17;
t62 = t37 + t65;
t9 = -t68 * t44 + t62;
t6 = t47 * t12 + t50 * t9;
t13 = -t44 * pkin(2) + t62;
t66 = -pkin(2) * t15 + qJ(3) * t13;
t5 = t50 * t12 - t47 * t9;
t42 = qJD(4) - t46;
t38 = t42 ^ 2;
t39 = -qJDD(4) + t45;
t63 = -t50 * t38 + t47 * t39;
t61 = 0.2e1 * t37 + t65;
t3 = t47 * t6 + t50 * t5;
t4 = -t47 * t5 + t50 * t6;
t60 = qJ(3) * t4 - t68 * t3;
t57 = t47 * t38 + t50 * t39;
t55 = 0.2e1 * t40 + t58;
t54 = qJ(3) * t57 - t68 * t63 + t6;
t53 = qJ(3) * t63 + t68 * t57 - t5;
t26 = pkin(1) * (-t48 * t44 + t51 * t45);
t1 = [0, 0, 0, 0, 0, qJDD(1), t64, t59, 0, 0, 0, 0, 0, 0, 0, t45, t16 + t26, -t17 - t67, 0, pkin(1) * (t51 * t16 + t48 * t17), 0, 0, 0, t45, 0, 0, t26 + t55, 0, t61 + t67, pkin(1) * (t48 * t13 - t51 * t15) + t66, 0, 0, 0, 0, 0, t39, pkin(1) * (t48 * t63 + t51 * t57) + t53, pkin(1) * (t48 * t57 - t51 * t63) + t54, 0, pkin(1) * (-t51 * t3 + t48 * t4) + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t16, -t17, 0, 0, 0, 0, 0, t45, 0, 0, t55, 0, t61, t66, 0, 0, 0, 0, 0, t39, t53, t54, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, -t44, t15, 0, 0, 0, 0, 0, 0, -t57, t63, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t5, -t6, 0, 0;];
tauJ_reg  = t1;
