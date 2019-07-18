% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:24
% EndTime: 2019-07-18 13:27:25
% DurationCPUTime: 0.17s
% Computational Cost: add. (317->34), mult. (438->48), div. (0->0), fcn. (272->6), ass. (0->30)
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t36 = t31 * g(1) + t28 * g(3);
t24 = qJD(2) + qJD(3);
t23 = qJDD(2) + qJDD(3);
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t15 = qJDD(2) * pkin(1) + t36;
t34 = -t28 * g(1) + t31 * g(3);
t16 = -qJD(2) ^ 2 * pkin(1) - t34;
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t8 = t30 * t15 - t27 * t16;
t6 = t23 * pkin(2) + t8;
t22 = t24 ^ 2;
t9 = t27 * t15 + t30 * t16;
t7 = -t22 * pkin(2) + t9;
t3 = -t26 * t7 + t29 * t6;
t19 = qJD(4) + t24;
t17 = t19 ^ 2;
t18 = qJDD(4) + t23;
t33 = t26 * t17 - t29 * t18;
t35 = -pkin(2) * t33 + t3;
t4 = t26 * t6 + t29 * t7;
t12 = -t29 * t17 - t26 * t18;
t32 = pkin(2) * t12 - t4;
t25 = g(2) + qJDD(1);
t2 = t26 * t4 + t29 * t3;
t1 = pkin(2) * t2;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t36, t34, 0, 0, 0, 0, 0, 0, 0, t23, pkin(1) * (-t27 * t22 + t30 * t23) + t8, pkin(1) * (-t30 * t22 - t27 * t23) - t9, 0, pkin(1) * (t27 * t9 + t30 * t8), 0, 0, 0, 0, 0, t18, pkin(1) * (t12 * t27 - t30 * t33) + t35, pkin(1) * (t30 * t12 + t27 * t33) + t32, 0, t1 + pkin(1) * (t27 * (-t26 * t3 + t29 * t4) + t30 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t8, -t9, 0, 0, 0, 0, 0, 0, 0, t18, t35, t32, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t3, -t4, 0, 0;];
tauJ_reg  = t5;
