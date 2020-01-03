% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PPRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (327->65), mult. (591->103), div. (0->0), fcn. (444->8), ass. (0->45)
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t45 = qJD(3) ^ 2;
t26 = t40 * t45 * t42;
t22 = qJDD(4) + t26;
t50 = t40 * t22;
t23 = qJDD(4) - t26;
t49 = t42 * t23;
t48 = t40 * qJDD(3);
t47 = qJD(3) * qJD(4);
t36 = sin(pkin(6));
t38 = cos(pkin(6));
t16 = -t36 * g(1) + t38 * g(2) + qJDD(2);
t21 = -t38 * g(1) - t36 * g(2);
t33 = -g(3) + qJDD(1);
t35 = sin(pkin(7));
t37 = cos(pkin(7));
t12 = t37 * t21 + t35 * t33;
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t46 = -t35 * t21 + t37 * t33;
t8 = t43 * t12 + t41 * t46;
t6 = -t45 * pkin(3) + qJDD(3) * pkin(5) + t8;
t3 = -t42 * t16 + t40 * t6;
t4 = t40 * t16 + t42 * t6;
t1 = t40 * t3 + t42 * t4;
t7 = -t41 * t12 + t43 * t46;
t28 = t42 * qJDD(3);
t15 = -0.2e1 * t40 * t47 + t28;
t14 = 0.2e1 * t42 * t47 + t48;
t44 = qJD(4) ^ 2;
t32 = t42 ^ 2;
t31 = t40 ^ 2;
t30 = t32 * t45;
t29 = t31 * t45;
t25 = -t30 - t44;
t24 = -t29 - t44;
t20 = t29 + t30;
t19 = -t41 * qJDD(3) - t43 * t45;
t18 = t43 * qJDD(3) - t41 * t45;
t17 = (t31 + t32) * qJDD(3);
t10 = -t40 * t24 - t49;
t9 = t42 * t25 - t50;
t5 = -qJDD(3) * pkin(3) - t45 * pkin(5) - t7;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t12 + t37 * t46, 0, 0, 0, 0, 0, 0, t37 * t18 + t35 * t19, -t35 * t18 + t37 * t19, 0, t35 * (-t41 * t7 + t43 * t8) + t37 * (t41 * t8 + t43 * t7), 0, 0, 0, 0, 0, 0, t35 * (-t41 * t15 + t43 * t9) + t37 * (t43 * t15 + t41 * t9), t35 * (t43 * t10 + t41 * t14) + t37 * (t41 * t10 - t43 * t14), t35 * (t43 * t17 - t41 * t20) + t37 * (t41 * t17 + t43 * t20), t35 * (t43 * t1 + t41 * t5) + t37 * (t41 * t1 - t43 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, t42 * t22 + t40 * t25, -t40 * t23 + t42 * t24, 0, -t42 * t3 + t40 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t7, -t8, 0, 0, t14 * t40, t42 * t14 + t40 * t15, t50 + t42 * (-t29 + t44), t15 * t42, t40 * (t30 - t44) + t49, 0, pkin(3) * t15 + pkin(5) * t9 - t42 * t5, -pkin(3) * t14 + pkin(5) * t10 + t40 * t5, pkin(3) * t20 + pkin(5) * t17 + t1, -pkin(3) * t5 + pkin(5) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t29 - t30, t48, t26, t28, qJDD(4), -t3, -t4, 0, 0;];
tauJ_reg = t2;
