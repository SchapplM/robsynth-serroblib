% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PPRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:52
% EndTime: 2019-12-31 16:19:53
% DurationCPUTime: 0.18s
% Computational Cost: add. (208->53), mult. (385->77), div. (0->0), fcn. (252->6), ass. (0->37)
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t38 = qJD(3) ^ 2;
t21 = t33 * t38 * t35;
t44 = t33 * (qJDD(4) + t21);
t43 = t35 * (qJDD(4) - t21);
t42 = t33 * qJDD(3);
t41 = qJD(3) * qJD(4);
t30 = sin(pkin(6));
t31 = cos(pkin(6));
t40 = -t31 * g(1) - t30 * g(2);
t28 = -g(3) + qJDD(1);
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t39 = -t30 * g(1) + t31 * g(2) + qJDD(2);
t10 = t36 * t28 + t34 * t39;
t6 = -t38 * pkin(3) + qJDD(3) * pkin(5) + t10;
t3 = t33 * t6 - t35 * t40;
t4 = t33 * t40 + t35 * t6;
t1 = t33 * t3 + t35 * t4;
t9 = -t34 * t28 + t36 * t39;
t23 = t35 * qJDD(3);
t14 = -0.2e1 * t33 * t41 + t23;
t13 = 0.2e1 * t35 * t41 + t42;
t37 = qJD(4) ^ 2;
t27 = t35 ^ 2;
t26 = t33 ^ 2;
t25 = t27 * t38;
t24 = t26 * t38;
t18 = t24 + t25;
t17 = -t34 * qJDD(3) - t36 * t38;
t16 = t36 * qJDD(3) - t34 * t38;
t15 = (t26 + t27) * qJDD(3);
t8 = -t43 - t33 * (-t24 - t37);
t7 = t35 * (-t25 - t37) - t44;
t5 = -qJDD(3) * pkin(3) - t38 * pkin(5) - t9;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t36 * t10 - t34 * t9, 0, 0, 0, 0, 0, 0, -t34 * t14 + t36 * t7, t34 * t13 + t36 * t8, t36 * t15 - t34 * t18, t36 * t1 + t34 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, t16, t17, 0, t34 * t10 + t36 * t9, 0, 0, 0, 0, 0, 0, t36 * t14 + t34 * t7, -t36 * t13 + t34 * t8, t34 * t15 + t36 * t18, t34 * t1 - t36 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t9, -t10, 0, 0, t13 * t33, t35 * t13 + t33 * t14, t44 + t35 * (-t24 + t37), t14 * t35, t33 * (t25 - t37) + t43, 0, pkin(3) * t14 + pkin(5) * t7 - t35 * t5, -pkin(3) * t13 + pkin(5) * t8 + t33 * t5, pkin(3) * t18 + pkin(5) * t15 + t1, -pkin(3) * t5 + pkin(5) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t24 - t25, t42, t21, t23, qJDD(4), -t3, -t4, 0, 0;];
tauJ_reg = t2;
