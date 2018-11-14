% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRP4
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:15
% EndTime: 2018-11-14 14:06:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->35), mult. (195->39), div. (0->0), fcn. (134->2), ass. (0->24)
t19 = sin(qJ(3));
t20 = cos(qJ(3));
t34 = -t19 * qJD(1) + t20 * qJD(2);
t33 = g(1) * t19 + g(2) * t20;
t12 = t20 * qJD(1) + t19 * qJD(2);
t32 = t12 * qJD(3);
t5 = qJD(3) * t34 + t20 * qJDD(1) + t19 * qJDD(2);
t31 = t5 * t20 - g(1);
t30 = t12 * t19;
t29 = qJD(3) * t19;
t28 = qJDD(3) * pkin(3);
t18 = qJDD(1) - g(1);
t26 = qJDD(2) + g(2);
t25 = t5 * t19 + t20 * t32 + g(2);
t24 = -t19 * qJDD(1) + t20 * qJDD(2);
t6 = t24 - t32;
t22 = t24 + t33;
t21 = qJD(3) ^ 2;
t10 = -t19 * qJDD(3) - t20 * t21;
t9 = -t20 * qJDD(3) + t19 * t21;
t8 = qJD(3) * pkin(3) + t34;
t4 = t6 + t28;
t1 = -t18 * t20 - t26 * t19;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, t10, t9, 0, -t6 * t19 + (-t20 * t34 - t30) * qJD(3) + t31, 0, 0, 0, 0, 0, 0, t10, t9, 0, -t4 * t19 + (-t20 * t8 - t30) * qJD(3) + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, -t9, t10, 0, t20 * t6 - t29 * t34 + t25, 0, 0, 0, 0, 0, 0, -t9, t10, 0, t20 * t4 - t8 * t29 + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t22, t1, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t22 + 0.2e1 * t28, t1, 0 (-t34 + t8) * t12 + (t4 + t33) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + g(3);];
tau_reg  = t2;
