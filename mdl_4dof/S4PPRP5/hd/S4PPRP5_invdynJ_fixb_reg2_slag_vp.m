% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRP5
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
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:07:16
% EndTime: 2018-11-14 14:07:16
% DurationCPUTime: 0.18s
% Computational Cost: add. (147->53), mult. (298->61), div. (0->0), fcn. (240->6), ass. (0->36)
t27 = sin(qJ(3));
t28 = cos(qJ(3));
t39 = qJD(1) * qJD(3);
t47 = qJDD(1) * t27 + t28 * t39;
t23 = pkin(5) + qJ(3);
t21 = sin(t23);
t22 = cos(t23);
t46 = g(1) * t22 - g(2) * t21;
t26 = cos(pkin(5));
t25 = sin(pkin(5));
t43 = t25 * t27;
t10 = -t28 * t26 + t43;
t6 = t10 * qJD(1);
t45 = qJD(4) + t6;
t42 = qJDD(3) * pkin(3);
t40 = qJDD(1) * t28;
t38 = qJDD(3) * qJ(4);
t37 = t25 * t40 + t47 * t26;
t36 = qJD(1) * t43;
t35 = t27 * t39;
t33 = (t35 - t40) * t26 + t47 * t25;
t32 = -t6 + t36;
t11 = t25 * t28 + t27 * t26;
t9 = t11 * qJD(3);
t31 = -qJD(3) * t9 - qJDD(3) * t10;
t7 = t11 * qJD(1);
t30 = g(1) * t21 + g(2) * t22 - t37;
t2 = qJDD(4) + t33 - t42;
t29 = qJD(3) * t7 - t33 - t46;
t24 = qJDD(2) + g(3);
t8 = t10 * qJD(3);
t5 = qJD(3) * qJ(4) + t7;
t4 = -qJD(3) * pkin(3) + t45;
t3 = -qJD(3) * t8 + qJDD(3) * t11;
t1 = t38 + (qJD(4) - t36) * qJD(3) + t37;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) + (t25 ^ 2 + t26 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t31, -t3, 0 (-t25 * t35 + t37) * t11 - t7 * t8 + t33 * t10 + t6 * t9 - g(1), 0, 0, 0, 0, 0, 0, t31, 0, t3, t1 * t11 + t10 * t2 + t4 * t9 - t5 * t8 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t29, qJD(3) * t32 + t30, 0, 0, 0, 0, 0, qJDD(3), 0, 0, -qJDD(4) + t29 + 0.2e1 * t42, 0, 0.2e1 * t38 + (0.2e1 * qJD(4) - t32) * qJD(3) - t30, t1 * qJ(4) - t2 * pkin(3) - t4 * t7 - g(1) * (pkin(3) * t22 + qJ(4) * t21) - g(2) * (-pkin(3) * t21 + qJ(4) * t22) + t45 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3), 0, -qJD(3) ^ 2, -qJD(3) * t5 + t2 + t46;];
tau_reg  = t12;
