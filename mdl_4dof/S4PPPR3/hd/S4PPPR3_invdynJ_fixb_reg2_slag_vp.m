% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:57
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:56:34
% EndTime: 2018-11-14 13:56:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (124->35), mult. (228->44), div. (0->0), fcn. (212->6), ass. (0->26)
t23 = sin(pkin(5));
t24 = cos(pkin(5));
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t9 = -t25 * t23 + t26 * t24;
t11 = -t23 * qJDD(1) + t24 * qJDD(2);
t12 = t24 * qJDD(1) + t23 * qJDD(2);
t29 = t26 * t11 - t25 * t12;
t28 = t25 * t11 + t26 * t12;
t13 = -t23 * qJD(1) + t24 * qJD(2);
t14 = t24 * qJD(1) + t23 * qJD(2);
t4 = t26 * t13 - t25 * t14;
t5 = t25 * t13 + t26 * t14;
t10 = -t26 * t23 - t25 * t24;
t6 = t9 * qJD(4);
t27 = -t6 * qJD(4) + t10 * qJDD(4);
t22 = qJDD(1) - g(2);
t21 = qJDD(3) + g(3);
t20 = pkin(5) + qJ(4);
t19 = cos(t20);
t18 = sin(t20);
t7 = t10 * qJD(4);
t3 = t7 * qJD(4) + t9 * qJDD(4);
t2 = -t5 * qJD(4) + t29;
t1 = t4 * qJD(4) + t28;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t23 + t12 * t24 - g(2), 0, 0, 0, 0, 0, 0, t27, -t3, 0, t1 * t9 + t2 * t10 - t4 * t6 + t5 * t7 - g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t24 + t12 * t23 - g(1), 0, 0, 0, 0, 0, 0, t3, t27, 0, -t1 * t10 + t2 * t9 + t4 * t7 + t5 * t6 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -g(1) * t19 + g(2) * t18 + t29, g(1) * t18 + g(2) * t19 - t28, 0, 0;];
tau_reg  = t8;
