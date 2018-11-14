% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:58:25
% EndTime: 2018-11-14 13:58:25
% DurationCPUTime: 0.09s
% Computational Cost: add. (61->18), mult. (144->23), div. (0->0), fcn. (90->2), ass. (0->15)
t14 = sin(qJ(3));
t15 = cos(qJ(3));
t24 = -t14 * qJD(1) + t15 * qJD(2);
t11 = t15 * qJD(1) + t14 * qJD(2);
t9 = t11 * qJD(3);
t23 = t11 * t14;
t16 = qJD(3) ^ 2;
t22 = t15 * t16;
t21 = qJD(3) * t14;
t8 = t24 * qJD(3);
t19 = t8 * t14;
t18 = t9 * t14 + t8 * t15;
t12 = t14 * t16;
t7 = qJD(3) * pkin(3) + t24;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t12, 0 (-t15 * t24 - t23) * qJD(3) + t18, 0, 0, 0, 0, 0, 0, -t22, t12, 0 (-t15 * t7 - t23) * qJD(3) + t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t22, 0, -t21 * t24 + t19, 0, 0, 0, 0, 0, 0, -t12, -t22, 0, -t7 * t21 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * pkin(3) + (-t24 + t7) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t1;
