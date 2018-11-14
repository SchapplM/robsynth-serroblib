% Calculate inertial parameters regressor of coriolis joint torque vector for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% 
% Output:
% taug_reg [3x(3*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S3PPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:10:05
% EndTime: 2018-11-14 10:10:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (28->11), mult. (68->19), div. (0->0), fcn. (44->2), ass. (0->10)
t5 = sin(qJ(3));
t7 = qJD(3) ^ 2;
t9 = t5 * t7;
t6 = cos(qJ(3));
t8 = t6 * t7;
t4 = t6 * qJD(1) + t5 * qJD(2);
t3 = -t5 * qJD(1) + t6 * qJD(2);
t2 = t4 * qJD(3);
t1 = t3 * qJD(3);
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t9, 0, t1 * t6 + t2 * t5 + (-t3 * t6 - t4 * t5) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t8, 0, t1 * t5 - t2 * t6 + (-t3 * t5 + t4 * t6) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t10;
