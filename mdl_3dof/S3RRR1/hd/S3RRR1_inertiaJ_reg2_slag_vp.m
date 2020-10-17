% Calculate inertial parameters regressor of joint inertia matrix for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% MM_reg [((3+1)*3/2)x(3*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S3RRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaJ_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:32:47
% EndTime: 2019-05-04 18:32:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (28->13), mult. (76->28), div. (0->0), fcn. (56->4), ass. (0->13)
t7 = sin(qJ(3));
t15 = t7 * pkin(2);
t8 = sin(qJ(2));
t14 = t8 * pkin(1);
t9 = cos(qJ(3));
t13 = t9 * t14;
t10 = cos(qJ(2));
t6 = t10 * pkin(1);
t4 = t6 + pkin(2);
t1 = -t7 * t14 + t9 * t4;
t5 = t9 * pkin(2);
t2 = t7 * t4 + t13;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t6, -0.2e1 * t14, 0 (t10 ^ 2 + t8 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t1, -0.2e1 * t2, 0, t1 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t6, -t14, 0, 0, 0, 0, 0, 0, 0, 1, t1 + t5, -t13 + (-pkin(2) - t4) * t7, 0 (t1 * t9 + t2 * t7) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t5, -0.2e1 * t15, 0 (t7 ^ 2 + t9 ^ 2) * pkin(2) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t5, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t3;
