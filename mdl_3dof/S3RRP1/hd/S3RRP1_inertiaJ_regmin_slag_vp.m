% Calculate minimal parameter regressor of joint inertia matrix for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% MM_reg [((3+1)*3/2)x9]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S3RRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_inertiaJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_inertiaJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:31:25
% EndTime: 2019-05-04 18:31:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (14->10), mult. (29->14), div. (0->0), fcn. (13->2), ass. (0->7)
t9 = cos(qJ(2)) * pkin(1);
t8 = 2 * pkin(2);
t7 = 2 * qJ(3);
t4 = sin(qJ(2)) * pkin(1);
t2 = pkin(2) + t9;
t1 = t4 + qJ(3);
t3 = [1, 0, 0, 1, 0.2e1 * t9, -0.2e1 * t4, 0.2e1 * t2, 0.2e1 * t1, t1 ^ 2 + t2 ^ 2; 0, 0, 0, 1, t9, -t4, t8 + t9, t7 + t4, t2 * pkin(2) + t1 * qJ(3); 0, 0, 0, 1, 0, 0, t8, t7, pkin(2) ^ 2 + qJ(3) ^ 2; 0, 0, 0, 0, 0, 0, -1, 0, -t2; 0, 0, 0, 0, 0, 0, -1, 0, -pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
