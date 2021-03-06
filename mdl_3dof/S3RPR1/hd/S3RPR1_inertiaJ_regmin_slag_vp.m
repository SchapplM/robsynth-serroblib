% Calculate minimal parameter regressor of joint inertia matrix for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% MM_reg [((3+1)*3/2)x9]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S3RPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_inertiaJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_inertiaJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:28:55
% EndTime: 2019-05-04 18:28:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (12->9), mult. (16->10), div. (0->0), fcn. (12->2), ass. (0->6)
t5 = -pkin(1) - pkin(2);
t4 = cos(qJ(3));
t3 = sin(qJ(3));
t2 = qJ(2) * t4 + t3 * t5;
t1 = qJ(2) * t3 - t4 * t5;
t6 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, -1, 0, -pkin(1), 0, -t4, t3; 0, 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, -1, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
