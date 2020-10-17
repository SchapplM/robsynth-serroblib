% Calculate minimal parameter regressor of joint inertia matrix for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% 
% Output:
% MM_reg [((3+1)*3/2)x7]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S3PRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_inertiaJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_inertiaJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:23:49
% EndTime: 2019-05-04 18:23:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (6->6), mult. (8->8), div. (0->0), fcn. (9->2), ass. (0->3)
t2 = cos(qJ(2));
t1 = sin(qJ(2));
t3 = [1, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2; 0, 0, t2, -t1, t2, t1, t2 * pkin(2) + t1 * qJ(3); 0, 1, 0, 0, 0.2e1 * pkin(2), 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2; 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, -1, 0, -pkin(2); 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
