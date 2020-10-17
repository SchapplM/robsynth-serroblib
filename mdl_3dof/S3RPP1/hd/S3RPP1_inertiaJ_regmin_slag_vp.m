% Calculate minimal parameter regressor of joint inertia matrix for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% MM_reg [((3+1)*3/2)x9]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S3RPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_inertiaJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_inertiaJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:27:37
% EndTime: 2019-05-04 18:27:37
% DurationCPUTime: 0.05s
% Computational Cost: add. (6->5), mult. (9->6), div. (0->0), fcn. (0->0), ass. (0->4)
t4 = (qJ(2) ^ 2);
t3 = 2 * qJ(2);
t1 = (pkin(1) + qJ(3));
t2 = [1, 0, 0, -2 * pkin(1), t3, pkin(1) ^ 2 + t4, t3, 2 * t1, t1 ^ 2 + t4; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t1; 0, 0, 0, 0, 0, 1, 0, 0, 1; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2); 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;
