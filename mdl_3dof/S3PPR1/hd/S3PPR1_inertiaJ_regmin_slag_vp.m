% Calculate minimal parameter regressor of joint inertia matrix for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% 
% Output:
% MM_reg [((3+1)*3/2)x5]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S3PPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_inertiaJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_inertiaJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:21:18
% EndTime: 2019-05-04 18:21:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
t2 = cos(qJ(3));
t1 = sin(qJ(3));
t3 = [1, 1, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 0, 0, -t1, -t2; 0, 0, 0, t2, -t1; 0, 0, 1, 0, 0;];
MM_reg  = t3;
