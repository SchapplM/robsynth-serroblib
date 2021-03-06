% Calculate minimal parameter regressor of joint inertia matrix for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% 
% Output:
% MM_reg [((3+1)*3/2)x7]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S3PRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_inertiaJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_inertiaJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:25:01
% EndTime: 2019-05-04 18:25:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (6->4), mult. (14->8), div. (0->0), fcn. (22->4), ass. (0->9)
t3 = sin(qJ(3));
t8 = t3 * pkin(2);
t5 = cos(qJ(3));
t7 = t5 * pkin(2);
t6 = cos(qJ(2));
t4 = sin(qJ(2));
t2 = -t3 * t6 - t5 * t4;
t1 = -t3 * t4 + t5 * t6;
t9 = [1, 0, 0, 0, 0, 0, 0; 0, 0, t6, -t4, 0, t1, t2; 0, 1, 0, 0, 1, 0.2e1 * t7, -0.2e1 * t8; 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 1, t7, -t8; 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
