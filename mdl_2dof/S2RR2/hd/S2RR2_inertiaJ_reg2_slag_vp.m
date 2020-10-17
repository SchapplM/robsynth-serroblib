% Calculate inertial parameters regressor of joint inertia matrix for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% MM_reg [((2+1)*2/2)x(2*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S2RR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_inertiaJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_inertiaJ_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:19:25
% EndTime: 2019-05-04 18:19:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (4->3), mult. (18->10), div. (0->0), fcn. (12->2), ass. (0->6)
t3 = sin(qJ(2));
t1 = t3 ^ 2;
t4 = cos(qJ(2));
t2 = t4 ^ 2;
t6 = t1 + t2;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t1, 0.2e1 * t3 * t4, 0, t2, 0, 0, 0, 0, 0.2e1 * t6 * pkin(1), t6 * pkin(1) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, t4, 0, -t3 * pkin(1), -t4 * pkin(1), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t5;
