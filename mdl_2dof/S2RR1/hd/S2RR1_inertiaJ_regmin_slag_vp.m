% Calculate minimal parameter regressor of joint inertia matrix for
% S2RR1
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
% MM_reg [((2+1)*2/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S2RR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_inertiaJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_inertiaJ_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:10
% EndTime: 2020-01-03 11:19:10
% DurationCPUTime: 0.06s
% Computational Cost: add. (2->2), mult. (5->5), div. (0->0), fcn. (7->2), ass. (0->3)
t2 = cos(qJ(2));
t1 = sin(qJ(2));
t3 = [1, 0, 0, t1 ^ 2, 0.2e1 * t1 * t2, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t1, -t2, 0, t1 * pkin(1), t2 * pkin(1); 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
