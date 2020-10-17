% Calculate minimal parameter regressor of joint inertia matrix for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x6]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:38:12
% EndTime: 2019-05-04 18:38:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (6->5), mult. (8->6), div. (0->0), fcn. (15->4), ass. (0->6)
t5 = cos(qJ(4));
t4 = sin(qJ(4));
t3 = cos(pkin(5));
t2 = sin(pkin(5));
t1 = t2 ^ 2 + t3 ^ 2;
t6 = [1, t1, t1, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 1, 0, 0, 0; 0, 0, -t3, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, -t2 * t4 - t3 * t5, -t2 * t5 + t3 * t4; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t5, -t4; 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
