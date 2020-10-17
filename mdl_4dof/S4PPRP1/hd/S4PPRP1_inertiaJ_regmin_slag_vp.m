% Calculate minimal parameter regressor of joint inertia matrix for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:43:14
% EndTime: 2019-05-04 18:43:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (6->6), mult. (8->8), div. (0->0), fcn. (9->2), ass. (0->3)
t2 = cos(qJ(3));
t1 = sin(qJ(3));
t3 = [1, 1, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t2, -t1, t2, t1, t2 * pkin(3) + t1 * qJ(4); 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, -1, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
