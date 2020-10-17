% Calculate minimal parameter regressor of joint inertia matrix for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:27
% DurationCPUTime: 0.07s
% Computational Cost: add. (7->7), mult. (13->12), div. (0->0), fcn. (21->4), ass. (0->6)
t3 = cos(qJ(4));
t5 = 0.2e1 * t3;
t4 = cos(qJ(3));
t2 = sin(qJ(3));
t1 = sin(qJ(4));
t6 = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t4, -t2, 0, 0, 0, 0, 0, t4 * t3, -t4 * t1; 0, 0, 1, 0, 0, t1 ^ 2, t1 * t5, 0, 0, 0, pkin(3) * t5, -0.2e1 * pkin(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t2, -t3 * t2; 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -t1 * pkin(5), -t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
