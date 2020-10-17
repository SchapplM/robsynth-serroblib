% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:04:38
% EndTime: 2019-05-04 19:04:39
% DurationCPUTime: 0.08s
% Computational Cost: add. (16->8), mult. (42->15), div. (0->0), fcn. (34->4), ass. (0->11)
t7 = sin(qJ(4));
t13 = t7 * pkin(3);
t12 = sin(qJ(3)) * pkin(2);
t9 = cos(qJ(4));
t11 = t9 * t12;
t6 = cos(qJ(3)) * pkin(2);
t4 = t6 + pkin(3);
t1 = -t7 * t12 + t9 * t4;
t5 = t9 * pkin(3);
t2 = -t7 * t4 - t11;
t3 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t6, -0.2e1 * t12, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t6, -t12, 1, t1 + t5, -t11 + (-pkin(3) - t4) * t7; 0, 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t5, -0.2e1 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, t5, -t13; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
