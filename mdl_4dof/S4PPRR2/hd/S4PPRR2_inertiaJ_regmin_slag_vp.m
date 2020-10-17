% Calculate minimal parameter regressor of joint inertia matrix for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:50:52
% EndTime: 2019-05-04 18:50:52
% DurationCPUTime: 0.07s
% Computational Cost: add. (16->7), mult. (36->14), div. (0->0), fcn. (54->6), ass. (0->13)
t7 = sin(qJ(4));
t12 = t7 * pkin(3);
t9 = cos(qJ(4));
t11 = t9 * pkin(3);
t10 = cos(qJ(3));
t8 = sin(qJ(3));
t6 = cos(pkin(6));
t5 = sin(pkin(6));
t4 = t10 * t5 + t8 * t6;
t3 = t10 * t6 - t8 * t5;
t2 = -t7 * t3 - t9 * t4;
t1 = t9 * t3 - t7 * t4;
t13 = [1, t5 ^ 2 + t6 ^ 2, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, t3, -t4, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 1, 0.2e1 * t11, -0.2e1 * t12; 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, t11, -t12; 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
