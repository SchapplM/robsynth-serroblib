% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t9 = sin(pkin(6));
t17 = pkin(1) * t9;
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t10 = cos(pkin(6));
t8 = t10 * pkin(1) + pkin(2);
t4 = -t11 * t8 - t12 * t17;
t16 = t11 * t17 - t12 * t8;
t14 = 2 * pkin(3);
t13 = 2 * qJ(4);
t2 = -pkin(3) + t16;
t1 = qJ(4) - t4;
t3 = [1, 0, 0 (t10 ^ 2 + t9 ^ 2) * pkin(1) ^ 2, 1, -0.2e1 * t16, 0.2e1 * t4, -0.2e1 * t2, 0.2e1 * t1, t1 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 1, -t16, t4, t14 - t16, t13 - t4, -t2 * pkin(3) + t1 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, t14, t13, pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, -1, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -1, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
