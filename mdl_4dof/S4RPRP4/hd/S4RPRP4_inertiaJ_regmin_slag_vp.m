% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRP4
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
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t9 = sin(qJ(3));
t18 = 0.2e1 * t9;
t10 = cos(qJ(3));
t17 = -0.2e1 * t10;
t7 = sin(pkin(6));
t3 = t7 * pkin(1) + pkin(5);
t16 = t9 * t3;
t5 = t9 ^ 2;
t15 = t10 ^ 2 + t5;
t14 = t10 * t3;
t8 = cos(pkin(6));
t4 = -t8 * pkin(1) - pkin(2);
t13 = -t9 * pkin(3) + t10 * qJ(4);
t12 = t10 * pkin(3) + t9 * qJ(4);
t1 = -t12 + t4;
t2 = [1, 0, 0, (t7 ^ 2 + t8 ^ 2) * pkin(1) ^ 2, t5, t10 * t18, 0, 0, 0, t4 * t17, t4 * t18, t1 * t17, 0.2e1 * t15 * t3, -0.2e1 * t1 * t9, t15 * t3 ^ 2 + t1 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, t9, t10, 0, -t16, -t14, -t16, t13, t14, t13 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, t10, 0, t9, t12; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
