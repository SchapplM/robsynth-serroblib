% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t11 = cos(pkin(8));
t20 = -0.2e1 * t11 * pkin(4) - (2 * pkin(3));
t10 = sin(pkin(8));
t19 = t10 ^ 2 + t11 ^ 2;
t18 = pkin(6) + qJ(4);
t17 = t19 * qJ(4);
t12 = sin(qJ(5));
t14 = cos(qJ(5));
t2 = t14 * t10 + t12 * t11;
t1 = t12 * t10 - t14 * t11;
t15 = cos(qJ(3));
t13 = sin(qJ(3));
t4 = t18 * t11;
t3 = t18 * t10;
t5 = [1, 1, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, t19 * t13 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t15, -t13, t15 * t11, -t15 * t10, t19 * t13, t15 * pkin(3) + t13 * t17, 0, 0, 0, 0, 0, -t15 * t1, -t15 * t2; 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t11, -0.2e1 * pkin(3) * t10, 0.2e1 * t17, t19 * qJ(4) ^ 2 + (pkin(3) ^ 2), t2 ^ 2, -0.2e1 * t2 * t1, 0, 0, 0, t1 * t20, t2 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t11, t10, 0, -pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t13, t1 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, 0, -t12 * t4 - t14 * t3, t12 * t3 - t14 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
