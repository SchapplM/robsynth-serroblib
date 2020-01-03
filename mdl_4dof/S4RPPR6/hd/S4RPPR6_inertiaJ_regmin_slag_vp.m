% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t13 = sin(pkin(6));
t14 = cos(pkin(6));
t22 = t13 ^ 2 + t14 ^ 2;
t18 = t13 * qJ(3) + pkin(1);
t21 = 0.2e1 * (pkin(2) + pkin(3)) * t14 + 0.2e1 * t18;
t20 = -0.2e1 * t13;
t19 = t22 * qJ(2) ^ 2;
t16 = cos(qJ(4));
t15 = sin(qJ(4));
t9 = t13 * qJ(2);
t7 = (-pkin(5) + qJ(2)) * t14;
t6 = -t13 * pkin(5) + t9;
t5 = -t14 * pkin(2) - t18;
t4 = 0.2e1 * t22 * qJ(2);
t3 = t13 * t16 - t14 * t15;
t2 = t13 * t15 + t14 * t16;
t1 = [1, 0, 0, 0.2e1 * pkin(1) * t14, pkin(1) * t20, t4, pkin(1) ^ 2 + t19, -0.2e1 * t5 * t14, t4, t5 * t20, t5 ^ 2 + t19, t3 ^ 2, -0.2e1 * t3 * t2, 0, 0, 0, t2 * t21, t3 * t21; 0, 0, 0, -t14, t13, 0, -pkin(1), -t14, 0, -t13, t5, 0, 0, 0, 0, 0, -t2, -t3; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, 0, -t15 * t7 + t16 * t6, -t15 * t6 - t16 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
