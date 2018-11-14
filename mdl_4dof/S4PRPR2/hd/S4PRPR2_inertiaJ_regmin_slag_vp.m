% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x8]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MM_reg = S4PRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t8 = sin(pkin(6));
t15 = pkin(2) * t8;
t13 = cos(qJ(2));
t12 = cos(qJ(4));
t11 = sin(qJ(2));
t10 = sin(qJ(4));
t9 = cos(pkin(6));
t7 = t9 * pkin(2) + pkin(3);
t6 = t9 * t11 + t8 * t13;
t5 = -t8 * t11 + t9 * t13;
t4 = -t10 * t7 - t12 * t15;
t3 = -t10 * t15 + t12 * t7;
t2 = -t10 * t5 - t12 * t6;
t1 = -t10 * t6 + t12 * t5;
t14 = [1, 0, 0, 0, t5 ^ 2 + t6 ^ 2, 0, 0, 0; 0, 0, t13, -t11 (t5 * t9 + t6 * t8) * pkin(2), 0, t1, t2; 0, 1, 0, 0 (t8 ^ 2 + t9 ^ 2) * pkin(2) ^ 2, 1, 0.2e1 * t3, 0.2e1 * t4; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 1, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t14;
