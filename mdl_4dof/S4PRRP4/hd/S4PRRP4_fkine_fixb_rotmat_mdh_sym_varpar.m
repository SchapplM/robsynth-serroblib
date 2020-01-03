% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4PRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:41
% EndTime: 2019-12-31 16:27:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (62->21), mult. (32->16), div. (0->0), fcn. (56->6), ass. (0->19)
t14 = sin(qJ(3));
t11 = pkin(6) + qJ(2);
t6 = sin(t11);
t23 = t6 * t14;
t7 = cos(t11);
t22 = t7 * t14;
t12 = sin(pkin(6));
t21 = t12 * pkin(1) + 0;
t13 = cos(pkin(6));
t20 = t13 * pkin(1) + 0;
t19 = qJ(1) + 0;
t8 = pkin(4) + t19;
t18 = t7 * pkin(2) + t6 * pkin(5) + t20;
t15 = cos(qJ(3));
t17 = pkin(3) * t15 + qJ(4) * t14;
t16 = t6 * pkin(2) - t7 * pkin(5) + t21;
t2 = t7 * t15;
t1 = t6 * t15;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t13, -t12, 0, 0; t12, t13, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; t7, -t6, 0, t20; t6, t7, 0, t21; 0, 0, 1, t8; 0, 0, 0, 1; t2, -t22, t6, t18; t1, -t23, -t7, t16; t14, t15, 0, t8; 0, 0, 0, 1; t2, t6, t22, t17 * t7 + t18; t1, -t7, t23, t17 * t6 + t16; t14, 0, -t15, t14 * pkin(3) - t15 * qJ(4) + t8; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
