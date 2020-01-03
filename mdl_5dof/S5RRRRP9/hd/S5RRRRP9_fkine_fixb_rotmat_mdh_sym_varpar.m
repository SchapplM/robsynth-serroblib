% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRRP9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:39
% EndTime: 2019-12-31 22:03:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (116->47), mult. (132->51), div. (0->0), fcn. (187->8), ass. (0->31)
t21 = qJ(3) + qJ(4);
t15 = sin(t21);
t23 = sin(qJ(2));
t39 = t23 * t15;
t22 = sin(qJ(3));
t24 = sin(qJ(1));
t38 = t24 * t22;
t11 = t24 * t23;
t26 = cos(qJ(2));
t37 = t24 * t26;
t27 = cos(qJ(1));
t12 = t27 * t23;
t36 = t27 * t26;
t20 = pkin(5) + 0;
t35 = t24 * pkin(1) + 0;
t34 = t27 * pkin(1) + t24 * pkin(6) + 0;
t25 = cos(qJ(3));
t13 = t25 * pkin(3) + pkin(2);
t28 = -pkin(8) - pkin(7);
t33 = t23 * t13 + t26 * t28 + t20;
t32 = pkin(2) * t26 + pkin(7) * t23;
t31 = -t27 * pkin(6) + t35;
t30 = pkin(3) * t38 - t28 * t12 + t13 * t36 + t34;
t29 = -t28 * t11 + t13 * t37 + (-pkin(3) * t22 - pkin(6)) * t27 + t35;
t16 = cos(t21);
t8 = t23 * t16;
t4 = t24 * t15 + t16 * t36;
t3 = t15 * t36 - t24 * t16;
t2 = -t27 * t15 + t16 * t37;
t1 = t15 * t37 + t27 * t16;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t24, 0, 0; t24, t27, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t36, -t12, t24, t34; t37, -t11, -t27, t31; t23, t26, 0, t20; 0, 0, 0, 1; t25 * t36 + t38, -t22 * t36 + t24 * t25, t12, t32 * t27 + t34; -t27 * t22 + t25 * t37, -t22 * t37 - t27 * t25, t11, t32 * t24 + t31; t23 * t25, -t23 * t22, -t26, t23 * pkin(2) - t26 * pkin(7) + t20; 0, 0, 0, 1; t4, -t3, t12, t30; t2, -t1, t11, t29; t8, -t39, -t26, t33; 0, 0, 0, 1; t4, t12, t3, t4 * pkin(4) + t3 * qJ(5) + t30; t2, t11, t1, t2 * pkin(4) + t1 * qJ(5) + t29; t8, -t26, t39, (pkin(4) * t16 + qJ(5) * t15) * t23 + t33; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
