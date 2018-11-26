% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPRP9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:01:46
% EndTime: 2018-11-23 16:01:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (136->55), mult. (144->53), div. (0->0), fcn. (203->8), ass. (0->36)
t23 = sin(pkin(9));
t27 = sin(qJ(1));
t45 = t27 * t23;
t26 = sin(qJ(3));
t44 = t27 * t26;
t28 = cos(qJ(3));
t43 = t27 * t28;
t21 = pkin(9) + qJ(5);
t14 = sin(t21);
t42 = t28 * t14;
t29 = cos(qJ(1));
t41 = t29 * t23;
t40 = t29 * t26;
t22 = pkin(6) + 0;
t39 = t27 * pkin(1) + 0;
t38 = pkin(2) + t22;
t17 = t27 * pkin(7);
t37 = t17 + t39;
t36 = t29 * pkin(1) + t27 * qJ(2) + 0;
t35 = t29 * pkin(7) + t36;
t34 = pkin(3) * t26 - qJ(4) * t28;
t33 = -t29 * qJ(2) + t39;
t24 = cos(pkin(9));
t12 = t24 * pkin(4) + pkin(3);
t25 = -pkin(8) - qJ(4);
t32 = t28 * t12 - t26 * t25 + t38;
t31 = pkin(4) * t41 + t12 * t44 + t25 * t43 + t35;
t30 = pkin(4) * t45 + (-t12 * t26 - t25 * t28 - qJ(2)) * t29 + t37;
t15 = cos(t21);
t13 = t29 * t28;
t8 = t28 * t15;
t4 = t27 * t14 - t15 * t40;
t3 = t14 * t40 + t27 * t15;
t2 = t29 * t14 + t15 * t44;
t1 = t14 * t44 - t29 * t15;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; 0, -t29, t27, t36; 0, -t27, -t29, t33; 1, 0, 0, t22; 0, 0, 0, 1; t44, t43, t29, t35; -t40, -t13, t27, t17 + t33; t28, -t26, 0, t38; 0, 0, 0, 1; t24 * t44 + t41, -t23 * t44 + t29 * t24, -t43, t34 * t27 + t35; -t24 * t40 + t45, t23 * t40 + t27 * t24, t13 (-qJ(2) - t34) * t29 + t37; t28 * t24, -t28 * t23, t26, t28 * pkin(3) + t26 * qJ(4) + t38; 0, 0, 0, 1; t2, -t1, -t43, t31; t4, t3, t13, t30; t8, -t42, t26, t32; 0, 0, 0, 1; t2, -t43, t1, t2 * pkin(5) + t1 * qJ(6) + t31; t4, t13, -t3, t4 * pkin(5) - t3 * qJ(6) + t30; t8, t26, t42 (pkin(5) * t15 + qJ(6) * t14) * t28 + t32; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
