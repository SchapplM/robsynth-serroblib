% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPRP8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:46:19
% EndTime: 2018-11-23 17:46:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (154->61), mult. (286->62), div. (0->0), fcn. (388->8), ass. (0->42)
t32 = sin(qJ(2));
t35 = cos(qJ(3));
t19 = t32 * t35;
t31 = sin(qJ(3));
t52 = t32 * t31;
t54 = pkin(3) * t19 + qJ(4) * t52;
t30 = sin(qJ(5));
t53 = t30 * t31;
t33 = sin(qJ(1));
t20 = t33 * t32;
t36 = cos(qJ(2));
t51 = t33 * t36;
t37 = cos(qJ(1));
t22 = t37 * t32;
t50 = t37 * t36;
t28 = pkin(6) + 0;
t49 = t32 * pkin(2) + t28;
t48 = pkin(5) * t30 + qJ(4);
t47 = t37 * pkin(1) + t33 * pkin(7) + 0;
t46 = t33 * pkin(1) - t37 * pkin(7) + 0;
t45 = t49 + t54;
t44 = pkin(2) * t50 + pkin(8) * t22 + t47;
t43 = -t36 * pkin(8) + t49;
t12 = t33 * t31 + t35 * t50;
t42 = t12 * pkin(3) + t44;
t41 = pkin(2) * t51 + pkin(8) * t20 + t46;
t10 = -t37 * t31 + t35 * t51;
t40 = t10 * pkin(3) + t41;
t11 = t31 * t50 - t33 * t35;
t39 = t11 * qJ(4) + t42;
t9 = t31 * t51 + t37 * t35;
t38 = t9 * qJ(4) + t40;
t34 = cos(qJ(5));
t29 = -qJ(6) - pkin(9);
t23 = t34 * pkin(5) + pkin(4);
t6 = (t34 * t35 + t53) * t32;
t5 = (-t30 * t35 + t31 * t34) * t32;
t4 = t11 * t30 + t12 * t34;
t3 = t11 * t34 - t12 * t30;
t2 = t10 * t34 + t9 * t30;
t1 = -t10 * t30 + t9 * t34;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t37, -t33, 0, 0; t33, t37, 0, 0; 0, 0, 1, t28; 0, 0, 0, 1; t50, -t22, t33, t47; t51, -t20, -t37, t46; t32, t36, 0, t28; 0, 0, 0, 1; t12, -t11, t22, t44; t10, -t9, t20, t41; t19, -t52, -t36, t43; 0, 0, 0, 1; t12, t22, t11, t39; t10, t20, t9, t38; t19, -t36, t52, t43 + t54; 0, 0, 0, 1; t4, t3, -t22, t12 * pkin(4) - pkin(9) * t22 + t39; t2, t1, -t20, t10 * pkin(4) - pkin(9) * t20 + t38; t6, t5, t36, pkin(4) * t19 + (-pkin(8) + pkin(9)) * t36 + t45; 0, 0, 0, 1; t4, t3, -t22, t48 * t11 + t12 * t23 + t29 * t22 + t42; t2, t1, -t20, t10 * t23 + t29 * t20 + t48 * t9 + t40; t6, t5, t36 (-pkin(8) - t29) * t36 + (pkin(5) * t53 + t23 * t35) * t32 + t45; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
