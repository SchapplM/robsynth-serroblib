% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:24:20
% EndTime: 2018-11-23 17:24:20
% DurationCPUTime: 0.16s
% Computational Cost: add. (169->60), mult. (217->65), div. (0->0), fcn. (293->10), ass. (0->42)
t33 = sin(qJ(1));
t36 = cos(qJ(2));
t18 = t33 * t36;
t32 = sin(qJ(2));
t52 = qJ(3) * t32;
t57 = pkin(2) * t18 + t33 * t52;
t31 = sin(qJ(4));
t56 = t32 * t31;
t55 = t33 * t32;
t29 = qJ(4) + qJ(5);
t22 = sin(t29);
t54 = t36 * t22;
t37 = cos(qJ(1));
t53 = t37 * t32;
t19 = t37 * t36;
t28 = pkin(6) + 0;
t51 = pkin(4) * t56;
t50 = t33 * pkin(1) + 0;
t49 = t32 * pkin(2) + t28;
t48 = t37 * pkin(1) + t33 * pkin(7) + 0;
t47 = t50 + t57;
t23 = cos(t29);
t5 = t32 * t22 + t36 * t23;
t35 = cos(qJ(4));
t46 = -t36 * t31 + t32 * t35;
t45 = t36 * t35 + t56;
t44 = -t37 * pkin(7) + t50;
t43 = pkin(2) * t19 + t37 * t52 + t48;
t42 = -t36 * qJ(3) + t49;
t20 = t35 * pkin(4) + pkin(3);
t38 = -pkin(9) - pkin(8);
t41 = t20 * t19 + t33 * t38 + t37 * t51 + t43;
t40 = t32 * t20 + (-pkin(4) * t31 - qJ(3)) * t36 + t49;
t39 = t33 * t51 + t20 * t18 + (-pkin(7) - t38) * t37 + t47;
t34 = cos(qJ(6));
t30 = sin(qJ(6));
t6 = t32 * t23 - t54;
t4 = t5 * t37;
t3 = t22 * t19 - t23 * t53;
t2 = t5 * t33;
t1 = -t23 * t55 + t33 * t54;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t37, -t33, 0, 0; t33, t37, 0, 0; 0, 0, 1, t28; 0, 0, 0, 1; t19, -t53, t33, t48; t18, -t55, -t37, t44; t32, t36, 0, t28; 0, 0, 0, 1; t19, t33, t53, t43; t18, -t37, t55, t44 + t57; t32, 0, -t36, t42; 0, 0, 0, 1; t45 * t37, t46 * t37, -t33, pkin(3) * t19 - t33 * pkin(8) + t43; t45 * t33, t46 * t33, t37, pkin(3) * t18 + (-pkin(7) + pkin(8)) * t37 + t47; t46, -t45, 0, t32 * pkin(3) + t42; 0, 0, 0, 1; t4, -t3, -t33, t41; t2, -t1, t37, t39; t6, -t5, 0, t40; 0, 0, 0, 1; -t33 * t30 + t4 * t34, -t4 * t30 - t33 * t34, t3, t4 * pkin(5) + t3 * pkin(10) + t41; t2 * t34 + t37 * t30, -t2 * t30 + t37 * t34, t1, t2 * pkin(5) + t1 * pkin(10) + t39; t6 * t34, -t6 * t30, t5, t6 * pkin(5) + t5 * pkin(10) + t40; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
