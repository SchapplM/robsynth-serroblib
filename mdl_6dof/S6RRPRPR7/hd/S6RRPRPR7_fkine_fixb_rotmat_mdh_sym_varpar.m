% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:04:27
% EndTime: 2018-11-23 17:04:27
% DurationCPUTime: 0.15s
% Computational Cost: add. (169->60), mult. (217->65), div. (0->0), fcn. (293->10), ass. (0->42)
t34 = sin(qJ(1));
t37 = cos(qJ(2));
t19 = t34 * t37;
t33 = sin(qJ(2));
t52 = qJ(3) * t33;
t57 = pkin(2) * t19 + t34 * t52;
t32 = sin(qJ(4));
t56 = t33 * t32;
t55 = t34 * t33;
t28 = qJ(4) + pkin(10);
t22 = sin(t28);
t54 = t37 * t22;
t38 = cos(qJ(1));
t53 = t38 * t33;
t20 = t38 * t37;
t29 = pkin(6) + 0;
t51 = pkin(4) * t56;
t50 = t34 * pkin(1) + 0;
t49 = t33 * pkin(2) + t29;
t48 = t38 * pkin(1) + t34 * pkin(7) + 0;
t47 = t50 + t57;
t23 = cos(t28);
t5 = t33 * t22 + t37 * t23;
t36 = cos(qJ(4));
t46 = -t37 * t32 + t33 * t36;
t45 = t37 * t36 + t56;
t44 = -t38 * pkin(7) + t50;
t43 = pkin(2) * t20 + t38 * t52 + t48;
t42 = -t37 * qJ(3) + t49;
t21 = t36 * pkin(4) + pkin(3);
t30 = -qJ(5) - pkin(8);
t41 = t21 * t20 + t34 * t30 + t38 * t51 + t43;
t40 = t33 * t21 + (-pkin(4) * t32 - qJ(3)) * t37 + t49;
t39 = t34 * t51 + t21 * t19 + (-pkin(7) - t30) * t38 + t47;
t35 = cos(qJ(6));
t31 = sin(qJ(6));
t6 = t33 * t23 - t54;
t4 = t5 * t38;
t3 = t22 * t20 - t23 * t53;
t2 = t5 * t34;
t1 = -t23 * t55 + t34 * t54;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t38, -t34, 0, 0; t34, t38, 0, 0; 0, 0, 1, t29; 0, 0, 0, 1; t20, -t53, t34, t48; t19, -t55, -t38, t44; t33, t37, 0, t29; 0, 0, 0, 1; t20, t34, t53, t43; t19, -t38, t55, t44 + t57; t33, 0, -t37, t42; 0, 0, 0, 1; t45 * t38, t46 * t38, -t34, pkin(3) * t20 - t34 * pkin(8) + t43; t45 * t34, t46 * t34, t38, pkin(3) * t19 + (-pkin(7) + pkin(8)) * t38 + t47; t46, -t45, 0, t33 * pkin(3) + t42; 0, 0, 0, 1; t4, -t3, -t34, t41; t2, -t1, t38, t39; t6, -t5, 0, t40; 0, 0, 0, 1; -t34 * t31 + t4 * t35, -t4 * t31 - t34 * t35, t3, t4 * pkin(5) + t3 * pkin(9) + t41; t2 * t35 + t38 * t31, -t2 * t31 + t38 * t35, t1, t2 * pkin(5) + t1 * pkin(9) + t39; t6 * t35, -t6 * t31, t5, t6 * pkin(5) + t5 * pkin(9) + t40; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
