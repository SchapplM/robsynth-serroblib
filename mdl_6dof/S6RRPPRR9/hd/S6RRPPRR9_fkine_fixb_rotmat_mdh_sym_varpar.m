% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2018-11-23 16:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPPRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:54:32
% EndTime: 2018-11-23 16:54:32
% DurationCPUTime: 0.18s
% Computational Cost: add. (441->75), mult. (469->69), div. (0->0), fcn. (516->14), ass. (0->52)
t68 = -pkin(3) - pkin(8);
t32 = sin(pkin(6));
t37 = sin(qJ(1));
t26 = t37 * t32;
t41 = cos(qJ(1));
t67 = t41 * t32;
t66 = pkin(9) - qJ(3);
t36 = sin(qJ(2));
t62 = pkin(6) - qJ(2);
t52 = cos(t62) / 0.2e1;
t61 = pkin(6) + qJ(2);
t56 = cos(t61);
t45 = t52 + t56 / 0.2e1;
t10 = t37 * t36 - t41 * t45;
t65 = t10 * qJ(3);
t12 = t41 * t36 + t37 * t45;
t64 = t12 * qJ(3);
t51 = sin(t61) / 0.2e1;
t55 = sin(t62);
t18 = t51 + t55 / 0.2e1;
t63 = t18 * qJ(3);
t60 = pkin(7) + 0;
t59 = t37 * pkin(1) + 0;
t33 = cos(pkin(6));
t58 = t33 * pkin(8) + t60;
t57 = t41 * pkin(1) + pkin(8) * t26 + 0;
t19 = t52 - t56 / 0.2e1;
t54 = t19 * pkin(2) + t58;
t40 = cos(qJ(2));
t50 = t51 - t55 / 0.2e1;
t13 = -t37 * t50 + t41 * t40;
t53 = t13 * pkin(2) + t57;
t49 = -pkin(8) * t67 + t59;
t11 = t37 * t40 + t41 * t50;
t6 = t11 * pkin(2);
t48 = t11 * qJ(4) + t59 + t6;
t47 = t33 * pkin(3) + t19 * qJ(4) + t54;
t46 = pkin(3) * t26 + t13 * qJ(4) + t53;
t44 = pkin(4) * t26 - t66 * t12 + t46;
t43 = t33 * pkin(4) + t66 * t18 + t47;
t42 = -t66 * t10 + (-pkin(4) + t68) * t67 + t48;
t39 = cos(qJ(5));
t38 = cos(qJ(6));
t35 = sin(qJ(5));
t34 = sin(qJ(6));
t9 = t19 * t35 + t33 * t39;
t8 = -t19 * t39 + t33 * t35;
t4 = t11 * t35 - t39 * t67;
t3 = t11 * t39 + t35 * t67;
t2 = t13 * t35 + t39 * t26;
t1 = -t13 * t39 + t35 * t26;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t41, -t37, 0, 0; t37, t41, 0, 0; 0, 0, 1, t60; 0, 0, 0, 1; t13, -t12, t26, t57; t11, -t10, -t67, t49; t19, t18, t33, t58; 0, 0, 0, 1; t26, -t13, t12, t53 + t64; -t67, -t11, t10, t49 + t6 + t65; t33, -t19, -t18, t54 - t63; 0, 0, 0, 1; t26, t12, t13, t46 + t64; -t67, t10, t11, t68 * t67 + t48 + t65; t33, -t18, t19, t47 - t63; 0, 0, 0, 1; t2, -t1, -t12, t44; t4, t3, -t10, t42; t9, -t8, t18, t43; 0, 0, 0, 1; -t12 * t34 + t2 * t38, -t12 * t38 - t2 * t34, t1, t2 * pkin(5) + t1 * pkin(10) + t44; -t10 * t34 + t4 * t38, -t10 * t38 - t4 * t34, -t3, t4 * pkin(5) - t3 * pkin(10) + t42; t18 * t34 + t9 * t38, t18 * t38 - t9 * t34, t8, t9 * pkin(5) + t8 * pkin(10) + t43; 0, 0, 0, 1;];
T_ges = t5;
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
