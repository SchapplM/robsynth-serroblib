% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2018-11-23 15:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRRPP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:20:45
% EndTime: 2018-11-23 15:20:45
% DurationCPUTime: 0.16s
% Computational Cost: add. (690->70), mult. (797->70), div. (0->0), fcn. (903->14), ass. (0->57)
t66 = pkin(6) + qJ(2);
t77 = sin(t66) / 0.2e1;
t67 = pkin(6) - qJ(2);
t60 = sin(t67);
t31 = t77 - t60 / 0.2e1;
t41 = sin(pkin(10));
t43 = cos(pkin(10));
t47 = cos(qJ(2));
t22 = t43 * t31 + t41 * t47;
t45 = sin(qJ(3));
t42 = sin(pkin(6));
t73 = cos(qJ(3));
t64 = t42 * t73;
t12 = t22 * t45 + t43 * t64;
t76 = t12 * pkin(9);
t24 = -t41 * t31 + t43 * t47;
t14 = t24 * t45 - t41 * t64;
t75 = t14 * pkin(9);
t59 = cos(t67) / 0.2e1;
t61 = cos(t66);
t32 = t59 - t61 / 0.2e1;
t68 = cos(pkin(6));
t25 = t32 * t45 - t68 * t73;
t74 = t25 * pkin(9);
t72 = cos(qJ(4));
t71 = t41 * t42;
t70 = t43 * t42;
t69 = pkin(9) - qJ(6);
t65 = qJ(1) + 0;
t63 = t43 * pkin(1) + pkin(7) * t71 + 0;
t62 = t68 * pkin(7) + t65;
t58 = t41 * pkin(1) - pkin(7) * t70 + 0;
t46 = sin(qJ(2));
t52 = t59 + t61 / 0.2e1;
t23 = t41 * t52 + t43 * t46;
t57 = t24 * pkin(2) + t23 * pkin(8) + t63;
t30 = t77 + t60 / 0.2e1;
t56 = t32 * pkin(2) - t30 * pkin(8) + t62;
t15 = t24 * t73 + t45 * t71;
t55 = t15 * pkin(3) + t57;
t26 = t32 * t73 + t68 * t45;
t54 = t26 * pkin(3) + t56;
t21 = t41 * t46 - t43 * t52;
t53 = t22 * pkin(2) + t21 * pkin(8) + t58;
t13 = t22 * t73 - t45 * t70;
t51 = t13 * pkin(3) + t53;
t44 = sin(qJ(4));
t5 = t15 * t44 - t23 * t72;
t6 = t15 * t72 + t23 * t44;
t50 = t6 * pkin(4) + t5 * qJ(5) + t55;
t8 = t26 * t44 + t30 * t72;
t9 = t26 * t72 - t30 * t44;
t49 = t9 * pkin(4) + t8 * qJ(5) + t54;
t3 = t13 * t44 - t21 * t72;
t4 = t13 * t72 + t21 * t44;
t48 = t4 * pkin(4) + t3 * qJ(5) + t51;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t43, -t41, 0, 0; t41, t43, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t24, -t23, t71, t63; t22, -t21, -t70, t58; t32, t30, t68, t62; 0, 0, 0, 1; t15, -t14, t23, t57; t13, -t12, t21, t53; t26, -t25, -t30, t56; 0, 0, 0, 1; t6, -t5, t14, t55 + t75; t4, -t3, t12, t51 + t76; t9, -t8, t25, t54 + t74; 0, 0, 0, 1; t6, t14, t5, t50 + t75; t4, t12, t3, t48 + t76; t9, t25, t8, t49 + t74; 0, 0, 0, 1; t6, t5, -t14, t6 * pkin(5) + t69 * t14 + t50; t4, t3, -t12, t4 * pkin(5) + t69 * t12 + t48; t9, t8, -t25, t9 * pkin(5) + t69 * t25 + t49; 0, 0, 0, 1;];
T_ges = t1;
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
