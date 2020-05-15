% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRP11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:48:47
% EndTime: 2018-11-23 17:48:47
% DurationCPUTime: 0.19s
% Computational Cost: add. (572->74), mult. (646->75), div. (0->0), fcn. (727->14), ass. (0->60)
t76 = pkin(4) + pkin(9);
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t67 = pkin(6) - qJ(2);
t56 = cos(t67) / 0.2e1;
t66 = pkin(6) + qJ(2);
t60 = cos(t66);
t47 = t56 + t60 / 0.2e1;
t18 = t42 * t41 - t45 * t47;
t75 = t18 * pkin(9);
t20 = t45 * t41 + t42 * t47;
t74 = t20 * pkin(9);
t55 = sin(t66) / 0.2e1;
t59 = sin(t67);
t24 = t55 + t59 / 0.2e1;
t73 = t24 * pkin(9);
t43 = cos(qJ(5));
t72 = t43 * pkin(5) + t76;
t71 = cos(qJ(3));
t37 = sin(pkin(6));
t70 = t42 * t37;
t69 = t45 * t37;
t68 = cos(pkin(6));
t65 = pkin(7) + 0;
t64 = t37 * t71;
t63 = t68 * pkin(8) + t65;
t39 = sin(qJ(5));
t62 = pkin(5) * t39 + qJ(4);
t61 = t45 * pkin(1) + pkin(8) * t70 + 0;
t26 = t56 - t60 / 0.2e1;
t58 = t26 * pkin(2) + t63;
t25 = t55 - t59 / 0.2e1;
t44 = cos(qJ(2));
t21 = -t42 * t25 + t45 * t44;
t57 = t21 * pkin(2) + t61;
t40 = sin(qJ(3));
t17 = t26 * t71 + t68 * t40;
t54 = t17 * pkin(3) + t58;
t12 = t21 * t71 + t40 * t70;
t53 = t12 * pkin(3) + t57;
t52 = t42 * pkin(1) - pkin(8) * t69 + 0;
t19 = t45 * t25 + t42 * t44;
t51 = t19 * pkin(2) + t52;
t10 = t19 * t71 - t40 * t69;
t50 = t10 * pkin(3) + t51;
t16 = t26 * t40 - t68 * t71;
t49 = t16 * qJ(4) + t54;
t11 = t21 * t40 - t42 * t64;
t48 = t11 * qJ(4) + t53;
t9 = t19 * t40 + t45 * t64;
t46 = t9 * qJ(4) + t50;
t38 = -qJ(6) - pkin(10);
t6 = t16 * t39 - t24 * t43;
t5 = t16 * t43 + t24 * t39;
t4 = t11 * t39 + t20 * t43;
t3 = t11 * t43 - t20 * t39;
t2 = t18 * t43 + t9 * t39;
t1 = -t18 * t39 + t9 * t43;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t45, -t42, 0, 0; t42, t45, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t21, -t20, t70, t61; t19, -t18, -t69, t52; t26, t24, t68, t63; 0, 0, 0, 1; t12, -t11, t20, t57 + t74; t10, -t9, t18, t51 + t75; t17, -t16, -t24, t58 - t73; 0, 0, 0, 1; t20, -t12, t11, t48 + t74; t18, -t10, t9, t46 + t75; -t24, -t17, t16, t49 - t73; 0, 0, 0, 1; t4, t3, t12, t12 * pkin(10) + t76 * t20 + t48; t2, t1, t10, t10 * pkin(10) + t76 * t18 + t46; t6, t5, t17, t17 * pkin(10) - t76 * t24 + t49; 0, 0, 0, 1; t4, t3, t12, t62 * t11 - t12 * t38 + t72 * t20 + t53; t2, t1, t10, -t10 * t38 + t72 * t18 + t62 * t9 + t50; t6, t5, t17, t62 * t16 - t17 * t38 - t72 * t24 + t54; 0, 0, 0, 1;];
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
