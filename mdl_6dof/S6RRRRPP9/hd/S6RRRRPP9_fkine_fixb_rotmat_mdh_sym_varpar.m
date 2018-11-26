% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPP9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:10:35
% EndTime: 2018-11-23 18:10:35
% DurationCPUTime: 0.18s
% Computational Cost: add. (690->70), mult. (797->70), div. (0->0), fcn. (903->14), ass. (0->57)
t66 = pkin(6) + qJ(2);
t77 = sin(t66) / 0.2e1;
t76 = pkin(5) + pkin(10);
t67 = pkin(6) - qJ(2);
t60 = sin(t67);
t31 = t77 - t60 / 0.2e1;
t45 = sin(qJ(1));
t46 = cos(qJ(2));
t47 = cos(qJ(1));
t24 = t47 * t31 + t45 * t46;
t43 = sin(qJ(3));
t41 = sin(pkin(6));
t72 = cos(qJ(3));
t64 = t41 * t72;
t12 = t24 * t43 + t47 * t64;
t75 = t12 * pkin(10);
t26 = -t45 * t31 + t47 * t46;
t14 = t26 * t43 - t45 * t64;
t74 = t14 * pkin(10);
t59 = cos(t67) / 0.2e1;
t61 = cos(t66);
t32 = t59 - t61 / 0.2e1;
t68 = cos(pkin(6));
t21 = t32 * t43 - t68 * t72;
t73 = t21 * pkin(10);
t71 = cos(qJ(4));
t70 = t45 * t41;
t69 = t47 * t41;
t65 = pkin(7) + 0;
t63 = t68 * pkin(8) + t65;
t62 = t47 * pkin(1) + pkin(8) * t70 + 0;
t58 = t45 * pkin(1) - pkin(8) * t69 + 0;
t30 = t77 + t60 / 0.2e1;
t57 = t32 * pkin(2) - t30 * pkin(9) + t63;
t44 = sin(qJ(2));
t52 = t59 + t61 / 0.2e1;
t25 = t47 * t44 + t45 * t52;
t56 = t26 * pkin(2) + t25 * pkin(9) + t62;
t22 = t32 * t72 + t68 * t43;
t55 = t22 * pkin(3) + t57;
t15 = t26 * t72 + t43 * t70;
t54 = t15 * pkin(3) + t56;
t23 = t45 * t44 - t47 * t52;
t53 = t24 * pkin(2) + t23 * pkin(9) + t58;
t13 = t24 * t72 - t43 * t69;
t51 = t13 * pkin(3) + t53;
t42 = sin(qJ(4));
t8 = t22 * t42 + t30 * t71;
t9 = t22 * t71 - t30 * t42;
t50 = t9 * pkin(4) + t8 * qJ(5) + t55;
t5 = t15 * t42 - t25 * t71;
t6 = t15 * t71 + t25 * t42;
t49 = t6 * pkin(4) + t5 * qJ(5) + t54;
t3 = t13 * t42 - t23 * t71;
t4 = t13 * t71 + t23 * t42;
t48 = t4 * pkin(4) + t3 * qJ(5) + t51;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t47, -t45, 0, 0; t45, t47, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t26, -t25, t70, t62; t24, -t23, -t69, t58; t32, t30, t68, t63; 0, 0, 0, 1; t15, -t14, t25, t56; t13, -t12, t23, t53; t22, -t21, -t30, t57; 0, 0, 0, 1; t6, -t5, t14, t54 + t74; t4, -t3, t12, t51 + t75; t9, -t8, t21, t55 + t73; 0, 0, 0, 1; t14, -t6, t5, t49 + t74; t12, -t4, t3, t48 + t75; t21, -t9, t8, t50 + t73; 0, 0, 0, 1; t14, t5, t6, t6 * qJ(6) + t76 * t14 + t49; t12, t3, t4, t4 * qJ(6) + t76 * t12 + t48; t21, t8, t9, t9 * qJ(6) + t76 * t21 + t50; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
