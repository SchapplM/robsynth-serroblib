% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:57
% EndTime: 2019-12-31 16:35:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (99->42), mult. (221->56), div. (0->0), fcn. (311->10), ass. (0->35)
t21 = sin(pkin(8));
t22 = sin(pkin(4));
t44 = t21 * t22;
t27 = sin(qJ(2));
t43 = t22 * t27;
t29 = cos(qJ(3));
t42 = t22 * t29;
t30 = cos(qJ(2));
t41 = t22 * t30;
t23 = cos(pkin(8));
t40 = t23 * t22;
t24 = cos(pkin(4));
t39 = t24 * t27;
t38 = t24 * t30;
t37 = qJ(1) + 0;
t36 = t23 * pkin(1) + pkin(5) * t44 + 0;
t35 = t24 * pkin(5) + t37;
t34 = t21 * pkin(1) - pkin(5) * t40 + 0;
t10 = -t21 * t39 + t23 * t30;
t9 = t21 * t38 + t23 * t27;
t33 = t10 * pkin(2) + t9 * pkin(6) + t36;
t32 = pkin(2) * t43 - pkin(6) * t41 + t35;
t7 = t21 * t27 - t23 * t38;
t8 = t21 * t30 + t23 * t39;
t31 = t8 * pkin(2) + t7 * pkin(6) + t34;
t28 = cos(qJ(4));
t26 = sin(qJ(3));
t25 = sin(qJ(4));
t12 = t24 * t26 + t27 * t42;
t11 = -t24 * t29 + t26 * t43;
t4 = t10 * t29 + t26 * t44;
t3 = t10 * t26 - t21 * t42;
t2 = -t26 * t40 + t8 * t29;
t1 = t8 * t26 + t29 * t40;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t21, 0, 0; t21, t23, 0, 0; 0, 0, 1, t37; 0, 0, 0, 1; t10, -t9, t44, t36; t8, -t7, -t40, t34; t43, t41, t24, t35; 0, 0, 0, 1; t4, -t3, t9, t33; t2, -t1, t7, t31; t12, -t11, -t41, t32; 0, 0, 0, 1; t9 * t25 + t4 * t28, -t4 * t25 + t9 * t28, t3, t4 * pkin(3) + t3 * pkin(7) + t33; t2 * t28 + t7 * t25, -t2 * t25 + t7 * t28, t1, t2 * pkin(3) + t1 * pkin(7) + t31; t12 * t28 - t25 * t41, -t12 * t25 - t28 * t41, t11, t12 * pkin(3) + t11 * pkin(7) + t32; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
