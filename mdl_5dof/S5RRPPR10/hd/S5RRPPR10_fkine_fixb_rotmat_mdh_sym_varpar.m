% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPPR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:42:54
% EndTime: 2019-12-31 19:42:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (98->46), mult. (179->51), div. (0->0), fcn. (250->8), ass. (0->30)
t23 = cos(pkin(8));
t25 = sin(qJ(2));
t13 = t25 * t23;
t22 = sin(pkin(8));
t41 = t25 * t22;
t42 = pkin(3) * t13 + qJ(4) * t41;
t26 = sin(qJ(1));
t15 = t26 * t25;
t28 = cos(qJ(2));
t40 = t26 * t28;
t29 = cos(qJ(1));
t16 = t29 * t25;
t39 = t29 * t28;
t38 = qJ(3) * t25;
t21 = pkin(5) + 0;
t37 = t25 * pkin(2) + t21;
t36 = t29 * pkin(1) + t26 * pkin(6) + 0;
t35 = t26 * pkin(1) - t29 * pkin(6) + 0;
t34 = pkin(2) * t39 + t29 * t38 + t36;
t33 = -t28 * qJ(3) + t37;
t32 = pkin(2) * t40 + t26 * t38 + t35;
t5 = t22 * t39 - t26 * t23;
t6 = t26 * t22 + t23 * t39;
t31 = t6 * pkin(3) + t5 * qJ(4) + t34;
t3 = t22 * t40 + t29 * t23;
t4 = -t29 * t22 + t23 * t40;
t30 = t4 * pkin(3) + t3 * qJ(4) + t32;
t27 = cos(qJ(5));
t24 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t39, -t16, t26, t36; t40, -t15, -t29, t35; t25, t28, 0, t21; 0, 0, 0, 1; t6, -t5, t16, t34; t4, -t3, t15, t32; t13, -t41, -t28, t33; 0, 0, 0, 1; t6, t16, t5, t31; t4, t15, t3, t30; t13, -t28, t41, t33 + t42; 0, 0, 0, 1; t5 * t24 + t6 * t27, -t6 * t24 + t5 * t27, -t16, t6 * pkin(4) - pkin(7) * t16 + t31; t3 * t24 + t4 * t27, -t4 * t24 + t3 * t27, -t15, t4 * pkin(4) - pkin(7) * t15 + t30; (t22 * t24 + t23 * t27) * t25, (t22 * t27 - t23 * t24) * t25, t28, pkin(4) * t13 + (pkin(7) - qJ(3)) * t28 + t37 + t42; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
