% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRPR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:17
% EndTime: 2019-12-31 21:32:17
% DurationCPUTime: 0.13s
% Computational Cost: add. (98->46), mult. (179->50), div. (0->0), fcn. (250->8), ass. (0->29)
t24 = sin(qJ(2));
t27 = cos(qJ(3));
t13 = t24 * t27;
t23 = sin(qJ(3));
t40 = t24 * t23;
t41 = pkin(3) * t13 + qJ(4) * t40;
t25 = sin(qJ(1));
t14 = t25 * t24;
t28 = cos(qJ(2));
t39 = t25 * t28;
t29 = cos(qJ(1));
t16 = t29 * t24;
t38 = t29 * t28;
t21 = pkin(5) + 0;
t37 = t24 * pkin(2) + t21;
t36 = t29 * pkin(1) + t25 * pkin(6) + 0;
t35 = t25 * pkin(1) - t29 * pkin(6) + 0;
t34 = pkin(2) * t38 + pkin(7) * t16 + t36;
t33 = -t28 * pkin(7) + t37;
t32 = pkin(2) * t39 + pkin(7) * t14 + t35;
t5 = t23 * t38 - t25 * t27;
t6 = t25 * t23 + t27 * t38;
t31 = t6 * pkin(3) + t5 * qJ(4) + t34;
t3 = t23 * t39 + t29 * t27;
t4 = -t29 * t23 + t27 * t39;
t30 = t4 * pkin(3) + t3 * qJ(4) + t32;
t26 = cos(qJ(5));
t22 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t25, 0, 0; t25, t29, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t38, -t16, t25, t36; t39, -t14, -t29, t35; t24, t28, 0, t21; 0, 0, 0, 1; t6, -t5, t16, t34; t4, -t3, t14, t32; t13, -t40, -t28, t33; 0, 0, 0, 1; t6, t16, t5, t31; t4, t14, t3, t30; t13, -t28, t40, t33 + t41; 0, 0, 0, 1; t5 * t22 + t6 * t26, -t6 * t22 + t5 * t26, -t16, t6 * pkin(4) - pkin(8) * t16 + t31; t3 * t22 + t4 * t26, -t4 * t22 + t3 * t26, -t14, t4 * pkin(4) - pkin(8) * t14 + t30; (t22 * t23 + t26 * t27) * t24, (-t22 * t27 + t23 * t26) * t24, t28, pkin(4) * t13 + (-pkin(7) + pkin(8)) * t28 + t37 + t41; 0, 0, 0, 1;];
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
