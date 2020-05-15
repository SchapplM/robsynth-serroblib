% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRPPR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:34
% EndTime: 2019-12-31 19:40:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (71->39), mult. (98->36), div. (0->0), fcn. (139->6), ass. (0->25)
t19 = sin(qJ(2));
t20 = sin(qJ(1));
t7 = t20 * t19;
t21 = cos(qJ(5));
t35 = t20 * t21;
t22 = cos(qJ(2));
t8 = t20 * t22;
t23 = cos(qJ(1));
t9 = t23 * t19;
t34 = t23 * t21;
t10 = t23 * t22;
t33 = qJ(3) * t19;
t17 = pkin(5) + 0;
t32 = t19 * pkin(2) + t17;
t31 = t23 * pkin(1) + t20 * pkin(6) + 0;
t30 = pkin(4) * t19 + pkin(7) * t22;
t29 = t20 * pkin(1) - t23 * pkin(6) + 0;
t28 = pkin(2) * t10 + t23 * t33 + t31;
t27 = -t22 * qJ(3) + t32;
t26 = pkin(2) * t8 + t20 * t33 + t29;
t25 = pkin(3) * t8 + t23 * qJ(4) + t26;
t24 = pkin(3) * t10 - t20 * qJ(4) + t28;
t18 = sin(qJ(5));
t12 = t19 * pkin(3);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t20, 0, 0; t20, t23, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; t10, -t9, t20, t31; t8, -t7, -t23, t29; t19, t22, 0, t17; 0, 0, 0, 1; t10, t20, t9, t28; t8, -t23, t7, t26; t19, 0, -t22, t27; 0, 0, 0, 1; t9, -t10, -t20, t24; t7, -t8, t23, t25; -t22, -t19, 0, t12 + t27; 0, 0, 0, 1; -t20 * t18 + t19 * t34, -t18 * t9 - t35, t10, t30 * t23 + t24; t23 * t18 + t19 * t35, -t18 * t7 + t34, t8, t30 * t20 + t25; -t22 * t21, t22 * t18, t19, t19 * pkin(7) + t12 + (-pkin(4) - qJ(3)) * t22 + t32; 0, 0, 0, 1;];
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
