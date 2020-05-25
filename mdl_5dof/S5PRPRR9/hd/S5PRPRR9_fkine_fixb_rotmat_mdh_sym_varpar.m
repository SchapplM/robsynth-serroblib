% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRPRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:32
% EndTime: 2019-12-31 17:39:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (112->29), mult. (60->20), div. (0->0), fcn. (98->8), ass. (0->21)
t30 = cos(qJ(4));
t29 = sin(qJ(4));
t28 = pkin(8) + qJ(2);
t16 = sin(pkin(8));
t27 = t16 * pkin(1) + 0;
t17 = cos(pkin(8));
t26 = t17 * pkin(1) + 0;
t25 = qJ(1) + 0;
t13 = pkin(5) + t25;
t24 = sin(t28);
t12 = cos(t28);
t23 = t12 * pkin(2) + t24 * qJ(3) + t26;
t22 = t12 * pkin(3) + t23;
t21 = t24 * pkin(2) - t12 * qJ(3) + t27;
t20 = t24 * pkin(3) + t21;
t19 = cos(qJ(5));
t18 = sin(qJ(5));
t11 = -pkin(6) + t13;
t2 = t12 * t29 - t24 * t30;
t1 = -t12 * t30 - t24 * t29;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t16, 0, 0; t16, t17, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t12, -t24, 0, t26; t24, t12, 0, t27; 0, 0, 1, t13; 0, 0, 0, 1; t12, 0, t24, t23; t24, 0, -t12, t21; 0, 1, 0, t13; 0, 0, 0, 1; -t1, -t2, 0, t22; -t2, t1, 0, t20; 0, 0, -1, t11; 0, 0, 0, 1; -t1 * t19, t1 * t18, t2, -t1 * pkin(4) + t2 * pkin(7) + t22; -t2 * t19, t2 * t18, -t1, -t2 * pkin(4) - t1 * pkin(7) + t20; -t18, -t19, 0, t11; 0, 0, 0, 1;];
T_ges = t3;
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
