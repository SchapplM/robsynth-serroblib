% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2018-11-23 15:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:41:30
% EndTime: 2018-11-23 15:41:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (98->51), mult. (97->48), div. (0->0), fcn. (143->8), ass. (0->31)
t12 = sin(pkin(9));
t34 = pkin(5) * t12;
t15 = sin(qJ(4));
t16 = sin(qJ(1));
t33 = t16 * t15;
t17 = cos(qJ(4));
t32 = t16 * t17;
t18 = cos(qJ(1));
t31 = t18 * t15;
t30 = t18 * t17;
t11 = pkin(6) + 0;
t29 = t16 * pkin(1) + 0;
t28 = pkin(2) + t11;
t27 = t18 * pkin(1) + t16 * qJ(2) + 0;
t26 = pkin(3) + t28;
t25 = t18 * qJ(3) + t27;
t13 = cos(pkin(9));
t1 = t13 * pkin(5) + pkin(4);
t14 = -pkin(8) - qJ(5);
t24 = t1 * t15 + t14 * t17;
t23 = pkin(4) * t15 - qJ(5) * t17;
t22 = -t18 * qJ(2) + t29;
t4 = t16 * qJ(3);
t21 = t22 + t4;
t20 = -t16 * pkin(7) + t25;
t8 = t18 * pkin(7);
t19 = t21 + t8;
t10 = pkin(9) + qJ(6);
t3 = cos(t10);
t2 = sin(t10);
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t16, 0, 0; t16, t18, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; 0, -t18, t16, t27; 0, -t16, -t18, t22; 1, 0, 0, t11; 0, 0, 0, 1; 0, t16, t18, t25; 0, -t18, t16, t21; 1, 0, 0, t28; 0, 0, 0, 1; t31, t30, -t16, t20; t33, t32, t18, t19; t17, -t15, 0, t26; 0, 0, 0, 1; -t16 * t12 + t13 * t31, -t12 * t31 - t16 * t13, -t30, t23 * t18 + t20; t18 * t12 + t13 * t33, -t12 * t33 + t18 * t13, -t32, t23 * t16 + t19; t17 * t13, -t17 * t12, t15, t17 * pkin(4) + t15 * qJ(5) + t26; 0, 0, 0, 1; -t16 * t2 + t3 * t31, -t16 * t3 - t2 * t31, -t30, t24 * t18 + (-pkin(7) - t34) * t16 + t25; t18 * t2 + t3 * t33, t18 * t3 - t2 * t33, -t32, t4 + t8 + (-qJ(2) + t34) * t18 + t24 * t16 + t29; t17 * t3, -t17 * t2, t15, t17 * t1 - t15 * t14 + t26; 0, 0, 0, 1;];
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
