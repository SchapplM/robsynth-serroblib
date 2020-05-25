% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2018-11-23 16:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPRR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:09:58
% EndTime: 2018-11-23 16:09:58
% DurationCPUTime: 0.14s
% Computational Cost: add. (104->58), mult. (122->57), div. (0->0), fcn. (172->8), ass. (0->29)
t16 = sin(qJ(5));
t36 = pkin(5) * t16;
t17 = sin(qJ(3));
t35 = t17 * t16;
t18 = sin(qJ(1));
t3 = t18 * t17;
t20 = cos(qJ(3));
t34 = t18 * t20;
t21 = cos(qJ(1));
t33 = t21 * t17;
t32 = t21 * t20;
t31 = qJ(4) * t20;
t14 = pkin(6) + 0;
t30 = t18 * pkin(1) + 0;
t29 = pkin(2) + t14;
t28 = t21 * pkin(1) + t18 * qJ(2) + 0;
t9 = t18 * pkin(7);
t27 = t21 * t31 + t30 + t9;
t26 = t21 * pkin(7) + t28;
t25 = t20 * pkin(3) + t17 * qJ(4) + t29;
t24 = pkin(3) * t3 + t26;
t23 = -t21 * qJ(2) + t30;
t22 = -pkin(9) - pkin(8);
t19 = cos(qJ(5));
t15 = qJ(5) + qJ(6);
t6 = cos(t15);
t5 = sin(t15);
t4 = t19 * pkin(5) + pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t18, 0, 0; t18, t21, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; 0, -t21, t18, t28; 0, -t18, -t21, t23; 1, 0, 0, t14; 0, 0, 0, 1; t3, t34, t21, t26; -t33, -t32, t18, t23 + t9; t20, -t17, 0, t29; 0, 0, 0, 1; t21, -t3, -t34, -t18 * t31 + t24; t18, t33, t32 (-pkin(3) * t17 - qJ(2)) * t21 + t27; 0, -t20, t17, t25; 0, 0, 0, 1; -t16 * t34 + t21 * t19, -t21 * t16 - t19 * t34, t3, t21 * pkin(4) + (pkin(8) * t17 - t31) * t18 + t24; t16 * t32 + t18 * t19, -t18 * t16 + t19 * t32, -t33, t18 * pkin(4) + (-qJ(2) + (-pkin(3) - pkin(8)) * t17) * t21 + t27; t35, t17 * t19, t20, t20 * pkin(8) + t25; 0, 0, 0, 1; t21 * t6 - t5 * t34, -t21 * t5 - t6 * t34, t3, t21 * t4 + (-t17 * t22 + (-qJ(4) - t36) * t20) * t18 + t24; t18 * t6 + t5 * t32, -t18 * t5 + t6 * t32, -t33, t18 * t4 + (t20 * t36 - qJ(2) + (-pkin(3) + t22) * t17) * t21 + t27; t17 * t5, t17 * t6, t20, pkin(5) * t35 - t20 * t22 + t25; 0, 0, 0, 1;];
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
