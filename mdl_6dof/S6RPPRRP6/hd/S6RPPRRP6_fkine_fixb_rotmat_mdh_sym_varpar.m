% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2018-11-23 15:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:46:16
% EndTime: 2018-11-23 15:46:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (89->41), mult. (110->35), div. (0->0), fcn. (160->6), ass. (0->31)
t19 = sin(qJ(4));
t20 = sin(qJ(1));
t41 = t20 * t19;
t21 = cos(qJ(5));
t40 = t20 * t21;
t22 = cos(qJ(4));
t39 = t20 * t22;
t18 = sin(qJ(5));
t38 = t22 * t18;
t23 = cos(qJ(1));
t37 = t23 * t19;
t36 = t23 * t21;
t35 = t23 * t22;
t17 = pkin(6) + 0;
t34 = pkin(2) + t17;
t33 = t23 * pkin(1) + t20 * qJ(2) + 0;
t32 = pkin(3) + t34;
t31 = t23 * qJ(3) + t33;
t30 = t20 * pkin(1) - t23 * qJ(2) + 0;
t29 = t22 * pkin(4) + t19 * pkin(8) + t32;
t28 = t20 * qJ(3) + t30;
t27 = -t20 * pkin(7) + t31;
t26 = t23 * pkin(7) + t28;
t25 = pkin(4) * t37 - pkin(8) * t35 + t27;
t24 = pkin(4) * t41 - pkin(8) * t39 + t26;
t7 = t22 * t21;
t4 = -t20 * t18 + t19 * t36;
t3 = t18 * t37 + t40;
t2 = t23 * t18 + t19 * t40;
t1 = t18 * t41 - t36;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t20, 0, 0; t20, t23, 0, 0; 0, 0, 1, t17; 0, 0, 0, 1; 0, -t23, t20, t33; 0, -t20, -t23, t30; 1, 0, 0, t17; 0, 0, 0, 1; 0, t20, t23, t31; 0, -t23, t20, t28; 1, 0, 0, t34; 0, 0, 0, 1; t37, t35, -t20, t27; t41, t39, t23, t26; t22, -t19, 0, t32; 0, 0, 0, 1; t4, -t3, -t35, t25; t2, -t1, -t39, t24; t7, -t38, t19, t29; 0, 0, 0, 1; t4, -t35, t3, t4 * pkin(5) + t3 * qJ(6) + t25; t2, -t39, t1, t2 * pkin(5) + t1 * qJ(6) + t24; t7, t19, t38 (pkin(5) * t21 + qJ(6) * t18) * t22 + t29; 0, 0, 0, 1;];
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
