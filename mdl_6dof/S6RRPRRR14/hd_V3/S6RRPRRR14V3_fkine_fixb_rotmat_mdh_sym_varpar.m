% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRRR14V3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:16
% EndTime: 2019-04-12 15:03:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (61->32), mult. (152->42), div. (0->0), fcn. (240->10), ass. (0->31)
t20 = sin(qJ(4));
t21 = sin(qJ(2));
t30 = t21 * t20;
t25 = cos(qJ(4));
t29 = t21 * t25;
t22 = sin(qJ(1));
t14 = t22 * t21;
t26 = cos(qJ(2));
t15 = t22 * t26;
t27 = cos(qJ(1));
t16 = t27 * t21;
t17 = t27 * t26;
t28 = qJ(3) * t21;
t24 = cos(qJ(5));
t23 = cos(qJ(6));
t19 = sin(qJ(5));
t18 = sin(qJ(6));
t13 = -t26 * qJ(3) + 0;
t12 = t27 * t28 + 0;
t11 = t22 * t28 + 0;
t10 = t25 * t17 + t22 * t20;
t9 = t20 * t17 - t22 * t25;
t8 = t25 * t15 - t27 * t20;
t7 = t20 * t15 + t27 * t25;
t6 = -t26 * t19 + t24 * t29;
t5 = t19 * t29 + t26 * t24;
t4 = t10 * t24 + t19 * t16;
t3 = t10 * t19 - t24 * t16;
t2 = t19 * t14 + t8 * t24;
t1 = -t24 * t14 + t8 * t19;
t31 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t22, 0, 0; t22, t27, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t16, t22, 0; t15, -t14, -t27, 0; t21, t26, 0, 0; 0, 0, 0, 1; t17, t22, t16, t12; t15, -t27, t14, t11; t21, 0, -t26, t13; 0, 0, 0, 1; t10, -t9, t16, t12; t8, -t7, t14, t11; t29, -t30, -t26, t13; 0, 0, 0, 1; t4, -t3, t9, t12; t2, -t1, t7, t11; t6, -t5, t30, t13; 0, 0, 0, 1; t9 * t18 + t4 * t23, -t4 * t18 + t9 * t23, t3, t12; t7 * t18 + t2 * t23, -t2 * t18 + t7 * t23, t1, t11; t18 * t30 + t6 * t23, -t6 * t18 + t23 * t30, t5, t13; 0, 0, 0, 1;];
T_ges = t31;
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
