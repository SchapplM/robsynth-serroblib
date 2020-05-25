% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:36
% EndTime: 2019-12-31 16:19:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (41->28), mult. (46->27), div. (0->0), fcn. (75->6), ass. (0->20)
t10 = sin(qJ(3));
t7 = sin(pkin(6));
t23 = t7 * t10;
t12 = cos(qJ(3));
t22 = t7 * t12;
t8 = cos(pkin(6));
t21 = t8 * t10;
t20 = t8 * t12;
t11 = cos(qJ(4));
t19 = t10 * t11;
t18 = t7 * pkin(1) + 0;
t6 = qJ(1) + 0;
t17 = t8 * pkin(1) + t7 * qJ(2) + 0;
t16 = pkin(2) + t6;
t15 = t8 * pkin(4) + t17;
t14 = pkin(3) * t10 - pkin(5) * t12;
t13 = -t8 * qJ(2) + t18;
t9 = sin(qJ(4));
t2 = t7 * pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t8, -t7, 0, 0; t7, t8, 0, 0; 0, 0, 1, t6; 0, 0, 0, 1; 0, -t8, t7, t17; 0, -t7, -t8, t13; 1, 0, 0, t6; 0, 0, 0, 1; t23, t22, t8, t15; -t21, -t20, t7, t13 + t2; t12, -t10, 0, t16; 0, 0, 0, 1; t7 * t19 + t8 * t9, t8 * t11 - t9 * t23, -t22, t14 * t7 + t15; -t8 * t19 + t7 * t9, t7 * t11 + t9 * t21, t20, t2 + (-qJ(2) - t14) * t8 + t18; t12 * t11, -t12 * t9, t10, t12 * pkin(3) + t10 * pkin(5) + t16; 0, 0, 0, 1;];
T_ges = t1;
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
