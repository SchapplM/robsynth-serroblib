% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:59
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S4PPRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:59:05
% EndTime: 2018-11-14 13:59:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->14), mult. (6->4), div. (0->0), fcn. (18->6), ass. (0->15)
t16 = pkin(1) + 0;
t9 = pkin(6) + qJ(3);
t8 = qJ(1) + 0;
t15 = qJ(2) + 0;
t11 = cos(pkin(6));
t14 = t11 * pkin(2) + t16;
t13 = pkin(4) + t15;
t10 = sin(pkin(6));
t12 = t10 * pkin(2) + t8;
t5 = qJ(4) + t9;
t4 = cos(t9);
t3 = sin(t9);
t2 = cos(t5);
t1 = sin(t5);
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, 0; 0, 0, 1, t8; 0, -1, 0, 0; 0, 0, 0, 1; t11, -t10, 0, t16; t10, t11, 0, t8; 0, 0, 1, t15; 0, 0, 0, 1; t4, -t3, 0, t14; t3, t4, 0, t12; 0, 0, 1, t13; 0, 0, 0, 1; t2, -t1, 0, pkin(3) * t4 + t14; t1, t2, 0, pkin(3) * t3 + t12; 0, 0, 1, pkin(5) + t13; 0, 0, 0, 1;];
T_ges = t6;
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
