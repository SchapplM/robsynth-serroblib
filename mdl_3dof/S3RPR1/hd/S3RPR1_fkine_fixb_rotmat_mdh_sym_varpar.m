% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% T_c_mdh [4x4x(3+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   4:  mdh base (link 0) -> mdh frame (4-1), link (4-1)
%   ...
%   3+1:  mdh base (link 0) -> mdh frame (3)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S3RPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:21
% EndTime: 2018-11-14 10:14:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (20->13), mult. (18->10), div. (0->0), fcn. (34->4), ass. (0->10)
t6 = pkin(3) + 0;
t10 = cos(qJ(1));
t8 = sin(qJ(1));
t12 = t10 * pkin(1) + t8 * qJ(2) + 0;
t11 = t8 * pkin(1) - t10 * qJ(2) + 0;
t9 = cos(qJ(3));
t7 = sin(qJ(3));
t2 = -t10 * t7 + t8 * t9;
t1 = -t10 * t9 - t8 * t7;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t10, -t8, 0, 0; t8, t10, 0, 0; 0, 0, 1, t6; 0, 0, 0, 1; t10, 0, t8, t12; t8, 0, -t10, t11; 0, 1, 0, t6; 0, 0, 0, 1; -t1, t2, 0, t10 * pkin(2) + t12; t2, t1, 0, t8 * pkin(2) + t11; 0, 0, -1, -pkin(4) + t6; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,3+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,3+1]); end % symbolisch
for i = 1:3+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
