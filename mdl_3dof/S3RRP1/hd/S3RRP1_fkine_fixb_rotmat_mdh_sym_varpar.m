% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S3RRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:05
% EndTime: 2018-11-14 10:15:05
% DurationCPUTime: 0.05s
% Computational Cost: add. (28->12), mult. (8->6), div. (0->0), fcn. (20->4), ass. (0->10)
t11 = pkin(3) + 0;
t7 = sin(qJ(1));
t10 = t7 * pkin(1) + 0;
t8 = cos(qJ(1));
t9 = t8 * pkin(1) + 0;
t6 = qJ(1) + qJ(2);
t3 = pkin(4) + t11;
t2 = cos(t6);
t1 = sin(t6);
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t8, -t7, 0, 0; t7, t8, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t2, -t1, 0, t9; t1, t2, 0, t10; 0, 0, 1, t3; 0, 0, 0, 1; t2, 0, t1, t2 * pkin(2) + t1 * qJ(3) + t9; t1, 0, -t2, t1 * pkin(2) - t2 * qJ(3) + t10; 0, 1, 0, t3; 0, 0, 0, 1;];
T_ges = t4;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,3+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,3+1]); end % symbolisch
for i = 1:3+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
