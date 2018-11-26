% Calculate homogenous joint transformation matrices for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% T_mdh [4x4x2]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S2RR2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:48:28
% EndTime: 2018-11-16 16:48:28
% DurationCPUTime: 0.02s
% Computational Cost: add. (5->5), mult. (0->0), div. (0->0), fcn. (8->4), ass. (0->5)
t8 = cos(qJ(1));
t7 = cos(qJ(2));
t6 = sin(qJ(1));
t5 = sin(qJ(2));
t1 = [t8, -t6, 0, 0; 0, 0, 1, 0; -t6, -t8, 0, 0; 0, 0, 0, 1; t7, -t5, 0, 0; 0, 0, -1, -pkin(1); t5, t7, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,2);             % numerisch
else,                         T_mdh = sym('xx', [4,4,2]); end % symbolisch

for i = 1:2
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
