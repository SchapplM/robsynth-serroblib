% Calculate homogenous joint transformation matrices for
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
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S6RRPRRR14V3_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:16
% EndTime: 2019-04-12 15:03:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (6->6), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t40 = cos(qJ(1));
t39 = cos(qJ(2));
t38 = cos(qJ(4));
t37 = cos(qJ(5));
t36 = cos(qJ(6));
t35 = sin(qJ(1));
t34 = sin(qJ(2));
t33 = sin(qJ(4));
t32 = sin(qJ(5));
t31 = sin(qJ(6));
t1 = [t40, -t35, 0, 0; t35, t40, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t39, -t34, 0, 0; 0, 0, -1, 0; t34, t39, 0, 0; 0, 0, 0, 1; 1, 0, 0, 0; 0, 0, -1, -qJ(3); 0, 1, 0, 0; 0, 0, 0, 1; t38, -t33, 0, 0; t33, t38, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t37, -t32, 0, 0; 0, 0, -1, 0; t32, t37, 0, 0; 0, 0, 0, 1; t36, -t31, 0, 0; 0, 0, -1, 0; t31, t36, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
