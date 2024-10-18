% Calculate homogenous joint transformation matrices for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PRRRR11_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:44:59
% EndTime: 2024-09-27 21:44:59
% DurationCPUTime: 0.00s
% Computational Cost: add. (7->7), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t64 = cos(qJ(2));
t63 = cos(qJ(3));
t62 = cos(qJ(4));
t61 = cos(qJ(5));
t60 = sin(qJ(2));
t59 = sin(qJ(3));
t58 = sin(qJ(4));
t57 = sin(qJ(5));
t56 = cos(pkin(5));
t55 = cos(pkin(10));
t54 = sin(pkin(5));
t53 = sin(pkin(10));
t1 = [t55, -t53, 0, 0; t53, t55, 0, 0; 0, 0, 1, qJ(1); t64, -t60, 0, pkin(1); t60, t64, 0, 0; 0, 0, 1, pkin(6); t63, -t59, 0, pkin(2); t56 * t59, t56 * t63, -t54, -t54 * pkin(7); t54 * t59, t54 * t63, t56, t56 * pkin(7); t62, -t58, 0, pkin(3); t58, t62, 0, 0; 0, 0, 1, pkin(8); t61, -t57, 0, pkin(4); t57, t61, 0, 0; 0, 0, 1, pkin(9);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
