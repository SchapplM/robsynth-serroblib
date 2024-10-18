% Calculate homogenous joint transformation matrices for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PRRRR12_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:06:58
% EndTime: 2024-09-28 18:06:58
% DurationCPUTime: 0.00s
% Computational Cost: add. (9->9), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t97 = cos(qJ(2));
t96 = cos(qJ(3));
t95 = cos(qJ(4));
t94 = cos(qJ(5));
t93 = sin(qJ(2));
t92 = sin(qJ(3));
t91 = sin(qJ(4));
t90 = sin(qJ(5));
t89 = cos(pkin(5));
t88 = cos(pkin(6));
t87 = cos(pkin(11));
t86 = sin(pkin(5));
t85 = sin(pkin(6));
t84 = sin(pkin(11));
t1 = [t87, -t84, 0, 0; t84, t87, 0, 0; 0, 0, 1, qJ(1); t97, -t93, 0, pkin(1); t89 * t93, t89 * t97, -t86, -t86 * pkin(7); t86 * t93, t86 * t97, t89, t89 * pkin(7); t96, -t92, 0, pkin(2); t92, t96, 0, 0; 0, 0, 1, pkin(8); t95, -t91, 0, pkin(3); t91, t95, 0, 0; 0, 0, 1, pkin(9); t94, -t90, 0, pkin(4); t88 * t90, t88 * t94, -t85, -t85 * pkin(10); t85 * t90, t85 * t94, t88, t88 * pkin(10);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
