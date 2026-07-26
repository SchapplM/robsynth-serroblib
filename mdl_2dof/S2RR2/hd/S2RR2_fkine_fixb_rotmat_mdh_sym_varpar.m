% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
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
% Tc_mdh [4x4x(2+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   3:  mdh base (link 0) -> mdh frame (3-1), link (3-1)
%   ...
%   2+1:  mdh base (link 0) -> mdh frame (2)
% T_c_stack [(2+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2026-07-26 14:41
% Revision: 05cfe754a3b58d1958e6f60508538bf5b2996c81 (2026-06-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S2RR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2026-07-26 14:41:22
% EndTime: 2026-07-26 14:41:22
% DurationCPUTime: 0.14s
% Computational Cost: add. (7->7), mult. (6->6), div. (0->0), fcn. (18->4), ass. (0->5)
t4 = cos(qJ(1));
t3 = cos(qJ(2));
t2 = sin(qJ(1));
t1 = sin(qJ(2));
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t4, -t2, 0, 0; 0, 0, 1, 0; -t2, -t4, 0, 0; t4 * t3, -t4 * t1, t2, t2 * pkin(1) + 0; t1, t3, 0, 0; -t2 * t3, t2 * t1, t4, t4 * pkin(1) + 0;];
Tc_stack = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,2+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,2+1]); end % symbolisch
for i = 1:2+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
