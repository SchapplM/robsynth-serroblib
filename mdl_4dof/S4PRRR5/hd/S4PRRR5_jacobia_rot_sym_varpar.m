% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4PRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4PRRR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:51
	% EndTime: 2019-12-29 12:29:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:56
	% EndTime: 2019-12-29 12:29:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:29:56
	% EndTime: 2019-12-29 12:29:56
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (257->15), mult. (243->34), div. (59->8), fcn. (349->9), ass. (0->26)
	t55 = qJ(2) + qJ(3);
	t54 = cos(t55);
	t53 = sin(t55);
	t56 = sin(pkin(7));
	t61 = t56 * t53;
	t49 = atan2(-t61, -t54);
	t47 = sin(t49);
	t64 = t47 * t54;
	t51 = t53 ^ 2;
	t63 = t51 / t54 ^ 2;
	t57 = cos(pkin(7));
	t62 = t54 * t57;
	t58 = sin(qJ(4));
	t59 = cos(qJ(4));
	t46 = t56 * t58 + t59 * t62;
	t44 = 0.1e1 / t46 ^ 2;
	t45 = -t56 * t59 + t58 * t62;
	t60 = t45 ^ 2 * t44 + 0.1e1;
	t48 = cos(t49);
	t43 = 0.1e1 / t60;
	t42 = (0.1e1 + t63) * t56 / (t56 ^ 2 * t63 + 0.1e1);
	t41 = -t47 * t61 - t48 * t54;
	t40 = 0.1e1 / t41 ^ 2;
	t38 = (-t58 / t46 + t59 * t45 * t44) * t57 * t53 * t43;
	t37 = (t54 / t41 - (-t56 * t64 + t48 * t53 + (-t48 * t61 + t64) * t42) * t53 * t40) * t57 / (t57 ^ 2 * t51 * t40 + 0.1e1);
	t1 = [0, t42, t42, 0; 0, t37, t37, 0; 0, t38, t38, t60 * t43;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end