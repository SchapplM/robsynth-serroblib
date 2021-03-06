% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S4PRRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(7));
	t6 = sin(pkin(7));
	t1 = [0, -t7 * t8, 0, 0; 0, -t6 * t8, 0, 0; 0, t9, 0, 0; 0, -t7 * t9, 0, 0; 0, -t6 * t9, 0, 0; 0, -t8, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
	t42 = sin(qJ(3));
	t43 = sin(qJ(2));
	t49 = t43 * t42;
	t44 = cos(qJ(3));
	t48 = t43 * t44;
	t45 = cos(qJ(2));
	t47 = t45 * t42;
	t46 = t45 * t44;
	t41 = cos(pkin(7));
	t40 = sin(pkin(7));
	t1 = [0, -t41 * t48, t40 * t44 - t41 * t47, 0; 0, -t40 * t48, -t40 * t47 - t41 * t44, 0; 0, t46, -t49, 0; 0, t41 * t49, -t40 * t42 - t41 * t46, 0; 0, t40 * t49, -t40 * t46 + t41 * t42, 0; 0, -t47, -t48, 0; 0, t41 * t45, 0, 0; 0, t40 * t45, 0, 0; 0, t43, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (41->12), mult. (40->18), div. (0->0), fcn. (69->6), ass. (0->16)
	t64 = qJ(3) + qJ(4);
	t62 = sin(t64);
	t67 = sin(qJ(2));
	t72 = t67 * t62;
	t63 = cos(t64);
	t71 = t67 * t63;
	t68 = cos(qJ(2));
	t70 = t68 * t62;
	t69 = t68 * t63;
	t66 = cos(pkin(7));
	t65 = sin(pkin(7));
	t61 = -t65 * t62 - t66 * t69;
	t60 = t65 * t63 - t66 * t70;
	t59 = t66 * t62 - t65 * t69;
	t58 = -t66 * t63 - t65 * t70;
	t1 = [0, -t66 * t71, t60, t60; 0, -t65 * t71, t58, t58; 0, t69, -t72, -t72; 0, t66 * t72, t61, t61; 0, t65 * t72, t59, t59; 0, -t70, -t71, -t71; 0, t66 * t68, 0, 0; 0, t65 * t68, 0, 0; 0, t67, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,4);
end