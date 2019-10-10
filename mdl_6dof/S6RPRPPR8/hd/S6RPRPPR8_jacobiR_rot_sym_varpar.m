% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; -t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0, 0; t16, 0, -t15, 0, 0, 0; 0, 0, -t11, 0, 0, 0; t15, 0, -t16, 0, 0, 0; t9, 0, t10, 0, 0, 0; 0, 0, -t13, 0, 0, 0; -t12, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t39 = sin(qJ(1));
	t40 = cos(qJ(3));
	t44 = t39 * t40;
	t38 = sin(qJ(3));
	t41 = cos(qJ(1));
	t43 = t41 * t38;
	t42 = t41 * t40;
	t37 = t39 * t38;
	t1 = [t43, 0, t44, 0, 0, 0; t37, 0, -t42, 0, 0, 0; 0, 0, -t38, 0, 0, 0; -t39, 0, 0, 0, 0, 0; t41, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t42, 0, t37, 0, 0, 0; -t44, 0, -t43, 0, 0, 0; 0, 0, t40, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t18 = t12 * t11;
	t13 = cos(qJ(3));
	t17 = t12 * t13;
	t14 = cos(qJ(1));
	t16 = t14 * t11;
	t15 = t14 * t13;
	t1 = [-t15, 0, t18, 0, 0, 0; -t17, 0, -t16, 0, 0, 0; 0, 0, t13, 0, 0, 0; -t16, 0, -t17, 0, 0, 0; -t18, 0, t15, 0, 0, 0; 0, 0, t11, 0, 0, 0; t12, 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t64 = sin(qJ(3));
	t65 = sin(qJ(1));
	t75 = t65 * t64;
	t66 = cos(qJ(6));
	t74 = t65 * t66;
	t63 = sin(qJ(6));
	t67 = cos(qJ(3));
	t73 = t67 * t63;
	t72 = t67 * t66;
	t68 = cos(qJ(1));
	t71 = t68 * t64;
	t70 = t68 * t66;
	t69 = t68 * t67;
	t62 = t65 * t63 - t66 * t69;
	t61 = t63 * t69 + t74;
	t60 = t68 * t63 + t65 * t72;
	t59 = t65 * t73 - t70;
	t1 = [t62, 0, t64 * t74, 0, 0, t59; -t60, 0, -t64 * t70, 0, 0, -t61; 0, 0, t72, 0, 0, -t64 * t63; t61, 0, -t63 * t75, 0, 0, t60; t59, 0, t63 * t71, 0, 0, t62; 0, 0, -t73, 0, 0, -t64 * t66; t71, 0, t65 * t67, 0, 0, 0; t75, 0, -t69, 0, 0, 0; 0, 0, -t64, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end