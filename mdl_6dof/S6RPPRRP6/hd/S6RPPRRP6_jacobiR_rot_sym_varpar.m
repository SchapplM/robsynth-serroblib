% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
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
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; -t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t5, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(4));
	t11 = sin(qJ(1));
	t16 = t11 * t10;
	t12 = cos(qJ(4));
	t15 = t11 * t12;
	t13 = cos(qJ(1));
	t14 = t13 * t10;
	t9 = t13 * t12;
	t1 = [-t16, 0, 0, t9, 0, 0; t14, 0, 0, t15, 0, 0; 0, 0, 0, -t10, 0, 0; -t15, 0, 0, -t14, 0, 0; t9, 0, 0, -t16, 0, 0; 0, 0, 0, -t12, 0, 0; -t13, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t65 = sin(qJ(4));
	t66 = sin(qJ(1));
	t76 = t66 * t65;
	t67 = cos(qJ(5));
	t75 = t66 * t67;
	t64 = sin(qJ(5));
	t68 = cos(qJ(4));
	t74 = t68 * t64;
	t73 = t68 * t67;
	t69 = cos(qJ(1));
	t72 = t69 * t65;
	t71 = t69 * t67;
	t70 = t69 * t68;
	t63 = -t66 * t64 + t65 * t71;
	t62 = -t64 * t72 - t75;
	t61 = -t69 * t64 - t65 * t75;
	t60 = t64 * t76 - t71;
	t1 = [t61, 0, 0, t67 * t70, t62, 0; t63, 0, 0, t66 * t73, -t60, 0; 0, 0, 0, -t65 * t67, -t74, 0; t60, 0, 0, -t64 * t70, -t63, 0; t62, 0, 0, -t66 * t74, t61, 0; 0, 0, 0, t65 * t64, -t73, 0; t66 * t68, 0, 0, t72, 0, 0; -t70, 0, 0, t76, 0, 0; 0, 0, 0, t68, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:18
	% EndTime: 2019-10-09 23:56:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t77 = sin(qJ(4));
	t78 = sin(qJ(1));
	t88 = t78 * t77;
	t79 = cos(qJ(5));
	t87 = t78 * t79;
	t76 = sin(qJ(5));
	t80 = cos(qJ(4));
	t86 = t80 * t76;
	t85 = t80 * t79;
	t81 = cos(qJ(1));
	t84 = t81 * t77;
	t83 = t81 * t79;
	t82 = t81 * t80;
	t75 = -t78 * t76 + t77 * t83;
	t74 = t76 * t84 + t87;
	t73 = t81 * t76 + t77 * t87;
	t72 = -t76 * t88 + t83;
	t1 = [-t73, 0, 0, t79 * t82, -t74, 0; t75, 0, 0, t78 * t85, t72, 0; 0, 0, 0, -t77 * t79, -t86, 0; t78 * t80, 0, 0, t84, 0, 0; -t82, 0, 0, t88, 0, 0; 0, 0, 0, t80, 0, 0; t72, 0, 0, t76 * t82, t75, 0; t74, 0, 0, t78 * t86, t73, 0; 0, 0, 0, -t77 * t76, t85, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end