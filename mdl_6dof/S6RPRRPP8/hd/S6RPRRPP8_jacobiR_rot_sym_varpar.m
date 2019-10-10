% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
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
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
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
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (16->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t65 = sin(qJ(3));
	t66 = sin(qJ(1));
	t76 = t66 * t65;
	t67 = cos(qJ(4));
	t75 = t66 * t67;
	t64 = sin(qJ(4));
	t68 = cos(qJ(3));
	t74 = t68 * t64;
	t73 = t68 * t67;
	t69 = cos(qJ(1));
	t72 = t69 * t65;
	t71 = t69 * t67;
	t70 = t69 * t68;
	t63 = -t66 * t64 + t65 * t71;
	t62 = t64 * t72 + t75;
	t61 = t69 * t64 + t65 * t75;
	t60 = -t64 * t76 + t71;
	t1 = [t63, 0, t66 * t73, t60, 0, 0; t61, 0, -t67 * t70, t62, 0, 0; 0, 0, -t65 * t67, -t74, 0, 0; -t62, 0, -t66 * t74, -t61, 0, 0; t60, 0, t64 * t70, t63, 0, 0; 0, 0, t65 * t64, -t73, 0, 0; -t70, 0, t76, 0, 0, 0; -t66 * t68, 0, -t72, 0, 0, 0; 0, 0, t68, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t75 = sin(qJ(3));
	t76 = sin(qJ(1));
	t86 = t76 * t75;
	t77 = cos(qJ(4));
	t85 = t76 * t77;
	t74 = sin(qJ(4));
	t78 = cos(qJ(3));
	t84 = t78 * t74;
	t83 = t78 * t77;
	t79 = cos(qJ(1));
	t82 = t79 * t75;
	t81 = t79 * t77;
	t80 = t79 * t78;
	t73 = t76 * t74 - t75 * t81;
	t72 = t74 * t82 + t85;
	t71 = t79 * t74 + t75 * t85;
	t70 = t74 * t86 - t81;
	t1 = [-t80, 0, t86, 0, 0, 0; -t76 * t78, 0, -t82, 0, 0, 0; 0, 0, t78, 0, 0, 0; t73, 0, -t76 * t83, t70, 0, 0; -t71, 0, t77 * t80, -t72, 0, 0; 0, 0, t75 * t77, t84, 0, 0; t72, 0, t76 * t84, t71, 0, 0; t70, 0, -t74 * t80, t73, 0, 0; 0, 0, -t75 * t74, t83, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (16->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t75 = sin(qJ(3));
	t76 = sin(qJ(1));
	t86 = t76 * t75;
	t77 = cos(qJ(4));
	t85 = t76 * t77;
	t74 = sin(qJ(4));
	t78 = cos(qJ(3));
	t84 = t78 * t74;
	t83 = t78 * t77;
	t79 = cos(qJ(1));
	t82 = t79 * t75;
	t81 = t79 * t77;
	t80 = t79 * t78;
	t73 = -t76 * t74 + t75 * t81;
	t72 = t74 * t82 + t85;
	t71 = t79 * t74 + t75 * t85;
	t70 = t74 * t86 - t81;
	t1 = [-t80, 0, t86, 0, 0, 0; -t76 * t78, 0, -t82, 0, 0, 0; 0, 0, t78, 0, 0, 0; t72, 0, t76 * t84, t71, 0, 0; t70, 0, -t74 * t80, -t73, 0, 0; 0, 0, -t75 * t74, t83, 0, 0; t73, 0, t76 * t83, -t70, 0, 0; t71, 0, -t77 * t80, t72, 0, 0; 0, 0, -t75 * t77, -t84, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end