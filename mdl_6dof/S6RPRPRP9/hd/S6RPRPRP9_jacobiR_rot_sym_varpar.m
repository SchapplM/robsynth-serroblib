% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
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
	% StartTime: 2019-10-10 00:42:47
	% EndTime: 2019-10-10 00:42:48
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
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
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
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
	t50 = sin(qJ(3));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(3));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(9));
	t48 = sin(pkin(9));
	t1 = [-t51 * t48 + t49 * t55, 0, t49 * t56, 0, 0, 0; t53 * t48 + t49 * t57, 0, -t49 * t54, 0, 0, 0; 0, 0, -t50 * t49, 0, 0, 0; -t48 * t55 - t51 * t49, 0, -t48 * t56, 0, 0, 0; -t48 * t57 + t53 * t49, 0, t48 * t54, 0, 0, 0; 0, 0, t50 * t48, 0, 0, 0; -t54, 0, t57, 0, 0, 0; -t56, 0, -t55, 0, 0, 0; 0, 0, t52, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (40->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t70 = sin(qJ(3));
	t71 = sin(qJ(1));
	t78 = t71 * t70;
	t69 = pkin(9) + qJ(5);
	t67 = sin(t69);
	t72 = cos(qJ(3));
	t77 = t72 * t67;
	t68 = cos(t69);
	t76 = t72 * t68;
	t73 = cos(qJ(1));
	t75 = t73 * t70;
	t74 = t73 * t72;
	t66 = -t71 * t67 + t68 * t75;
	t65 = t67 * t75 + t71 * t68;
	t64 = t73 * t67 + t68 * t78;
	t63 = -t67 * t78 + t73 * t68;
	t1 = [t66, 0, t71 * t76, 0, t63, 0; t64, 0, -t68 * t74, 0, t65, 0; 0, 0, -t70 * t68, 0, -t77, 0; -t65, 0, -t71 * t77, 0, -t64, 0; t63, 0, t67 * t74, 0, t66, 0; 0, 0, t70 * t67, 0, -t76, 0; -t74, 0, t78, 0, 0, 0; -t71 * t72, 0, -t75, 0, 0, 0; 0, 0, t72, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (40->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t82 = sin(qJ(3));
	t83 = sin(qJ(1));
	t90 = t83 * t82;
	t81 = pkin(9) + qJ(5);
	t79 = sin(t81);
	t84 = cos(qJ(3));
	t89 = t84 * t79;
	t80 = cos(t81);
	t88 = t84 * t80;
	t85 = cos(qJ(1));
	t87 = t85 * t82;
	t86 = t85 * t84;
	t78 = -t83 * t79 + t80 * t87;
	t77 = t79 * t87 + t83 * t80;
	t76 = t85 * t79 + t80 * t90;
	t75 = t79 * t90 - t85 * t80;
	t1 = [t78, 0, t83 * t88, 0, -t75, 0; t76, 0, -t80 * t86, 0, t77, 0; 0, 0, -t82 * t80, 0, -t89, 0; -t86, 0, t90, 0, 0, 0; -t83 * t84, 0, -t87, 0, 0, 0; 0, 0, t84, 0, 0, 0; t77, 0, t83 * t89, 0, t76, 0; t75, 0, -t79 * t86, 0, -t78, 0; 0, 0, -t82 * t79, 0, t88, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end