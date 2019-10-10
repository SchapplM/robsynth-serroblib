% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:59
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
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
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->8), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t47 = sin(qJ(2));
	t48 = sin(qJ(1));
	t55 = t48 * t47;
	t45 = sin(pkin(9));
	t49 = cos(qJ(2));
	t54 = t49 * t45;
	t46 = cos(pkin(9));
	t53 = t49 * t46;
	t50 = cos(qJ(1));
	t52 = t50 * t47;
	t51 = t50 * t49;
	t1 = [t50 * t45 - t48 * t53, -t46 * t52, 0, 0, 0, 0; t48 * t45 + t46 * t51, -t46 * t55, 0, 0, 0, 0; 0, t53, 0, 0, 0, 0; t50 * t46 + t48 * t54, t45 * t52, 0, 0, 0, 0; -t45 * t51 + t48 * t46, t45 * t55, 0, 0, 0, 0; 0, -t54, 0, 0, 0, 0; -t55, t51, 0, 0, 0, 0; t52, t48 * t49, 0, 0, 0, 0; 0, t47, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t68 = sin(qJ(2));
	t69 = sin(qJ(1));
	t76 = t69 * t68;
	t67 = pkin(9) + qJ(4);
	t65 = sin(t67);
	t70 = cos(qJ(2));
	t75 = t70 * t65;
	t66 = cos(t67);
	t74 = t70 * t66;
	t71 = cos(qJ(1));
	t73 = t71 * t68;
	t72 = t71 * t70;
	t64 = t69 * t65 + t66 * t72;
	t63 = -t65 * t72 + t69 * t66;
	t62 = t71 * t65 - t69 * t74;
	t61 = t71 * t66 + t69 * t75;
	t1 = [t62, -t66 * t73, 0, t63, 0, 0; t64, -t66 * t76, 0, -t61, 0, 0; 0, t74, 0, -t68 * t65, 0, 0; t61, t65 * t73, 0, -t64, 0, 0; t63, t65 * t76, 0, t62, 0, 0; 0, -t75, 0, -t68 * t66, 0, 0; -t76, t72, 0, 0, 0, 0; t73, t69 * t70, 0, 0, 0, 0; 0, t68, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t80 = sin(qJ(2));
	t81 = sin(qJ(1));
	t88 = t81 * t80;
	t79 = pkin(9) + qJ(4);
	t77 = sin(t79);
	t82 = cos(qJ(2));
	t87 = t82 * t77;
	t78 = cos(t79);
	t86 = t82 * t78;
	t83 = cos(qJ(1));
	t85 = t83 * t80;
	t84 = t83 * t82;
	t76 = t81 * t77 + t78 * t84;
	t75 = t77 * t84 - t81 * t78;
	t74 = -t83 * t77 + t81 * t86;
	t73 = t83 * t78 + t81 * t87;
	t1 = [-t88, t84, 0, 0, 0, 0; t85, t81 * t82, 0, 0, 0, 0; 0, t80, 0, 0, 0, 0; t74, t78 * t85, 0, t75, 0, 0; -t76, t78 * t88, 0, t73, 0, 0; 0, -t86, 0, t80 * t77, 0, 0; -t73, -t77 * t85, 0, t76, 0, 0; t75, -t77 * t88, 0, t74, 0, 0; 0, t87, 0, t80 * t78, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t79 = sin(qJ(2));
	t80 = sin(qJ(1));
	t87 = t80 * t79;
	t78 = pkin(9) + qJ(4);
	t76 = sin(t78);
	t81 = cos(qJ(2));
	t86 = t81 * t76;
	t77 = cos(t78);
	t85 = t81 * t77;
	t82 = cos(qJ(1));
	t84 = t82 * t79;
	t83 = t82 * t81;
	t75 = t80 * t76 + t77 * t83;
	t74 = t76 * t83 - t80 * t77;
	t73 = -t82 * t76 + t80 * t85;
	t72 = -t82 * t77 - t80 * t86;
	t1 = [-t87, t83, 0, 0, 0, 0; t84, t80 * t81, 0, 0, 0, 0; 0, t79, 0, 0, 0, 0; t72, -t76 * t84, 0, t75, 0, 0; t74, -t76 * t87, 0, t73, 0, 0; 0, t86, 0, t79 * t77, 0, 0; -t73, -t77 * t84, 0, -t74, 0, 0; t75, -t77 * t87, 0, t72, 0, 0; 0, t85, 0, -t79 * t76, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end