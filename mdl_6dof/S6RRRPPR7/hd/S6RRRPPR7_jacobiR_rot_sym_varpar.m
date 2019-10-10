% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
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
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
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
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t60 = sin(qJ(2));
	t61 = sin(qJ(1));
	t71 = t61 * t60;
	t62 = cos(qJ(3));
	t70 = t61 * t62;
	t59 = sin(qJ(3));
	t63 = cos(qJ(2));
	t69 = t63 * t59;
	t68 = t63 * t62;
	t64 = cos(qJ(1));
	t67 = t64 * t60;
	t66 = t64 * t62;
	t65 = t64 * t63;
	t58 = t61 * t59 + t62 * t65;
	t57 = -t59 * t65 + t70;
	t56 = t64 * t59 - t61 * t68;
	t55 = t61 * t69 + t66;
	t1 = [t56, -t60 * t66, t57, 0, 0, 0; t58, -t60 * t70, -t55, 0, 0, 0; 0, t68, -t60 * t59, 0, 0, 0; t55, t59 * t67, -t58, 0, 0, 0; t57, t59 * t71, t56, 0, 0, 0; 0, -t69, -t60 * t62, 0, 0, 0; -t71, t65, 0, 0, 0, 0; t67, t61 * t63, 0, 0, 0, 0; 0, t60, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:46
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t76 = sin(qJ(2));
	t77 = sin(qJ(1));
	t87 = t77 * t76;
	t78 = cos(qJ(3));
	t86 = t77 * t78;
	t75 = sin(qJ(3));
	t79 = cos(qJ(2));
	t85 = t79 * t75;
	t84 = t79 * t78;
	t80 = cos(qJ(1));
	t83 = t80 * t76;
	t82 = t80 * t78;
	t81 = t80 * t79;
	t74 = t77 * t75 + t78 * t81;
	t73 = t75 * t81 - t86;
	t72 = -t80 * t75 + t77 * t84;
	t71 = -t77 * t85 - t82;
	t1 = [-t72, -t76 * t82, -t73, 0, 0, 0; t74, -t76 * t86, t71, 0, 0, 0; 0, t84, -t76 * t75, 0, 0, 0; -t87, t81, 0, 0, 0, 0; t83, t77 * t79, 0, 0, 0, 0; 0, t76, 0, 0, 0, 0; t71, -t75 * t83, t74, 0, 0, 0; t73, -t75 * t87, t72, 0, 0, 0; 0, t85, t76 * t78, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:46
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->19), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
	t85 = sin(qJ(2));
	t82 = sin(pkin(10));
	t83 = cos(pkin(10));
	t84 = sin(qJ(3));
	t87 = cos(qJ(3));
	t92 = t82 * t84 + t83 * t87;
	t91 = t92 * t85;
	t86 = sin(qJ(1));
	t88 = cos(qJ(2));
	t99 = t86 * t88;
	t89 = cos(qJ(1));
	t98 = t89 * t88;
	t76 = -t84 * t99 - t89 * t87;
	t77 = -t89 * t84 + t87 * t99;
	t97 = t76 * t83 + t77 * t82;
	t78 = t84 * t98 - t86 * t87;
	t79 = t86 * t84 + t87 * t98;
	t96 = t78 * t82 + t79 * t83;
	t95 = t76 * t82 - t77 * t83;
	t94 = t78 * t83 - t79 * t82;
	t93 = t82 * t87 - t83 * t84;
	t90 = t93 * t85;
	t1 = [t95, -t89 * t91, -t94, 0, 0, 0; t96, -t86 * t91, t97, 0, 0, 0; 0, t92 * t88, t90, 0, 0, 0; t97, t89 * t90, t96, 0, 0, 0; t94, t86 * t90, -t95, 0, 0, 0; 0, -t93 * t88, t91, 0, 0, 0; t86 * t85, -t98, 0, 0, 0, 0; -t89 * t85, -t99, 0, 0, 0, 0; 0, -t85, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:46
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (94->27), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->24)
	t101 = pkin(10) + qJ(6);
	t100 = cos(t101);
	t102 = sin(qJ(3));
	t105 = cos(qJ(3));
	t107 = cos(qJ(1));
	t104 = sin(qJ(1));
	t106 = cos(qJ(2));
	t114 = t104 * t106;
	t93 = t102 * t114 + t107 * t105;
	t94 = -t107 * t102 + t105 * t114;
	t99 = sin(t101);
	t117 = t93 * t100 - t94 * t99;
	t103 = sin(qJ(2));
	t110 = t100 * t105 + t102 * t99;
	t116 = t110 * t103;
	t113 = t107 * t106;
	t111 = t94 * t100 + t93 * t99;
	t95 = t102 * t113 - t104 * t105;
	t96 = t104 * t102 + t105 * t113;
	t89 = t96 * t100 + t95 * t99;
	t88 = t95 * t100 - t96 * t99;
	t109 = t100 * t102 - t105 * t99;
	t90 = t109 * t103;
	t1 = [-t111, -t107 * t116, -t88, 0, 0, t88; t89, -t104 * t116, -t117, 0, 0, t117; 0, t110 * t106, -t90, 0, 0, t90; -t117, -t107 * t90, t89, 0, 0, -t89; t88, -t104 * t90, t111, 0, 0, -t111; 0, t109 * t106, t116, 0, 0, -t116; t104 * t103, -t113, 0, 0, 0, 0; -t107 * t103, -t114, 0, 0, 0, 0; 0, -t103, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end