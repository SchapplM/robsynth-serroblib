% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:47
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
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
	% StartTime: 2019-10-10 11:47:48
	% EndTime: 2019-10-10 11:47:49
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
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
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
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
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
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (50->24), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->23)
	t93 = sin(qJ(5));
	t94 = sin(qJ(3));
	t97 = cos(qJ(5));
	t98 = cos(qJ(3));
	t102 = t93 * t94 + t97 * t98;
	t95 = sin(qJ(2));
	t109 = t102 * t95;
	t100 = cos(qJ(1));
	t96 = sin(qJ(1));
	t99 = cos(qJ(2));
	t107 = t96 * t99;
	t87 = t100 * t98 + t94 * t107;
	t88 = -t100 * t94 + t98 * t107;
	t105 = -t87 * t97 + t88 * t93;
	t106 = t100 * t99;
	t104 = t87 * t93 + t88 * t97;
	t89 = t94 * t106 - t96 * t98;
	t90 = t98 * t106 + t96 * t94;
	t82 = t89 * t97 - t90 * t93;
	t83 = t89 * t93 + t90 * t97;
	t103 = t93 * t98 - t94 * t97;
	t85 = t103 * t95;
	t1 = [-t104, -t100 * t109, -t82, 0, t82, 0; t83, -t96 * t109, t105, 0, -t105, 0; 0, t102 * t99, t85, 0, -t85, 0; t105, t100 * t85, t83, 0, -t83, 0; t82, t96 * t85, t104, 0, -t104, 0; 0, -t103 * t99, t109, 0, -t109, 0; t96 * t95, -t106, 0, 0, 0, 0; -t100 * t95, -t107, 0, 0, 0, 0; 0, -t95, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:47:49
	% EndTime: 2019-10-10 11:47:49
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (50->25), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->23)
	t103 = sin(qJ(5));
	t107 = cos(qJ(5));
	t104 = sin(qJ(3));
	t108 = cos(qJ(3));
	t110 = cos(qJ(1));
	t106 = sin(qJ(1));
	t109 = cos(qJ(2));
	t118 = t106 * t109;
	t97 = t104 * t118 + t110 * t108;
	t98 = -t110 * t104 + t108 * t118;
	t116 = t98 * t103 - t97 * t107;
	t105 = sin(qJ(2));
	t112 = t103 * t104 + t107 * t108;
	t120 = t112 * t105;
	t117 = t110 * t109;
	t115 = t97 * t103 + t98 * t107;
	t100 = t106 * t104 + t108 * t117;
	t99 = t104 * t117 - t106 * t108;
	t93 = t100 * t107 + t99 * t103;
	t114 = t100 * t103 - t99 * t107;
	t113 = t103 * t108 - t104 * t107;
	t95 = t113 * t105;
	t1 = [-t115, -t110 * t120, t114, 0, -t114, 0; t93, -t106 * t120, t116, 0, -t116, 0; 0, t112 * t109, t95, 0, -t95, 0; t116, t110 * t95, t93, 0, -t93, 0; -t114, t106 * t95, t115, 0, -t115, 0; 0, -t113 * t109, t120, 0, -t120, 0; t106 * t105, -t117, 0, 0, 0, 0; -t110 * t105, -t118, 0, 0, 0, 0; 0, -t105, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end