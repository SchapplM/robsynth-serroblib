% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
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
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(9));
	t7 = sin(pkin(9));
	t1 = [-t9 * t8, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t14 = pkin(9) + qJ(3);
	t12 = sin(t14);
	t15 = sin(qJ(1));
	t20 = t15 * t12;
	t13 = cos(t14);
	t19 = t15 * t13;
	t16 = cos(qJ(1));
	t18 = t16 * t12;
	t17 = t16 * t13;
	t1 = [-t19, 0, -t18, 0, 0, 0; t17, 0, -t20, 0, 0, 0; 0, 0, t13, 0, 0, 0; t20, 0, -t17, 0, 0, 0; -t18, 0, -t19, 0, 0, 0; 0, 0, -t12, 0, 0, 0; t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t54 = sin(pkin(10));
	t56 = sin(qJ(1));
	t61 = t56 * t54;
	t55 = cos(pkin(10));
	t60 = t56 * t55;
	t57 = cos(qJ(1));
	t59 = t57 * t54;
	t58 = t57 * t55;
	t53 = pkin(9) + qJ(3);
	t52 = cos(t53);
	t51 = sin(t53);
	t1 = [-t52 * t60 + t59, 0, -t51 * t58, 0, 0, 0; t52 * t58 + t61, 0, -t51 * t60, 0, 0, 0; 0, 0, t52 * t55, 0, 0, 0; t52 * t61 + t58, 0, t51 * t59, 0, 0, 0; -t52 * t59 + t60, 0, t51 * t61, 0, 0, 0; 0, 0, -t52 * t54, 0, 0, 0; -t56 * t51, 0, t57 * t52, 0, 0, 0; t57 * t51, 0, t56 * t52, 0, 0, 0; 0, 0, t51, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (24->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t70 = sin(pkin(10));
	t72 = sin(qJ(1));
	t77 = t72 * t70;
	t71 = cos(pkin(10));
	t76 = t72 * t71;
	t73 = cos(qJ(1));
	t75 = t73 * t70;
	t74 = t73 * t71;
	t69 = pkin(9) + qJ(3);
	t68 = cos(t69);
	t67 = sin(t69);
	t1 = [-t68 * t76 + t75, 0, -t67 * t74, 0, 0, 0; t68 * t74 + t77, 0, -t67 * t76, 0, 0, 0; 0, 0, t68 * t71, 0, 0, 0; -t72 * t67, 0, t73 * t68, 0, 0, 0; t73 * t67, 0, t72 * t68, 0, 0, 0; 0, 0, t67, 0, 0, 0; -t68 * t77 - t74, 0, -t67 * t75, 0, 0, 0; t68 * t75 - t76, 0, -t67 * t77, 0, 0, 0; 0, 0, t68 * t70, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (73->22), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->26)
	t87 = pkin(9) + qJ(3);
	t85 = sin(t87);
	t88 = sin(pkin(10));
	t89 = cos(pkin(10));
	t90 = sin(qJ(6));
	t92 = cos(qJ(6));
	t97 = t88 * t92 - t89 * t90;
	t95 = t97 * t85;
	t91 = sin(qJ(1));
	t103 = t91 * t88;
	t102 = t91 * t89;
	t93 = cos(qJ(1));
	t101 = t93 * t88;
	t100 = t93 * t89;
	t86 = cos(t87);
	t80 = t86 * t103 + t100;
	t81 = t86 * t102 - t101;
	t99 = t80 * t92 - t81 * t90;
	t98 = -t80 * t90 - t81 * t92;
	t96 = t88 * t90 + t89 * t92;
	t94 = t96 * t85;
	t83 = t86 * t100 + t103;
	t82 = t86 * t101 - t102;
	t79 = t82 * t90 + t83 * t92;
	t78 = t82 * t92 - t83 * t90;
	t1 = [t98, 0, -t93 * t94, 0, 0, t78; t79, 0, -t91 * t94, 0, 0, t99; 0, 0, t96 * t86, 0, 0, t95; -t99, 0, -t93 * t95, 0, 0, -t79; t78, 0, -t91 * t95, 0, 0, t98; 0, 0, t97 * t86, 0, 0, -t94; t91 * t85, 0, -t93 * t86, 0, 0, 0; -t93 * t85, 0, -t91 * t86, 0, 0, 0; 0, 0, -t85, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end