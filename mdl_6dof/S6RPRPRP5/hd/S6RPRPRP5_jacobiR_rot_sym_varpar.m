% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
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
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
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
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
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
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t77 = pkin(9) + qJ(3);
	t73 = sin(t77);
	t78 = sin(qJ(1));
	t85 = t78 * t73;
	t76 = pkin(10) + qJ(5);
	t74 = cos(t76);
	t84 = t78 * t74;
	t75 = cos(t77);
	t83 = t78 * t75;
	t79 = cos(qJ(1));
	t82 = t79 * t73;
	t81 = t79 * t74;
	t80 = t79 * t75;
	t72 = sin(t76);
	t71 = t78 * t72 + t74 * t80;
	t70 = -t72 * t80 + t84;
	t69 = t79 * t72 - t74 * t83;
	t68 = t72 * t83 + t81;
	t1 = [t69, 0, -t73 * t81, 0, t70, 0; t71, 0, -t73 * t84, 0, -t68, 0; 0, 0, t75 * t74, 0, -t73 * t72, 0; t68, 0, t72 * t82, 0, -t71, 0; t70, 0, t72 * t85, 0, t69, 0; 0, 0, -t75 * t72, 0, -t73 * t74, 0; -t85, 0, t80, 0, 0, 0; t82, 0, t83, 0, 0, 0; 0, 0, t73, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t92 = pkin(9) + qJ(3);
	t88 = sin(t92);
	t93 = sin(qJ(1));
	t100 = t93 * t88;
	t91 = pkin(10) + qJ(5);
	t89 = cos(t91);
	t99 = t93 * t89;
	t90 = cos(t92);
	t98 = t93 * t90;
	t94 = cos(qJ(1));
	t97 = t94 * t88;
	t96 = t94 * t89;
	t95 = t94 * t90;
	t87 = sin(t91);
	t86 = t93 * t87 + t89 * t95;
	t85 = t87 * t95 - t99;
	t84 = -t94 * t87 + t89 * t98;
	t83 = -t87 * t98 - t96;
	t1 = [-t84, 0, -t88 * t96, 0, -t85, 0; t86, 0, -t88 * t99, 0, t83, 0; 0, 0, t90 * t89, 0, -t88 * t87, 0; -t100, 0, t95, 0, 0, 0; t97, 0, t98, 0, 0, 0; 0, 0, t88, 0, 0, 0; t83, 0, -t87 * t97, 0, t86, 0; t85, 0, -t87 * t100, 0, t84, 0; 0, 0, t90 * t87, 0, t88 * t89, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end