% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
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
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
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
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
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
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(3) + pkin(10);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t21 = t18 * t15;
	t16 = cos(t17);
	t19 = cos(qJ(1));
	t20 = t19 * t16;
	t14 = t19 * t15;
	t13 = t18 * t16;
	t1 = [t14, 0, t13, 0, 0, 0; t21, 0, -t20, 0, 0, 0; 0, 0, -t15, 0, 0, 0; t20, 0, -t21, 0, 0, 0; t13, 0, t14, 0, 0, 0; 0, 0, -t16, 0, 0, 0; -t18, 0, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (37->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t74 = sin(qJ(5));
	t75 = sin(qJ(1));
	t81 = t75 * t74;
	t76 = cos(qJ(5));
	t80 = t75 * t76;
	t77 = cos(qJ(1));
	t79 = t77 * t74;
	t78 = t77 * t76;
	t73 = qJ(3) + pkin(10);
	t72 = cos(t73);
	t71 = sin(t73);
	t70 = t71 * t78 - t81;
	t69 = t71 * t79 + t80;
	t68 = t71 * t80 + t79;
	t67 = -t71 * t81 + t78;
	t1 = [t70, 0, t72 * t80, 0, t67, 0; t68, 0, -t72 * t78, 0, t69, 0; 0, 0, -t71 * t76, 0, -t72 * t74, 0; -t69, 0, -t72 * t81, 0, -t68, 0; t67, 0, t72 * t79, 0, t70, 0; 0, 0, t71 * t74, 0, -t72 * t76, 0; -t77 * t72, 0, t75 * t71, 0, 0, 0; -t75 * t72, 0, -t77 * t71, 0, 0, 0; 0, 0, t72, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (83->19), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t94 = qJ(3) + pkin(10);
	t91 = cos(t94);
	t95 = qJ(5) + qJ(6);
	t92 = sin(t95);
	t103 = t91 * t92;
	t93 = cos(t95);
	t102 = t91 * t93;
	t96 = sin(qJ(1));
	t101 = t96 * t92;
	t100 = t96 * t93;
	t97 = cos(qJ(1));
	t99 = t97 * t92;
	t98 = t97 * t93;
	t90 = sin(t94);
	t89 = t90 * t98 - t101;
	t88 = t90 * t99 + t100;
	t87 = t90 * t100 + t99;
	t86 = -t90 * t101 + t98;
	t1 = [t89, 0, t91 * t100, 0, t86, t86; t87, 0, -t91 * t98, 0, t88, t88; 0, 0, -t90 * t93, 0, -t103, -t103; -t88, 0, -t91 * t101, 0, -t87, -t87; t86, 0, t91 * t99, 0, t89, t89; 0, 0, t90 * t92, 0, -t102, -t102; -t97 * t91, 0, t96 * t90, 0, 0, 0; -t96 * t91, 0, -t97 * t90, 0, 0, 0; 0, 0, t91, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end