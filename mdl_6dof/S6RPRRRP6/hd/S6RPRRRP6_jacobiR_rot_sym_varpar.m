% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
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
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(10));
	t7 = sin(pkin(10));
	t1 = [-t9 * t8, 0, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t14 = pkin(10) + qJ(3);
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
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t71 = sin(qJ(4));
	t72 = sin(qJ(1));
	t78 = t72 * t71;
	t73 = cos(qJ(4));
	t77 = t72 * t73;
	t74 = cos(qJ(1));
	t76 = t74 * t71;
	t75 = t74 * t73;
	t70 = pkin(10) + qJ(3);
	t69 = cos(t70);
	t68 = sin(t70);
	t67 = t69 * t75 + t78;
	t66 = -t69 * t76 + t77;
	t65 = -t69 * t77 + t76;
	t64 = t69 * t78 + t75;
	t1 = [t65, 0, -t68 * t75, t66, 0, 0; t67, 0, -t68 * t77, -t64, 0, 0; 0, 0, t69 * t73, -t68 * t71, 0, 0; t64, 0, t68 * t76, -t67, 0, 0; t66, 0, t68 * t78, t65, 0, 0; 0, 0, -t69 * t71, -t68 * t73, 0, 0; -t72 * t68, 0, t74 * t69, 0, 0, 0; t74 * t68, 0, t72 * t69, 0, 0, 0; 0, 0, t68, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (81->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t92 = pkin(10) + qJ(3);
	t88 = sin(t92);
	t93 = qJ(4) + qJ(5);
	t90 = sin(t93);
	t101 = t88 * t90;
	t91 = cos(t93);
	t100 = t88 * t91;
	t94 = sin(qJ(1));
	t99 = t94 * t90;
	t98 = t94 * t91;
	t95 = cos(qJ(1));
	t97 = t95 * t90;
	t96 = t95 * t91;
	t89 = cos(t92);
	t87 = t89 * t96 + t99;
	t86 = -t89 * t97 + t98;
	t85 = -t89 * t98 + t97;
	t84 = t89 * t99 + t96;
	t1 = [t85, 0, -t88 * t96, t86, t86, 0; t87, 0, -t88 * t98, -t84, -t84, 0; 0, 0, t89 * t91, -t101, -t101, 0; t84, 0, t88 * t97, -t87, -t87, 0; t86, 0, t88 * t99, t85, t85, 0; 0, 0, -t89 * t90, -t100, -t100, 0; -t94 * t88, 0, t95 * t89, 0, 0, 0; t95 * t88, 0, t94 * t89, 0, 0, 0; 0, 0, t88, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (81->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t93 = pkin(10) + qJ(3);
	t89 = sin(t93);
	t94 = qJ(4) + qJ(5);
	t91 = sin(t94);
	t102 = t89 * t91;
	t92 = cos(t94);
	t101 = t89 * t92;
	t95 = sin(qJ(1));
	t100 = t95 * t91;
	t99 = t95 * t92;
	t96 = cos(qJ(1));
	t98 = t96 * t91;
	t97 = t96 * t92;
	t90 = cos(t93);
	t88 = t90 * t97 + t100;
	t87 = -t90 * t98 + t99;
	t86 = -t90 * t99 + t98;
	t85 = t90 * t100 + t97;
	t1 = [t86, 0, -t89 * t97, t87, t87, 0; t88, 0, -t89 * t99, -t85, -t85, 0; 0, 0, t90 * t92, -t102, -t102, 0; t85, 0, t89 * t98, -t88, -t88, 0; t87, 0, t89 * t100, t86, t86, 0; 0, 0, -t90 * t91, -t101, -t101, 0; -t95 * t89, 0, t96 * t90, 0, 0, 0; t96 * t89, 0, t95 * t90, 0, 0, 0; 0, 0, t89, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end