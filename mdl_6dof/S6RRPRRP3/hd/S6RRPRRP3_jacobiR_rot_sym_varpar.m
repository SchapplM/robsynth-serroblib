% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(2) + pkin(10);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t23 = t18 * t15;
	t16 = cos(t17);
	t22 = t18 * t16;
	t19 = cos(qJ(1));
	t21 = t19 * t15;
	t20 = t19 * t16;
	t1 = [-t22, -t21, 0, 0, 0, 0; t20, -t23, 0, 0, 0, 0; 0, t16, 0, 0, 0, 0; t23, -t20, 0, 0, 0, 0; -t21, -t22, 0, 0, 0, 0; 0, -t15, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
	t74 = sin(qJ(4));
	t75 = sin(qJ(1));
	t81 = t75 * t74;
	t76 = cos(qJ(4));
	t80 = t75 * t76;
	t77 = cos(qJ(1));
	t79 = t77 * t74;
	t78 = t77 * t76;
	t73 = qJ(2) + pkin(10);
	t72 = cos(t73);
	t71 = sin(t73);
	t70 = t72 * t78 + t81;
	t69 = -t72 * t79 + t80;
	t68 = -t72 * t80 + t79;
	t67 = t72 * t81 + t78;
	t1 = [t68, -t71 * t78, 0, t69, 0, 0; t70, -t71 * t80, 0, -t67, 0, 0; 0, t72 * t76, 0, -t71 * t74, 0, 0; t67, t71 * t79, 0, -t70, 0, 0; t69, t71 * t81, 0, t68, 0, 0; 0, -t72 * t74, 0, -t71 * t76, 0, 0; -t75 * t71, t77 * t72, 0, 0, 0, 0; t77 * t71, t75 * t72, 0, 0, 0, 0; 0, t71, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (81->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t95 = qJ(2) + pkin(10);
	t91 = sin(t95);
	t96 = qJ(4) + qJ(5);
	t93 = sin(t96);
	t104 = t91 * t93;
	t94 = cos(t96);
	t103 = t91 * t94;
	t97 = sin(qJ(1));
	t102 = t97 * t93;
	t101 = t97 * t94;
	t98 = cos(qJ(1));
	t100 = t98 * t93;
	t99 = t98 * t94;
	t92 = cos(t95);
	t90 = t92 * t99 + t102;
	t89 = -t92 * t100 + t101;
	t88 = -t92 * t101 + t100;
	t87 = t92 * t102 + t99;
	t1 = [t88, -t91 * t99, 0, t89, t89, 0; t90, -t91 * t101, 0, -t87, -t87, 0; 0, t92 * t94, 0, -t104, -t104, 0; t87, t91 * t100, 0, -t90, -t90, 0; t89, t91 * t102, 0, t88, t88, 0; 0, -t92 * t93, 0, -t103, -t103, 0; -t97 * t91, t98 * t92, 0, 0, 0, 0; t98 * t91, t97 * t92, 0, 0, 0, 0; 0, t91, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (81->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t96 = qJ(2) + pkin(10);
	t92 = sin(t96);
	t97 = qJ(4) + qJ(5);
	t94 = sin(t97);
	t105 = t92 * t94;
	t95 = cos(t97);
	t104 = t92 * t95;
	t98 = sin(qJ(1));
	t103 = t98 * t94;
	t102 = t98 * t95;
	t99 = cos(qJ(1));
	t101 = t99 * t94;
	t100 = t99 * t95;
	t93 = cos(t96);
	t91 = t93 * t100 + t103;
	t90 = -t93 * t101 + t102;
	t89 = -t93 * t102 + t101;
	t88 = t93 * t103 + t100;
	t1 = [t89, -t92 * t100, 0, t90, t90, 0; t91, -t92 * t102, 0, -t88, -t88, 0; 0, t93 * t95, 0, -t105, -t105, 0; t88, t92 * t101, 0, -t91, -t91, 0; t90, t92 * t103, 0, t89, t89, 0; 0, -t93 * t94, 0, -t104, -t104, 0; -t98 * t92, t99 * t93, 0, 0, 0, 0; t99 * t92, t98 * t93, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end