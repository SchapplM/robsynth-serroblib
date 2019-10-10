% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
	t82 = qJ(2) + pkin(10);
	t78 = sin(t82);
	t83 = sin(qJ(1));
	t90 = t83 * t78;
	t81 = qJ(4) + pkin(11);
	t79 = cos(t81);
	t89 = t83 * t79;
	t80 = cos(t82);
	t88 = t83 * t80;
	t84 = cos(qJ(1));
	t87 = t84 * t78;
	t86 = t84 * t79;
	t85 = t84 * t80;
	t77 = sin(t81);
	t76 = t83 * t77 + t79 * t85;
	t75 = -t77 * t85 + t89;
	t74 = t84 * t77 - t79 * t88;
	t73 = t77 * t88 + t86;
	t1 = [t74, -t78 * t86, 0, t75, 0, 0; t76, -t78 * t89, 0, -t73, 0, 0; 0, t80 * t79, 0, -t78 * t77, 0, 0; t73, t77 * t87, 0, -t76, 0, 0; t75, t77 * t90, 0, t74, 0, 0; 0, -t80 * t77, 0, -t78 * t79, 0, 0; -t90, t85, 0, 0, 0, 0; t87, t88, 0, 0, 0, 0; 0, t78, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (115->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t95 = qJ(4) + pkin(11) + qJ(6);
	t91 = sin(t95);
	t96 = qJ(2) + pkin(10);
	t93 = sin(t96);
	t104 = t93 * t91;
	t92 = cos(t95);
	t103 = t93 * t92;
	t97 = sin(qJ(1));
	t102 = t97 * t93;
	t94 = cos(t96);
	t101 = t97 * t94;
	t98 = cos(qJ(1));
	t100 = t98 * t93;
	t99 = t98 * t94;
	t90 = t91 * t97 + t92 * t99;
	t89 = -t91 * t99 + t92 * t97;
	t88 = -t101 * t92 + t91 * t98;
	t87 = t101 * t91 + t92 * t98;
	t1 = [t88, -t92 * t100, 0, t89, 0, t89; t90, -t92 * t102, 0, -t87, 0, -t87; 0, t94 * t92, 0, -t104, 0, -t104; t87, t91 * t100, 0, -t90, 0, -t90; t89, t91 * t102, 0, t88, 0, t88; 0, -t94 * t91, 0, -t103, 0, -t103; -t102, t99, 0, 0, 0, 0; t100, t101, 0, 0, 0, 0; 0, t93, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end