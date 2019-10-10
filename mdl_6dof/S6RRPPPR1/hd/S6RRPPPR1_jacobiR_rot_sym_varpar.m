% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
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
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
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
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(2) + pkin(9);
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
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (23->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t57 = sin(pkin(10));
	t59 = sin(qJ(1));
	t64 = t59 * t57;
	t58 = cos(pkin(10));
	t63 = t59 * t58;
	t60 = cos(qJ(1));
	t62 = t60 * t57;
	t61 = t60 * t58;
	t56 = qJ(2) + pkin(9);
	t55 = cos(t56);
	t54 = sin(t56);
	t1 = [-t55 * t63 + t62, -t54 * t61, 0, 0, 0, 0; t55 * t61 + t64, -t54 * t63, 0, 0, 0, 0; 0, t55 * t58, 0, 0, 0, 0; t55 * t64 + t61, t54 * t62, 0, 0, 0, 0; -t55 * t62 + t63, t54 * t64, 0, 0, 0, 0; 0, -t55 * t57, 0, 0, 0, 0; -t59 * t54, t60 * t55, 0, 0, 0, 0; t60 * t54, t59 * t55, 0, 0, 0, 0; 0, t54, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (24->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t73 = sin(pkin(10));
	t75 = sin(qJ(1));
	t80 = t75 * t73;
	t74 = cos(pkin(10));
	t79 = t75 * t74;
	t76 = cos(qJ(1));
	t78 = t76 * t73;
	t77 = t76 * t74;
	t72 = qJ(2) + pkin(9);
	t71 = cos(t72);
	t70 = sin(t72);
	t1 = [-t71 * t79 + t78, -t70 * t77, 0, 0, 0, 0; t71 * t77 + t80, -t70 * t79, 0, 0, 0, 0; 0, t71 * t74, 0, 0, 0, 0; -t75 * t70, t76 * t71, 0, 0, 0, 0; t76 * t70, t75 * t71, 0, 0, 0, 0; 0, t70, 0, 0, 0, 0; -t71 * t80 - t77, -t70 * t78, 0, 0, 0, 0; t71 * t78 - t79, -t70 * t80, 0, 0, 0, 0; 0, t71 * t73, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (73->22), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->26)
	t91 = sin(pkin(10));
	t92 = cos(pkin(10));
	t93 = sin(qJ(6));
	t95 = cos(qJ(6));
	t100 = t91 * t95 - t92 * t93;
	t90 = qJ(2) + pkin(9);
	t88 = sin(t90);
	t98 = t100 * t88;
	t94 = sin(qJ(1));
	t106 = t94 * t91;
	t105 = t94 * t92;
	t96 = cos(qJ(1));
	t104 = t96 * t91;
	t103 = t96 * t92;
	t89 = cos(t90);
	t83 = t89 * t106 + t103;
	t84 = t89 * t105 - t104;
	t102 = t83 * t95 - t84 * t93;
	t101 = -t83 * t93 - t84 * t95;
	t99 = t91 * t93 + t92 * t95;
	t97 = t99 * t88;
	t86 = t89 * t103 + t106;
	t85 = t89 * t104 - t105;
	t82 = t85 * t93 + t86 * t95;
	t81 = t85 * t95 - t86 * t93;
	t1 = [t101, -t96 * t97, 0, 0, 0, t81; t82, -t94 * t97, 0, 0, 0, t102; 0, t99 * t89, 0, 0, 0, t98; -t102, -t96 * t98, 0, 0, 0, -t82; t81, -t94 * t98, 0, 0, 0, t101; 0, t100 * t89, 0, 0, 0, -t97; t94 * t88, -t96 * t89, 0, 0, 0, 0; -t96 * t88, -t94 * t89, 0, 0, 0, 0; 0, -t88, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end