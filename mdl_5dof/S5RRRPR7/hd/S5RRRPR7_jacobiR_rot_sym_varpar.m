% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRPR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:18
	% EndTime: 2019-12-29 20:05:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
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
	t1 = [-t14, -t15, 0, 0, 0; t13, -t16, 0, 0, 0; 0, t11, 0, 0, 0; t16, -t13, 0, 0, 0; -t15, -t14, 0, 0, 0; 0, -t9, 0, 0, 0; t12, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t22 = qJ(2) + qJ(3);
	t20 = sin(t22);
	t23 = sin(qJ(1));
	t28 = t23 * t20;
	t21 = cos(t22);
	t27 = t23 * t21;
	t24 = cos(qJ(1));
	t26 = t24 * t20;
	t25 = t24 * t21;
	t1 = [-t27, -t26, -t26, 0, 0; t25, -t28, -t28, 0, 0; 0, t21, t21, 0, 0; t28, -t25, -t25, 0, 0; -t26, -t27, -t27, 0, 0; 0, -t20, -t20, 0, 0; t24, 0, 0, 0, 0; t23, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t87 = qJ(2) + qJ(3);
	t86 = cos(t87);
	t88 = sin(pkin(9));
	t98 = t86 * t88;
	t90 = sin(qJ(1));
	t97 = t90 * t88;
	t89 = cos(pkin(9));
	t96 = t90 * t89;
	t91 = cos(qJ(1));
	t95 = t91 * t88;
	t94 = t91 * t89;
	t85 = sin(t87);
	t93 = t85 * t96;
	t92 = t85 * t94;
	t84 = t91 * t86;
	t83 = t90 * t86;
	t82 = t86 * t89;
	t81 = t85 * t95;
	t80 = t85 * t97;
	t1 = [-t86 * t96 + t95, -t92, -t92, 0, 0; t86 * t94 + t97, -t93, -t93, 0, 0; 0, t82, t82, 0, 0; t86 * t97 + t94, t81, t81, 0, 0; -t86 * t95 + t96, t80, t80, 0, 0; 0, -t98, -t98, 0, 0; -t90 * t85, t84, t84, 0, 0; t91 * t85, t83, t83, 0, 0; 0, t85, t85, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (77->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
	t100 = pkin(9) + qJ(5);
	t96 = sin(t100);
	t101 = qJ(2) + qJ(3);
	t99 = cos(t101);
	t108 = t99 * t96;
	t102 = sin(qJ(1));
	t98 = sin(t101);
	t107 = t102 * t98;
	t94 = t102 * t99;
	t103 = cos(qJ(1));
	t106 = t103 * t98;
	t95 = t103 * t99;
	t97 = cos(t100);
	t105 = t97 * t107;
	t104 = t97 * t106;
	t93 = t99 * t97;
	t92 = t96 * t106;
	t91 = t96 * t107;
	t90 = t102 * t96 + t95 * t97;
	t89 = t102 * t97 - t95 * t96;
	t88 = t103 * t96 - t94 * t97;
	t87 = t103 * t97 + t94 * t96;
	t1 = [t88, -t104, -t104, 0, t89; t90, -t105, -t105, 0, -t87; 0, t93, t93, 0, -t98 * t96; t87, t92, t92, 0, -t90; t89, t91, t91, 0, t88; 0, -t108, -t108, 0, -t98 * t97; -t107, t95, t95, 0, 0; t106, t94, t94, 0, 0; 0, t98, t98, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end