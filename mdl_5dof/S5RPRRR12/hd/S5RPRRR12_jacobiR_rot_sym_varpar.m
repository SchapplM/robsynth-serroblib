% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:13:40
	% EndTime: 2019-12-31 19:13:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:13:40
	% EndTime: 2019-12-31 19:13:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:13:40
	% EndTime: 2019-12-31 19:13:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t5, 0, 0, 0, 0; -t6, 0, 0, 0, 0; 0, 0, 0, 0, 0; t6, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:13:40
	% EndTime: 2019-12-31 19:13:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0; t16, 0, -t15, 0, 0; 0, 0, -t11, 0, 0; t15, 0, -t16, 0, 0; t9, 0, t10, 0, 0; 0, 0, -t13, 0, 0; -t12, 0, 0, 0, 0; t14, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:13:40
	% EndTime: 2019-12-31 19:13:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t24 = qJ(3) + qJ(4);
	t22 = sin(t24);
	t25 = sin(qJ(1));
	t28 = t25 * t22;
	t23 = cos(t24);
	t26 = cos(qJ(1));
	t27 = t26 * t23;
	t21 = t26 * t22;
	t20 = t25 * t23;
	t1 = [t21, 0, t20, t20, 0; t28, 0, -t27, -t27, 0; 0, 0, -t22, -t22, 0; t27, 0, -t28, -t28, 0; t20, 0, t21, t21, 0; 0, 0, -t23, -t23, 0; -t25, 0, 0, 0, 0; t26, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:13:41
	% EndTime: 2019-12-31 19:13:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (50->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t94 = qJ(3) + qJ(4);
	t92 = sin(t94);
	t97 = cos(qJ(5));
	t106 = t92 * t97;
	t95 = sin(qJ(5));
	t96 = sin(qJ(1));
	t105 = t96 * t95;
	t104 = t96 * t97;
	t98 = cos(qJ(1));
	t103 = t98 * t92;
	t102 = t98 * t95;
	t101 = t98 * t97;
	t93 = cos(t94);
	t100 = t93 * t105;
	t99 = t93 * t101;
	t91 = t96 * t92;
	t90 = t92 * t95;
	t89 = t93 * t102;
	t88 = t93 * t104;
	t87 = t92 * t101 - t105;
	t86 = t92 * t102 + t104;
	t85 = t92 * t104 + t102;
	t84 = -t92 * t105 + t101;
	t1 = [t87, 0, t88, t88, t84; t85, 0, -t99, -t99, t86; 0, 0, -t106, -t106, -t93 * t95; -t86, 0, -t100, -t100, -t85; t84, 0, t89, t89, t87; 0, 0, t90, t90, -t93 * t97; -t98 * t93, 0, t91, t91, 0; -t96 * t93, 0, -t103, -t103, 0; 0, 0, t93, t93, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end